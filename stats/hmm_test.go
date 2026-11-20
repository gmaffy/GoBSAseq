package stats

import (
	"io"
	"math"
	"math/rand"
	"sort"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"gonum.org/v1/gonum/stat"
)

// f2Model is the worked example: F2 bulk of 20 plants (U = 40), no
// sequencing error, so results match the hand calculation.
func f2Model(t *testing.T) *poolModel {
	m, err := newPoolModel("F2", 20)
	if err != nil {
		t.Fatal(err)
	}
	if m.U != 40 {
		t.Fatalf("U = %d, want 40", m.U)
	}
	m.err = 0
	return &m
}

// Single SNP, 28 high-allele reads of 40: LOD 0.87 (vs 1.43 for the
// read-level LOD that ignores pool sampling).
func TestHMMLODSingleSNP(t *testing.T) {
	m := f2Model(t)
	L := m.chromLikelihoods(hmmBulkData{h: []int{28}, d: []int{40}}, markerBins(1), []float64{0})
	if got := m.lodOneBulk(L[0]); math.Abs(got-0.87) > 0.005 {
		t.Fatalf("LOD = %.4f, want 0.87", got)
	}
}

// Three linked SNPs, per-chromosome switch probability 0.008 between them:
// LOD 1.22 at the middle SNP (hand calculation by explicit forward sums).
func TestHMMLODThreeLinkedSNPs(t *testing.T) {
	m := f2Model(t)
	dM := -0.5 * math.Log(1-2*0.008) // Haldane: r = 0.008
	L := m.chromLikelihoods(hmmBulkData{h: []int{27, 28, 26}, d: []int{40, 40, 38}}, markerBins(3), []float64{0, dM, dM})
	if got := m.lodOneBulk(L[1]); math.Abs(got-1.22) > 0.005 {
		t.Fatalf("LOD = %.4f, want 1.22", got)
	}
}

// Two-bulk differential LOD. When many linked SNPs pin each pool's frequency,
// both hypotheses fit the pools, and the LOD then reflects only how far apart
// the bulks are: segregation distortion (both bulks at 0.7) scores like no
// shift (both at 0.5), while opposite shifts (0.7 vs 0.3) score far higher.
// (The no-shift LOD is not 0: pinning a pool exactly beats drawing it from
// Binomial(U, f); the genome-wide null simulation absorbs that offset.)
func TestHMMTwoBulkDifferential(t *testing.T) {
	m := f2Model(t)
	two := newTwoBulkLOD(m, m)
	n := 20
	dM := make([]float64, n)
	for i := 1; i < n; i++ {
		dM[i] = 0.002
	}
	reads := func(f float64) hmmBulkData {
		d := hmmBulkData{h: make([]int, n), d: make([]int, n)}
		for i := range n {
			d.d[i], d.h[i] = 50, int(math.Round(50*f))
		}
		return d
	}
	lod := func(fh, fl float64) float64 {
		return two.lod(m.chromLikelihoods(reads(fh), markerBins(n), dM)[10], m.chromLikelihoods(reads(fl), markerBins(n), dM)[10])
	}
	none, dist, qtl := lod(0.5, 0.5), lod(0.7, 0.7), lod(0.7, 0.3)
	t.Logf("LOD no shift %.2f, distortion %.2f, QTL %.2f", none, dist, qtl)
	if math.Abs(dist-none) > 0.3 || qtl-none < 2 {
		t.Fatalf("no shift %.2f, distortion %.2f (want ≈ no shift), QTL %.2f (want ≥ no shift + 2)", none, dist, qtl)
	}
}

// Population model: mean and variance of the pool frequency match the
// population's p0 and p0(1−p0)/(c·n).
func TestPoolModelMoments(t *testing.T) {
	for _, pop := range []string{"F2", "F3", "RIL", "BC1H", "BC2L", "BC1F2H"} {
		m, err := newPoolModel(pop, 30)
		if err != nil {
			t.Fatal(err)
		}
		mean, sq := 0.0, 0.0
		for j, p := range m.prior {
			mean += p * m.freq[j]
			sq += p * m.freq[j] * m.freq[j]
		}
		p0, _ := ExpectedAF(pop)
		want := p0 * (1 - p0) / (PopulationVarianceScale(pop) * 30)
		if math.Abs(mean-p0) > 1e-9 || math.Abs((sq-mean*mean)/want-1) > 0.05 {
			t.Errorf("%s: mean %.4f (want %.4f), var %.6f (want %.6f)", pop, mean, p0, sq-mean*mean, want)
		}
	}
}

// The genome-wide null runs end to end on a small genome and gives finite,
// ordered thresholds.
func TestGenomeNullThresholdsRuns(t *testing.T) {
	color.Output = io.Discard
	var raw []BSAstats
	for _, c := range []string{"c1", "c2"} {
		for i := range 1500 {
			raw = append(raw, BSAstats{CHROM: c, POS: int64(1000 + 10_000*i), HighDepth: 40, LowDepth: 40, Depth: 40})
		}
	}
	cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high", NullReps: 20}
	th, err := genomeNullThresholds(cfg, "2p2b", raw)
	zt, ht := th.Z, th.HMM
	if err != nil {
		t.Fatal(err)
	}
	if !(ht.LODP99 >= ht.LODP95 && ht.LODP95 > 0 && zt.CompositeZP99 >= zt.CompositeZP95) || math.IsInf(ht.LODP99, 0) {
		t.Fatalf("thresholds HMM %+v, Z %+v", ht, zt)
	}
	if tb := th.TwoBulk; !(tb.GstatP99 >= tb.GstatP95 && tb.GstatP95 > 0 && tb.DeltaSIP99 > 0 && tb.DeltaSIMp99 < 0 && tb.ED4P99 > 0) {
		t.Fatalf("per-statistic thresholds %+v", tb)
	}
	if zt.NumStats != 2 {
		t.Fatalf("NumStats = %d, want 2", zt.NumStats)
	}
}

// The HMM uses only each position's representative view: a duplicated view
// at a position does not change the LOD.
func TestHMMIgnoresDuplicateViews(t *testing.T) {
	cfg := utils.AnalysisConfig{Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high"}
	setup, err := newHMMSetup(cfg, "2p2b")
	if err != nil {
		t.Fatal(err)
	}
	var raw []BSAstats
	for i := range 200 {
		h := 30
		if i >= 90 && i < 110 {
			h = 40 // a shift in the high bulk
		}
		raw = append(raw, BSAstats{CHROM: "c", POS: int64(1 + 20_000*i), HighBulkH: h, HighBulkL: 50 - h, LowBulkH: 25, LowBulkL: 25, HighDepth: 50, LowDepth: 50})
	}
	single := append([]BSAstats(nil), raw...)
	computeHMMLOD(setup, single, newMarkerOrder(single))
	dup := append(append([]BSAstats(nil), raw...), raw[100])     // second view at position 100, same depth
	dup[len(dup)-1].HighBulkH, dup[len(dup)-1].HighBulkL = 50, 0 // a different contrast
	computeHMMLOD(setup, dup, newMarkerOrder(dup))
	for i := range single {
		if single[i].HMMLOD != dup[i].HMMLOD {
			t.Fatalf("marker %d: LOD %v without the duplicate view, %v with it", i, single[i].HMMLOD, dup[i].HMMLOD)
		}
	}
}

// extremeQuantile recovers the 99% quantile of a Gumbel (light) and a
// Fréchet-type GEV (heavy tail, shape 0.2) from 200 maxima with lower error
// than the empirical quantile.
func TestExtremeQuantile(t *testing.T) {
	const mu, beta = 3.0, 0.5
	for _, shape := range []float64{0, 0.2} {
		draw := func(u float64) float64 { // GEV(μ, β, ξ) inverse CDF
			if shape == 0 {
				return mu - beta*math.Log(-math.Log(u))
			}
			return mu + beta*(math.Pow(-math.Log(u), -shape)-1)/shape
		}
		want := draw(0.99)
		rng := rand.New(rand.NewSource(3))
		var errG, errE float64
		for range 300 {
			x := make([]float64, 200)
			for i := range x {
				x[i] = draw(rng.Float64())
			}
			g := extremeQuantile(x, 0.99)
			sort.Float64s(x)
			e := stat.Quantile(0.99, stat.Empirical, x, nil)
			errG += (g - want) * (g - want)
			errE += (e - want) * (e - want)
		}
		if !(errG < errE) || math.Sqrt(errG/300) > 0.15*want {
			t.Errorf("shape %.1f: RMSE GEV %.3f, empirical %.3f (true q99 %.2f)", shape, math.Sqrt(errG/300), math.Sqrt(errE/300), want)
		}
	}
}
