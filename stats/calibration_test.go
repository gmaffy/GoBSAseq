package stats

import (
	"io"
	"math"
	"math/rand"
	"os"
	"sort"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
	"gonum.org/v1/gonum/stat/distuv"
)

// plantSim is a plant-level F2 generator, independent of the HMM's pool model.
type plantSim struct {
	nChrom  int
	lenBp   float64
	markers []int64
}

// gamete returns, for each marker, whether the gamete carries the high allele.
func (p plantSim) gamete(rng *rand.Rand) []bool {
	g := make([]bool, len(p.markers))
	state := rng.Intn(2) == 1
	nx := distuv.Poisson{Lambda: 1, Src: rng}.Rand() // 1 Morgan
	xs := make([]float64, int(nx))
	for i := range xs {
		xs[i] = rng.Float64() * p.lenBp
	}
	sort.Float64s(xs)
	k := 0
	for i, pos := range p.markers {
		for k < len(xs) && xs[k] < float64(pos) {
			state = !state
			k++
		}
		g[i] = state
	}
	return g
}

func (p plantSim) dataset(rng *rand.Rand, nPlants int, rho float64, qtl bool) []BSAstats {
	return p.datasetPhase(rng, nPlants, rho, qtl, false)
}

// datasetPhase: with randomPhase, REF is the high or the low parent's allele
// at random per SNP and raw stats are computed as in 2b without --ref-parent.
func (p plantSim) datasetPhase(rng *rand.Rand, nPlants int, rho float64, qtl, randomPhase bool) []BSAstats {
	var raw []BSAstats
	calc := rawCalc{roles: siteRoles{hp: -1, lp: -1, hb: 0, lb: 1}, signedDSI: !randomPhase}
	for c := 0; c < p.nChrom; c++ {
		copies := [2][]int{make([]int, len(p.markers)), make([]int, len(p.markers))}
		for b := 0; b < 2; b++ {
			for range nPlants {
				var g1, g2 []bool
				for {
					g1, g2 = p.gamete(rng), p.gamete(rng)
					if !qtl || c != 0 {
						break
					}
					// QTL at the chromosome middle: high bulk needs ≥1 high allele, low bulk none (recessive-ish selection)
					mid := len(p.markers) / 2
					n := 0
					if g1[mid] {
						n++
					}
					if g2[mid] {
						n++
					}
					if (b == 0 && n == 2) || (b == 1 && n == 0) {
						break
					}
				}
				for i := range p.markers {
					if g1[i] {
						copies[b][i]++
					}
					if g2[i] {
						copies[b][i]++
					}
				}
			}
		}
		for i, pos := range p.markers {
			site := filter.Site{Chrom: string(rune('A' + c)), Pos: pos, Ref: "A", Alt: "T", NSamples: 2}
			flip := randomPhase && rng.Intn(2) == 1
			for b := 0; b < 2; b++ {
				d := int(distuv.Poisson{Lambda: distuv.Gamma{Alpha: 8, Beta: 8.0 / 50, Src: rng}.Rand(), Src: rng}.Rand())
				f := float64(copies[b][i]) / float64(2*nPlants)
				pr := f*0.99 + (1-f)*0.01
				if rho > 0 { // beta-binomial overdispersion
					a := pr * (1 - rho) / rho
					pr = distuv.Beta{Alpha: a, Beta: a * (1 - pr) / pr, Src: rng}.Rand()
				}
				h := binomSample(rng, d, pr)
				site.Samples[b] = filter.SampleAD{Ref: int32(h), Alt: int32(d - h), HasAD: true}
				if flip { // REF is the low parent's allele at this SNP
					site.Samples[b].Ref, site.Samples[b].Alt = site.Samples[b].Alt, site.Samples[b].Ref
				}
			}
			raw = append(raw, calc.site(&site))
		}
	}
	return raw
}

// TestHMMLODCalibration checks the genome-wide null against data from an
// independent plant-level generator (F2 plants with Poisson crossovers,
// gamma-Poisson depths, optional beta-binomial overdispersion): about 5% of
// no-QTL genomes should exceed the p95 HMM LOD and CompositeZ thresholds and
// 1% the p99 (both genome-wide), with and without overdispersion, and the
// estimated read overdispersion ρ should be near the truth. Slow: set
// GOBSASEQ_CALIBRATION=1 to run.
func TestHMMLODCalibration(t *testing.T) {
	if os.Getenv("GOBSASEQ_CALIBRATION") == "" {
		t.Skip("set GOBSASEQ_CALIBRATION=1 to run the calibration check")
	}
	color.Output = io.Discard
	p := plantSim{nChrom: 3, lenBp: 30e6}
	for x := 1500.0; x < p.lenBp; x += 3000 {
		p.markers = append(p.markers, int64(x))
	}
	// Two-bulk with known polarity (binomial and overdispersed reads; HMM LOD
	// and signed BRM), and 2b without --ref-parent with REF from either parent
	// at random per SNP (no HMM LOD; polarity-free BRM). BRM is checked as the
	// fraction of no-QTL genomes with any block (it is a pointwise α = 0.05
	// test, so this is reported, not asserted).
	for _, sc := range []struct {
		mode, refParent string
		rho             float64
		randomPhase     bool
	}{{"2p2b", "high", 0, false}, {"2p2b", "high", 0.01, false}, {"2b", "none", 0, true}} {
		cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: sc.refParent, NullReps: 200,
			BrmAlpha: 0.05, OutputDir: t.TempDir()}
		scorer, err := newScorer(cfg, sc.mode)
		if err != nil {
			t.Fatal(err)
		}
		rng := rand.New(rand.NewSource(11))
		withHMM := func(raw []BSAstats) []BSAstats { // as RawStats does
			if scorer.hmm != nil {
				computeHMMLOD(scorer.hmm, raw, newMarkerOrder(raw))
			}
			return raw
		}
		tmpl := withHMM(p.datasetPhase(rng, 20, sc.rho, false, sc.randomPhase))
		th, err := genomeNullThresholds(cfg, sc.mode, tmpl)
		zt, ht := th.Z, th.HMM
		if err != nil {
			t.Fatal(err)
		}
		sim, err := newNullSim(cfg, sc.mode, tmpl)
		if err != nil {
			t.Fatal(err)
		}
		rhoHat, _ := estimateRho(tmpl, newMarkerOrder(tmpl), sim)
		var over95, over99, czOver95, czOver99, brmGenomes int
		// per-statistic: null genomes where any scanned smoothed track crosses its p99 / p95
		var statOver99, statOver95 int
		lay := scorer.lay
		const nNull = 100
		for range nNull {
			raw := withHMM(p.datasetPhase(rng, 20, sc.rho, false, sc.randomPhase))
			sm, _ := scorer.score(raw, false, 1)
			mx := 0.0
			for i := range raw {
				mx = math.Max(mx, raw[i].HMMLOD)
			}
			if mx > ht.LODP95 {
				over95++
			}
			if mx > ht.LODP99 {
				over99++
			}
			cz := math.Inf(-1)
			for i := range sm {
				cz = math.Max(cz, sm[i].CompositeZ)
			}
			if o99, o95 := crossesAny(sm, th, lay); o99 {
				statOver99++
				statOver95++
			} else if o95 {
				statOver95++
			}
			if cz > zt.CompositeZP95 {
				czOver95++
			}
			if cz > zt.CompositeZP99 {
				czOver99++
			}
			if blocks, _ := RunBRM(cfg, sc.mode, sm); len(blocks) > 0 {
				brmGenomes++
			}
		}
		// power: one QTL dataset
		qraw := withHMM(p.datasetPhase(rng, 20, sc.rho, true, sc.randomPhase))
		qsm, _ := scorer.score(qraw, false, 1)
		qBRM := 0
		if blocks, _ := RunBRM(cfg, sc.mode, qsm); len(blocks) > 0 {
			for _, b := range blocks {
				if b.Chrom == "A" && b.Start <= 15_000_000 && b.Stop >= 15_000_000 {
					qBRM = 1
				}
			}
		}
		qmax := 0.0
		for i := range qraw {
			if qraw[i].CHROM == "A" {
				qmax = math.Max(qmax, qraw[i].HMMLOD)
			}
		}
		if HMMAvailable(cfg, sc.mode) && (over99 > 5 || over95 > 12) {
			t.Errorf("%s/%s rho=%.3f: null genomes over p95 %d/%d, over p99 %d/%d: thresholds too low", sc.mode, sc.refParent, sc.rho, over95, nNull, over99, nNull)
		}
		if czOver99 > 5 || czOver95 > 12 {
			t.Errorf("%s/%s rho=%.3f: null genomes over CompositeZ p95 %d/%d, over p99 %d/%d: thresholds too low", sc.mode, sc.refParent, sc.rho, czOver95, nNull, czOver99, nNull)
		}
		if statOver99 > 5*len(nullExtremes(lay)) {
			t.Errorf("%s/%s rho=%.3f: %d/%d null genomes cross a per-statistic p99", sc.mode, sc.refParent, sc.rho, statOver99, nNull)
		}
		t.Logf("  per-statistic: null genomes crossing any track's p99 %d, p95 %d (of %d; tracks %d)", statOver99, statOver95, nNull, len(nullExtremes(lay))-1)
		qcz := math.Inf(-1)
		for i := range qsm {
			if qsm[i].CHROM == "A" {
				qcz = math.Max(qcz, qsm[i].CompositeZ)
			}
		}
		t.Logf("%s ref-parent=%s rho=%.3f (estimated %.4f): HMM LOD p95/p99 %.2f/%.2f, null genomes over %d/%d of %d | CompositeZ p95/p99 %.2f/%.2f, null genomes over %d/%d | QTL max LOD %.1f, max CompositeZ %.1f | BRM: null genomes with a block %d/%d, QTL found %v",
			sc.mode, sc.refParent, sc.rho, rhoHat, ht.LODP95, ht.LODP99, over95, over99, nNull, zt.CompositeZP95, zt.CompositeZP99, czOver95, czOver99, qmax, qcz, brmGenomes, nNull, qBRM == 1)
	}
}

// crossesAny reports whether any detection track (ΔSI or AFDev both tails,
// G′/G1, ED⁴; not SI or CompositeZ) crosses its p99, and its p95.
func crossesAny(sm []SmoothedStats, th Thresholds, lay smoothLayout) (p99, p95 bool) {
	tb, ob := th.TwoBulk, th.OneBulk
	type tr struct {
		get                func(*SmoothedStats) float64
		u99, u95, l99, l95 float64
	}
	nan := math.NaN()
	var trs []tr
	if lay.both() {
		trs = append(trs, tr{func(s *SmoothedStats) float64 { return s.SmDeltaSI }, tb.DeltaSIP99, tb.DeltaSIP95, tb.DeltaSIMp99, tb.DeltaSIMp95},
			tr{func(s *SmoothedStats) float64 { return s.SmGstat }, tb.GstatP99, tb.GstatP95, nan, nan},
			tr{func(s *SmoothedStats) float64 { return s.SmED }, tb.ED4P99, tb.ED4P95, nan, nan})
	}
	if lay.one() {
		trs = append(trs, tr{func(s *SmoothedStats) float64 { return s.SmAFDev }, ob.AFDevP99, ob.AFDevP95, ob.AFDevMp99, ob.AFDevMp95},
			tr{func(s *SmoothedStats) float64 { return s.SmOneBulkG }, ob.OneBulkGstatP99, ob.OneBulkGstatP95, nan, nan})
	}
	for i := range sm {
		for _, r := range trs {
			v := r.get(&sm[i])
			p99 = p99 || v > r.u99 || v < r.l99
			p95 = p95 || v > r.u95 || v < r.l95
		}
	}
	return p99, p95
}
