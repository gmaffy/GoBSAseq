package stats

import (
	"math"
	"math/rand"
	"strconv"
	"strings"
	"testing"

	"github.com/gmaffy/GoBSAseq/utils"
)

// The pool AF sampler must reproduce the analytic mean (ExpectedAF) and
// variance p(1-p)/(scale*n) (PopulationVarianceScale) for every population.
func TestSamplePoolAFMatchesAnalyticMoments(t *testing.T) {
	const n, reps = 20, 200000
	for _, pop := range []string{"F2", "F3", "F4", "RIL", "BC1H", "BC1L", "BC2H", "BC1F2H", "BC1F2L", "BC2F2H", "BC3F2L"} {
		p0, err := ExpectedAF(pop)
		if err != nil {
			t.Fatal(err)
		}
		wantVar := p0 * (1 - p0) / (PopulationVarianceScale(pop) * n)
		rng := rand.New(rand.NewSource(1))
		var sum, ss float64
		for range reps {
			f := samplePoolAF(rng, pop, n, p0)
			sum += f
			ss += f * f
		}
		mean := sum / reps
		v := ss/reps - mean*mean
		if math.Abs(mean-p0) > 0.002 {
			t.Errorf("%s: mean %.4f, want %.4f", pop, mean, p0)
		}
		if math.Abs(v-wantVar)/wantVar > 0.03 {
			t.Errorf("%s: variance %.6f, want %.6f", pop, v, wantVar)
		}
	}
}

func TestBC1F2HasDonorHomozygotes(t *testing.T) {
	rng := rand.New(rand.NewSource(2))
	donor := 0
	const reps = 100000
	for range reps {
		if samplePoolAF(rng, "BC1F2H", 1, 0.75) == 0 {
			donor++
		}
	}
	if got := float64(donor) / reps; math.Abs(got-0.125) > 0.005 {
		t.Fatalf("BC1F2 donor-homozygote frequency %.4f, want 0.125", got)
	}
}

func TestBulkSizesPerMode(t *testing.T) {
	cfg := utils.AnalysisConfig{HighBulkSize: 20, LowBulkSize: 16}
	cases := map[string]int{"2phb": 20, "hphb": 20, "lphb": 20, "2plb": 16, "hplb": 16, "lplb": 16, "2p2b": 0, "2b": 0}
	for mode, want := range cases {
		_, _, one, err := cfg.BulkSizes(mode)
		if err != nil {
			t.Fatalf("%s: %v", mode, err)
		}
		if one != want {
			t.Errorf("%s: one-bulk size %d, want %d", mode, one, want)
		}
	}
	if _, _, _, err := (utils.AnalysisConfig{HighBulkSize: 20}).BulkSizes("hplb"); err == nil {
		t.Error("expected error for missing low bulk size in hplb")
	}
	if _, _, _, err := (utils.AnalysisConfig{LowBulkSize: 20}).BulkSizes("2p2b"); err == nil {
		t.Error("expected error for missing high bulk size in 2p2b")
	}
}

func TestCalculateThresholdsRejectsMissingBulkSize(t *testing.T) {
	cfg := utils.AnalysisConfig{Population: "F2", OutputDir: t.TempDir(), HighBulkSize: 20}
	raw := []BSAstats{{CHROM: "chr1", POS: 1, Depth: 30, LowDepth: 30}}
	if _, err := CalculateThresholds(cfg, "lplb", raw); err == nil {
		t.Fatal("expected error when the low bulk size is missing")
	}
}

// samplePoolAF draws the high-parent allele frequency of a pool of bulkSize
// diploid plants by genotype dosage (F2/F3/F4 selfing, RIL, BCn, BCnF2): a
// plant-level reference for PopulationVarianceScale, which the HMM pool model
// and the marker weights use.
func samplePoolAF(rng *rand.Rand, population string, bulkSize int, p0 float64) float64 {
	p := strings.ToUpper(strings.TrimSpace(population))
	selfGenerations := map[string]int{"F2": 1, "F3": 2, "F4": 3}
	if g, ok := selfGenerations[p]; ok || p == "RIL" {
		het := 0.0
		if p != "RIL" {
			het = math.Pow(0.5, float64(g))
		}
		pHom := (1 - het) / 2
		copies := 0
		for range bulkSize {
			u := rng.Float64()
			switch {
			case u < pHom:
			case u < pHom+het:
				copies++
			default:
				copies += 2
			}
		}
		return float64(copies) / float64(2*bulkSize)
	}
	// BCn[F2]H/L: recurrent homozygote, heterozygote q = 0.5^n (BCnF2: q/2,
	// donor homozygote q/4); SI is the recurrent allele for H, donor for L.
	s := p[2:]
	isHigh := strings.HasSuffix(s, "H")
	s = strings.TrimSuffix(strings.TrimSuffix(s, "H"), "L")
	isF2 := strings.HasSuffix(s, "F2")
	n, _ := strconv.Atoi(strings.TrimSuffix(s, "F2"))
	q := math.Pow(0.5, float64(n))
	pHet, pDonor := q, 0.0
	if isF2 {
		pHet, pDonor = q/2, q/4
	}
	recurrent := 0
	for range bulkSize {
		u := rng.Float64()
		switch {
		case u < pDonor:
		case u < pDonor+pHet:
			recurrent++
		default:
			recurrent += 2
		}
	}
	f := float64(recurrent) / float64(2*bulkSize)
	if !isHigh {
		f = 1 - f
	}
	return f
}
