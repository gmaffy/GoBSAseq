package stats

import (
	"math"
	"math/rand"
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
	cfg := utils.AnalysisConfig{Population: "F2", Rep: 10, OutputDir: t.TempDir(), HighBulkSize: 20}
	sm := []SmoothedStats{{CHROM: "chr1", POS: 1, Depth: 30}}
	if _, err := CalculateThresholds(cfg, "lplb", nil, sm); err == nil {
		t.Fatal("expected error when the low bulk size is missing")
	}
}
