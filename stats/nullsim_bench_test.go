package stats

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"testing"
)

// syntheticLayout builds a null layout for a genome of nChrom chromosomes of
// chromLen bp with one marker every spacing bp (2p2b, F2, depth ~60, σ = 1 Mb,
// the default 2 Mb window).
func syntheticLayout(nChrom int, chromLen, spacing int64) *nullLayout {
	rng := rand.New(rand.NewSource(1))
	var raw []BSAstats
	var sm []SmoothedStats
	lens := map[string]float64{}
	for c := 0; c < nChrom; c++ {
		name := fmt.Sprintf("chr%02d", c+1)
		lens[name] = float64(chromLen)
		for p := int64(1000); p < chromLen; p += spacing {
			h, l := 30+rng.Intn(30), 30+rng.Intn(30)
			raw = append(raw, BSAstats{CHROM: name, POS: p, HighBulkH: h / 2, HighBulkL: h - h/2, LowBulkH: l / 2, LowBulkL: l - l/2, Depth: h + l})
			sm = append(sm, SmoothedStats{CHROM: name, POS: p})
		}
	}
	design, _ := parseCrossDesign("F2")
	return buildNullLayout(raw, sm, "2p2b", 1e6, lens, design, 0.5, 50, 50, 0)
}

// go test ./stats -run XXX -bench NullGenome -benchtime 5x
func BenchmarkNullGenome(b *testing.B) {
	cases := []struct {
		name            string
		nChrom          int
		chromLen, space int64
	}{
		{"150k_markers_150Mb", 5, 30_000_000, 1_000},
		{"cucumber_350Mb_7chr_350k", 7, 50_000_000, 1_000},
		{"tomato_900Mb_12chr_900k", 12, 75_000_000, 1_000},
		{"tomato_900Mb_12chr_3M", 12, 75_000_000, 300},
	}
	for _, c := range cases {
		L := syntheticLayout(c.nChrom, c.chromLen, c.space)
		b.Run(c.name, func(b *testing.B) {
			w := &nullWorker{}
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				L.simulateGenome(rand.New(rand.NewSource(int64(i))), w)
			}
		})
	}
}

// weightedMedian on (value, count) pairs must equal the median of the
// expanded vector, including ties and even totals.
func TestWeightedMedianMatchesExpanded(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	var s robustZScratch
	for trial := 0; trial < 300; trial++ {
		n := 1 + rng.Intn(40)
		vals := make([]float64, n)
		counts := make([]int, n)
		var expanded []float64
		for i := range vals {
			if trial%2 == 0 {
				vals[i] = rng.NormFloat64()
			} else {
				vals[i] = float64(rng.Intn(4)) // many ties
			}
			counts[i] = 1 + rng.Intn(6)
			for range counts[i] {
				expanded = append(expanded, vals[i])
			}
		}
		sort.Float64s(expanded)
		m := len(expanded)
		want := expanded[m/2]
		if m%2 == 0 {
			want = (expanded[m/2-1] + expanded[m/2]) / 2
		}
		if got := s.weightedMedian(vals, counts, m); got != want {
			t.Fatalf("trial %d: got %v, want %v", trial, got, want)
		}
	}
}

// Null scales come from QTL-free simulated genomes: ΔSI is centred near 0,
// and Z-scores of the threshold genomes use the same fixed transform.
func TestEstimateScales(t *testing.T) {
	L := syntheticLayout(3, 20_000_000, 20_000)
	L.brmScale = PopulationVarianceScale("F2")
	sc := L.estimateScales(50, 1)
	if m := sc[tDeltaSI].Med; math.Abs(m) > 0.02 {
		t.Errorf("null median of ΔSI = %v, want ≈ 0", m)
	}
	for _, tr := range L.activeTracks() {
		if !(sc[tr].Scale > 0) {
			t.Errorf("track %d: scale %v, want > 0", tr, sc[tr].Scale)
		}
	}
	// Same seed → same scales (reproducibility).
	if again := L.estimateScales(50, 1); again != sc {
		t.Error("estimateScales is not reproducible for a fixed seed")
	}
}

// The tabulated beta-binomial BF must equal betaBinomialLogBF bit for bit.
func TestLgammaTablesMatchBetaBinomialLogBF(t *testing.T) {
	tab := newLgammaTables(160)
	for hs := 0; hs <= 80; hs += 3 {
		for hf := 0; hf <= 80; hf += 7 {
			for ls := 0; ls <= 80; ls += 5 {
				for lf := 0; lf <= 80; lf += 11 {
					if got, want := tab.betaBinomialLogBF(hs, hf, ls, lf), betaBinomialLogBF(hs, hf, ls, lf); got != want {
						t.Fatalf("(%d,%d,%d,%d): got %v want %v", hs, hf, ls, lf, got, want)
					}
				}
			}
		}
	}
}

// The tabulated read sampler must reproduce Binomial(depth, c/2n): compare
// sample mean and variance with depth·p and depth·p·(1-p).
func TestReadSamplerMoments(t *testing.T) {
	rng := rand.New(rand.NewSource(5))
	s := newReadSampler(100, []int{7, 60, 250, 900}) // 900 > maxTableDepth: gonum fallback
	for _, depth := range []int{7, 60, 250, 900} {
		for _, c := range []int{0, 1, 13, 50, 87, 100} {
			p := float64(c) / 100
			const n = 40000
			var sum, sum2 float64
			for range n {
				h, l := s.draw(rng, depth, c, false)
				if h+l != depth || h < 0 {
					t.Fatalf("bad draw %d+%d for depth %d", h, l, depth)
				}
				sum += float64(h)
				sum2 += float64(h) * float64(h)
			}
			mean := sum / n
			variance := sum2/n - mean*mean
			wantM, wantV := float64(depth)*p, float64(depth)*p*(1-p)
			if d := mean - wantM; d*d > 25*wantV/n+1e-12 { // 5 standard errors
				t.Errorf("depth %d c %d: mean %.4f, want %.4f", depth, c, mean, wantM)
			}
			if wantV > 0 && (variance < 0.95*wantV || variance > 1.05*wantV) {
				t.Errorf("depth %d c %d: variance %.4f, want %.4f", depth, c, variance, wantV)
			}
			// flip: high reads come from the low-parent allele frequency.
			if h, _ := s.draw(rng, depth, 0, true); h != depth {
				t.Errorf("flip of c=0 should give all %d reads, got %d", depth, h)
			}
		}
	}
}
