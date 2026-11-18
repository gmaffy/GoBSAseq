package stats

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"testing"
)

func TestFindPeakIntersectionsUsesLocalThresholds(t *testing.T) {
	sm := []SmoothedStats{
		{CHROM: "chr1", POS: 100, SmGstat: 0.1},
		{CHROM: "chr1", POS: 200, SmGstat: 0.6},
	}
	th := []Thresholds{
		{TwoBulk: TwoBulkThresholds{GstatP99: 0.2}},
		{TwoBulk: TwoBulkThresholds{GstatP99: 0.8}},
	}
	peaks := FindPeakIntersections("Gstat", sm, th,
		func(s SmoothedStats) float64 { return s.SmGstat },
		func(t Thresholds) float64 { return t.TwoBulk.GstatP99 }, UpperTail)
	if len(peaks) != 0 {
		t.Fatalf("expected no peak below its local threshold, got %#v", peaks)
	}
}

func TestMergeCompositeBRMPreservesDisconnectedLoci(t *testing.T) {
	peaks := []PeakIntersection{
		{Chrom: "chr1", Start: 100, End: 200, PeakPos: 150, PeakValue: 5, ThresholdLevel: "p99"},
		{Chrom: "chr1", Start: 500, End: 600, PeakPos: 550, PeakValue: 4, ThresholdLevel: "p99"},
	}
	merged := MergeCompositeBRM(peaks, nil, 0)
	if len(merged) != 2 {
		t.Fatalf("expected two disconnected QTLs, got %#v", merged)
	}
}

func TestMergeCompositeBRMBridgesNearbyLoci(t *testing.T) {
	peaks := []PeakIntersection{
		{Chrom: "chr1", Start: 100, End: 200, PeakPos: 150, PeakValue: 5, ThresholdLevel: "p99"},
		{Chrom: "chr1", Start: 250, End: 350, PeakPos: 300, PeakValue: 4, ThresholdLevel: "p99"},
	}
	merged := MergeCompositeBRM(peaks, nil, 50)
	if len(merged) != 1 {
		t.Fatalf("expected one merged QTL, got %#v", merged)
	}
	if merged[0].Start != 100 || merged[0].Stop != 350 {
		t.Fatalf("unexpected merged span: %#v", merged[0])
	}
}

func TestConsolidateQTLsMergesWithinGap(t *testing.T) {
	peaks := []PeakIntersection{
		{Chrom: "chr1", Start: 100, End: 200, PeakPos: 150, PeakValue: 2, ThresholdLevel: "p99", Stat: "Gstat"},
		{Chrom: "chr1", Start: 260, End: 360, PeakPos: 310, PeakValue: 3, ThresholdLevel: "p99", Stat: "LOD"},
	}
	qtls := ConsolidateQTLs(peaks, 100)
	if len(qtls) != 1 {
		t.Fatalf("expected one consolidated QTL, got %#v", qtls)
	}
}

func TestCrossSimulationMatchesExpectedFrequencies(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	pos := []float64{1, 5e6, 15e6, 29e6}
	cases := []struct {
		pop       string
		wantAF    float64
		fixedOnly bool // every individual homozygous at every marker
	}{
		{"F2", 0.5, false},
		{"RIL", 0.5, false},
		{"BC1H", 0.75, false},
		{"BC2L", 0.125, false},
		{"BC1F2H", 0.75, false},
	}
	for _, c := range cases {
		d, err := parseCrossDesign(c.pop)
		if err != nil {
			t.Fatal(err)
		}
		exp, _ := ExpectedAF(c.pop)
		if math.Abs(exp-c.wantAF) > 1e-12 {
			t.Fatalf("%s: ExpectedAF=%v, test expects %v", c.pop, exp, c.wantAF)
		}
		const n = 4000
		diff := make([]int, len(pos)+1)
		cnt := make([]int, len(pos))
		poolHighAlleleCount(rng, d, n, pos, 30e6, morgansPerChromosome, diff, cnt)
		for i, ci := range cnt {
			f := float64(ci) / (2 * n)
			if math.Abs(f-c.wantAF) > 0.03 {
				t.Errorf("%s marker %d: pool AF %.3f, want ≈ %.3f", c.pop, i, f, c.wantAF)
			}
		}
	}

	// F2 individuals have dosage 0, 1 or 2; RILs are (almost) fully homozygous.
	f2, _ := parseCrossDesign("F2")
	ril, _ := parseCrossDesign("RIL")
	het := 0
	for range 2000 {
		ind := simulateIndividual(rng, ril, 30e6, morgansPerChromosome)
		if ind[0].alleleAt(15e6) != ind[1].alleleAt(15e6) {
			het++
		}
		_ = simulateIndividual(rng, f2, 30e6, morgansPerChromosome)
	}
	if het > 40 { // expected ~1/128 of 2000 ≈ 16
		t.Errorf("RIL heterozygosity too high: %d / 2000", het)
	}
}

func TestGameteRecombinationRate(t *testing.T) {
	// Recombination fraction between two markers d Morgans apart must follow
	// Haldane: r = (1 - exp(-2d)) / 2.
	rng := rand.New(rand.NewSource(2))
	H, L := haplotype{start: 1}, haplotype{start: 0}
	const length = 100e6 // 1 Morgan over 100 Mb -> 0.01 M per Mb
	a, b := 10e6, 40e6 // 0.3 M apart
	rec := 0
	const n = 20000
	for range n {
		g := gamete(rng, H, L, length, 1.0)
		if g.alleleAt(a) != g.alleleAt(b) {
			rec++
		}
	}
	want := (1 - math.Exp(-2*0.3)) / 2
	if got := float64(rec) / n; math.Abs(got-want) > 0.015 {
		t.Errorf("recombination fraction %.3f, want ≈ %.3f", got, want)
	}
}

// TestNullThresholdCalibration checks that thresholds from one set of null
// genomes hold their genome-wide error rate on an independent set.
func TestNullThresholdCalibration(t *testing.T) {
	if testing.Short() {
		t.Skip("slow")
	}
	var raw []BSAstats
	var sm []SmoothedStats
	for c := 0; c < 3; c++ {
		chrom := fmt.Sprintf("chr%d", c+1)
		for p := int64(5000); p < 20e6; p += 20000 {
			raw = append(raw, BSAstats{CHROM: chrom, POS: p, HighBulkH: 30, HighBulkL: 30, LowBulkH: 30, LowBulkL: 30, Depth: 60})
			sm = append(sm, SmoothedStats{CHROM: chrom, POS: p})
		}
	}
	d, _ := parseCrossDesign("F2")
	L := buildNullLayout(raw, sm, "2p2b", 1e6, nil, d, 0.5, 30, 30, 0)
	L.brmScale = PopulationVarianceScale("F2")

	calib := simulateNullExtremes(L, 400, 11)
	test := simulateNullExtremes(L, 400, 12)
	for _, c := range []struct {
		name string
		f    func(nullExtremes) float64
	}{
		{"CompositeZ", func(e nullExtremes) float64 { return e.CompositeZMax }},
		{"BRM |ΔSI|/σ", func(e nullExtremes) float64 { return e.BRMZMax }},
	} {
		xs := make([]float64, len(calib))
		for i, e := range calib {
			xs[i] = c.f(e)
		}
		sort.Float64s(xs)
		thr95 := xs[int(0.95*float64(len(xs)))]
		exceed := 0
		for _, e := range test {
			if c.f(e) > thr95 {
				exceed++
			}
		}
		// Expected 5% of 400 = 20; allow sampling noise.
		if exceed < 8 || exceed > 36 {
			t.Errorf("independent null genomes exceeding P95 %s: %d / 400, want ≈ 20", c.name, exceed)
		}
	}
}

// BRM's σ is pool sampling only and must not depend on read depth; with F2
// bulks of 50 at af = 0.5 it is sqrt(2 · 0.25 / (2 · 50)) ≈ 0.0707.
func TestBRMSigma(t *testing.T) {
	if s := brmSigmaTwoBulk(0.5, 50, 50, PopulationVarianceScale("F2")); math.Abs(s-math.Sqrt(0.005)) > 1e-12 {
		t.Errorf("two-bulk σ = %v, want %v", s, math.Sqrt(0.005))
	}
	if s := brmSigmaOneBulk(0.5, 50, PopulationVarianceScale("RIL")); math.Abs(s-math.Sqrt(0.005)) > 1e-12 {
		t.Errorf("one-bulk σ = %v, want %v", s, math.Sqrt(0.005))
	}
}
