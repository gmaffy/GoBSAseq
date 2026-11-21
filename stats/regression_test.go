package stats

import (
	"math"
	"math/rand"
	"testing"
)

func TestFindPeakIntersectionsGenomeWideThreshold(t *testing.T) {
	sm := []SmoothedStats{
		{CHROM: "chr1", POS: 100, SmGstat: 0.1},
		{CHROM: "chr1", POS: 200, SmGstat: 0.6},
		{CHROM: "chr1", POS: 300, SmGstat: 0.9},
		{CHROM: "chr1", POS: 400, SmGstat: 0.3},
	}
	g := func(s SmoothedStats) float64 { return s.SmGstat }
	peaks := FindPeakIntersections("Gstat", sm, 0.5, g, UpperTail)
	if len(peaks) != 1 || peaks[0].PeakPos != 300 || peaks[0].Threshold != 0.5 ||
		math.Abs(peaks[0].Start-180) > 1e-9 || math.Abs(peaks[0].End-(300+100*0.4/0.6)) > 1e-9 {
		t.Fatalf("peaks %#v", peaks)
	}
	if p := FindPeakIntersections("Gstat", sm, math.NaN(), g, UpperTail); p != nil {
		t.Fatalf("NaN threshold (no tail) gave peaks %#v", p)
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

func TestPoolSamplerRespectsDiploidAndRILGenotypes(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for range 100 {
		f2 := samplePoolAF(rng, "F2", 1, 0.5)
		if f2 != 0 && f2 != 0.5 && f2 != 1 {
			t.Fatalf("F2 dosage frequency must be 0, 0.5, or 1; got %g", f2)
		}
		ril := samplePoolAF(rng, "RIL", 1, 0.5)
		if ril != 0 && ril != 1 {
			t.Fatalf("RIL genotype must be fixed; got %g", ril)
		}
	}
}
