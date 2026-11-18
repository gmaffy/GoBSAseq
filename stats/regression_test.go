package stats

import (
	"math/rand"
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
