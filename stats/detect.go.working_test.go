package stats

import "testing"

func TestConsolidateQTLs_chr04LikeOverlaps(t *testing.T) {
	peaks := []PeakIntersection{
		{Stat: "Gstat", Chrom: "Cp4.1LG04", Start: 5193652, End: 7359659, PeakPos: 5193652, PeakValue: 36.57, ThresholdLevel: "p95"},
		{Stat: "Gstat", Chrom: "Cp4.1LG04", Start: 5193652, End: 5193652, PeakPos: 5193652, PeakValue: 36.57, ThresholdLevel: "p95"},
		{Stat: "LOD", Chrom: "Cp4.1LG04", Start: 5215748, End: 7343499, PeakPos: 5215748, PeakValue: 1.0, ThresholdLevel: "p95"},
		{Stat: "BBLogBF", Chrom: "Cp4.1LG04", Start: 5215748, End: 7343499, PeakPos: 5215748, PeakValue: 2.0, ThresholdLevel: "p95"},
	}

	qtls := ConsolidateQTLs(peaks)
	if len(qtls) != 1 {
		t.Fatalf("expected 1 consolidated QTL, got %d: %+v", len(qtls), qtls)
	}
	if qtls[0].Gstat == nil || qtls[0].LOD == nil || qtls[0].BBLogBF == nil {
		t.Fatalf("expected merged evidence, got %+v", qtls[0])
	}
}

func TestConsolidateQTLs_floatStartGap(t *testing.T) {
	peaks := []PeakIntersection{
		{Stat: "Gstat", Chrom: "Cp4.1LG04", Start: 5193652, End: 5193652, PeakPos: 5193652, PeakValue: 36.57, ThresholdLevel: "p95"},
		{Stat: "Gstat", Chrom: "Cp4.1LG04", Start: 5193652.000001, End: 7359659, PeakPos: 5193652, PeakValue: 36.57, ThresholdLevel: "p95"},
		{Stat: "LOD", Chrom: "Cp4.1LG04", Start: 5215748, End: 7343499, PeakPos: 5215748, PeakValue: 1.0, ThresholdLevel: "p95"},
	}

	qtls := ConsolidateQTLs(peaks)
	if len(qtls) != 1 {
		t.Fatalf("expected 1 consolidated QTL, got %d: %+v", len(qtls), qtls)
	}
}

func TestConsolidateQTLs_bridgedByMiddlePeak(t *testing.T) {
	peaks := []PeakIntersection{
		{Stat: "Gstat", Chrom: "Cp4.1LG04", Start: 5193652, End: 5200000, PeakPos: 5193652, PeakValue: 1, ThresholdLevel: "p95"},
		{Stat: "ED4", Chrom: "Cp4.1LG04", Start: 7000000, End: 7100000, PeakPos: 7050000, PeakValue: 2, ThresholdLevel: "p95"},
		{Stat: "LOD", Chrom: "Cp4.1LG04", Start: 5215748, End: 7343499, PeakPos: 5215748, PeakValue: 3, ThresholdLevel: "p95"},
	}

	qtls := ConsolidateQTLs(peaks)
	if len(qtls) != 2 {
		t.Fatalf("expected 2 consolidated QTLs, got %d: %+v", len(qtls), qtls)
	}
}
