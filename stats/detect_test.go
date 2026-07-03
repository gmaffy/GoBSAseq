package stats

import "testing"

func TestDetectCompositeZQTLs_upperTail(t *testing.T) {
	smoothed := []SmoothedStats{
		{CHROM: "chr1", POS: 100, CompositeZ: 0.5},
		{CHROM: "chr1", POS: 200, CompositeZ: 3.0},
		{CHROM: "chr1", POS: 300, CompositeZ: 2.5},
		{CHROM: "chr1", POS: 400, CompositeZ: 0.5},
	}
	thresholds := make([]Thresholds, len(smoothed))
	for i := range thresholds {
		thresholds[i].Z.CompositeZP99 = 2.0
		thresholds[i].Z.CompositeZP95 = 1.5
	}

	peaks := DetectCompositeZQTLs(smoothed, thresholds, "2b")
	if len(peaks) != 1 {
		t.Fatalf("expected 1 CompositeZ peak, got %d: %+v", len(peaks), peaks)
	}
	if peaks[0].ThresholdLevel != "p99" {
		t.Fatalf("expected p99 threshold level, got %q", peaks[0].ThresholdLevel)
	}
	if peaks[0].PeakPos != 200 || peaks[0].PeakValue != 3.0 {
		t.Fatalf("unexpected peak: %+v", peaks[0])
	}
}

func TestMergeCompositeBRM_unionsIntervals(t *testing.T) {
	composite := []PeakIntersection{
		{
			Chrom: "chr1", Start: 1000, End: 2000,
			PeakPos: 1500, PeakValue: 3.5, ThresholdLevel: "p99",
		},
	}
	brm := []BRMBlock{
		{Chrom: "chr1", Start: 1800, Stop: 3000, PeakPos: 2200, Peak: 0.8, Threshold: 0.05},
	}

	merged := MergeCompositeBRM(composite, brm)
	if len(merged) != 1 {
		t.Fatalf("expected 1 merged interval, got %d", len(merged))
	}
	m := merged[0]
	if m.Source != "CompositeZ+BRM" {
		t.Fatalf("expected CompositeZ+BRM source, got %q", m.Source)
	}
	if m.Start != 1000 || m.Stop != 3000 {
		t.Fatalf("expected union interval 1000-3000, got %d-%d", m.Start, m.Stop)
	}
	if m.CompositeZPeak != 3.5 || m.BRMPeak != 0.8 {
		t.Fatalf("unexpected peaks: %+v", m)
	}
}
