package stats

import (
	"testing"
)

func TestFindZSegmentsBridgesShortGapsAndKeepsSingleRegion(t *testing.T) {
	variants := []SmoothedStats{
		{CHROM: "1", POS: 100, CompositeZ: 2.0},
		{CHROM: "1", POS: 101, CompositeZ: 2.2},
		{CHROM: "1", POS: 102, CompositeZ: 0.8},
		{CHROM: "1", POS: 103, CompositeZ: 2.1},
		{CHROM: "1", POS: 104, CompositeZ: 2.3},
	}
	thresholds := make([]Thresholds, len(variants))
	for i := range thresholds {
		thresholds[i] = Thresholds{}
	}

	records := findZSegments(variants, thresholds, true, func(Thresholds) float64 { return 1.0 })
	if len(records) != 1 {
		t.Fatalf("expected one segmented QTL, got %d", len(records))
	}
	if records[0].Start != 100 || records[0].Stop != 104 {
		t.Fatalf("expected run boundaries 100-104, got %d-%d", records[0].Start, records[0].Stop)
	}
}

func TestSelectBestQTLPrefersStrongestValleyOrPeak(t *testing.T) {
	records := []QTLRecord{
		{Chrom: "1", Start: 100, Stop: 110, PeakPos: 105, PeakVal: 3.2},
		{Chrom: "1", Start: 200, Stop: 210, PeakPos: 205, PeakVal: -4.5},
	}

	best := selectBestQTL(records)
	if best.PeakPos != 205 || best.PeakVal != -4.5 {
		t.Fatalf("expected the strongest valley to be selected, got pos=%d val=%.2f", best.PeakPos, best.PeakVal)
	}
	if best.Start != 200 || best.Stop != 210 {
		t.Fatalf("expected the valley interval to be preserved, got %d-%d", best.Start, best.Stop)
	}
}

func TestSelectChromQTLPoolPrefersP99OverP95(t *testing.T) {
	variants := []SmoothedStats{
		{CHROM: "1", POS: 100, CompositeZ: 4.0},
		{CHROM: "1", POS: 101, CompositeZ: 4.2},
		{CHROM: "1", POS: 102, CompositeZ: 4.1},
		{CHROM: "1", POS: 103, CompositeZ: 4.4},
		{CHROM: "1", POS: 104, CompositeZ: 4.3},
		{CHROM: "1", POS: 105, CompositeZ: 2.4},
		{CHROM: "1", POS: 106, CompositeZ: 2.5},
		{CHROM: "1", POS: 107, CompositeZ: 2.6},
		{CHROM: "1", POS: 108, CompositeZ: 2.7},
		{CHROM: "1", POS: 109, CompositeZ: 2.8},
	}
	thresholds := make([]Thresholds, len(variants))
	for i := range thresholds {
		thresholds[i] = Thresholds{Z: ZThresholds{CompositeZP99: 3.5, CompositeZP95: 2.0}}
	}

	pool, stat, thr := selectChromQTLPool(variants, thresholds)
	if stat != "Stouffer" || thr != "P99" {
		t.Fatalf("expected P99 tier to be selected, got stat=%s threshold=%s", stat, thr)
	}
	if len(pool) != 1 {
		t.Fatalf("expected only the P99 region to be retained, got %d regions", len(pool))
	}
	if pool[0].PeakPos != 103 || pool[0].NQTLs != 1 {
		t.Fatalf("expected the P99 region to be kept with N_QTLs=1, got peak=%d nqtls=%d", pool[0].PeakPos, pool[0].NQTLs)
	}
}
