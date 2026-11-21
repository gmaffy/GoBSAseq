package stats

import (
	"math"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// Consolidation keeps each statistic's strongest peak, whatever the order
// the regions arrive in, including lower-tail (negative) peaks.
func TestConsolidateKeepsStrongestPeak(t *testing.T) {
	peaks := []PeakIntersection{
		{Stat: "Gstat", Chrom: "c", Start: 100, End: 300, PeakPos: 200, PeakValue: 9, ThresholdLevel: "p99"},
		{Stat: "Gstat", Chrom: "c", Start: 250, End: 500, PeakPos: 400, PeakValue: 5, ThresholdLevel: "p99"},
		{Stat: "DeltaSI-", Chrom: "c", Start: 150, End: 350, PeakPos: 260, PeakValue: -0.5, ThresholdLevel: "p99"},
		{Stat: "DeltaSI-", Chrom: "c", Start: 320, End: 450, PeakPos: 420, PeakValue: -0.2, ThresholdLevel: "p99"},
	}
	for _, order := range [][]int{{0, 1, 2, 3}, {3, 2, 1, 0}, {1, 3, 0, 2}} {
		in := make([]PeakIntersection, len(order))
		for i, j := range order {
			in[i] = peaks[j]
		}
		q := ConsolidateQTLs(in, 0)
		if len(q) != 1 || q[0].Peaks[slotGstat].Pos != 200 || q[0].Peaks[slotDeltaSI].Pos != 260 || q[0].Peaks[slotDeltaSI].Value != -0.5 {
			t.Fatalf("order %v: %+v (Gstat %+v, ΔSI %+v)", order, q, q[0].Peaks[slotGstat], q[0].Peaks[slotDeltaSI])
		}
	}
	// merging two consolidated QTLs keeps the stronger peak too
	m := mergeConsolidatedQTL(ConsolidatedQTL{Peaks: [nSlots]*Peak{slotGstat: {Pos: 1, Value: 3}}}, ConsolidatedQTL{Peaks: [nSlots]*Peak{slotGstat: {Pos: 2, Value: 7}}})
	if m.Peaks[slotGstat].Pos != 2 {
		t.Fatalf("merged peak %+v", m.Peaks[slotGstat])
	}
}

// Regions with fewer than minMarkers markers beyond the line are dropped,
// before the p95 fallback is decided.
func TestMinQTLMarkers(t *testing.T) {
	sm := make([]SmoothedStats, 30)
	for i := range sm {
		sm[i] = SmoothedStats{CHROM: "c", POS: int64(100 * (i + 1))}
	}
	for i := 5; i < 8; i++ { // 3 markers above p99
		sm[i].SmGstat = 10
	}
	for i := 15; i < 22; i++ { // 7 markers above p95 only
		sm[i].SmGstat = 6
	}
	g := func(s SmoothedStats) float64 { return s.SmGstat }
	chroms := byChromosome(sm)
	if p := scanGenome(chroms, "Gstat", 8, 5, g, UpperTail, 1, ""); len(p) != 1 || p[0].ThresholdLevel != "p99" || p[0].NMarkers != 3 {
		t.Fatalf("minMarkers 1: %+v", p)
	}
	p := scanGenome(chroms, "Gstat", 8, 5, g, UpperTail, 5, "")
	if len(p) != 1 || p[0].ThresholdLevel != "p95" || p[0].NMarkers != 7 || p[0].PeakPos != 1600 {
		t.Fatalf("minMarkers 5: the 3-marker p99 region must be dropped and p95 used: %+v", p)
	}
}

// The p95 fallback is decided per statistic for the whole genome: a p99
// region anywhere means no p95 regions on other chromosomes.
func TestGenomeWideFallback(t *testing.T) {
	var sm []SmoothedStats
	for _, c := range []string{"a", "b"} {
		for i := range 30 {
			s := SmoothedStats{CHROM: c, POS: int64(100 * (i + 1))}
			if i >= 10 && i < 20 {
				s.SmGstat = map[string]float64{"a": 10, "b": 6}[c] // a: above p99; b: above p95 only
			}
			sm = append(sm, s)
		}
	}
	g := func(s SmoothedStats) float64 { return s.SmGstat }
	p := scanGenome(byChromosome(sm), "Gstat", 8, 5, g, UpperTail, 1, "")
	if len(p) != 1 || p[0].Chrom != "a" || p[0].ThresholdLevel != "p99" {
		t.Fatalf("with a p99 region on a: %+v", p)
	}
	p = scanGenome(byChromosome(sm), "Gstat", 20, 5, g, UpperTail, 1, "") // no p99 anywhere
	if len(p) != 2 || p[0].ThresholdLevel != "p95" || p[1].ThresholdLevel != "p95" {
		t.Fatalf("no p99 anywhere: %+v", p)
	}
}

// Opposite directions stay separate QTLs even when they overlap; unsigned
// evidence joins the signed QTL it overlaps most; lone unsigned evidence is
// direction "none".
func TestConsolidateByDirection(t *testing.T) {
	peaks := []PeakIntersection{
		{Stat: "DeltaSI+", Direction: "+", Chrom: "c", Start: 100, End: 300, PeakPos: 200, PeakValue: 0.4, ThresholdLevel: "p99"},
		{Stat: "DeltaSI-", Direction: "−", Chrom: "c", Start: 320, End: 500, PeakPos: 400, PeakValue: -0.5, ThresholdLevel: "p99"},
		{Stat: "Gstat", Chrom: "c", Start: 150, End: 330, PeakPos: 210, PeakValue: 9, ThresholdLevel: "p99"}, // mostly over +
		{Stat: "ED4", Chrom: "c", Start: 380, End: 520, PeakPos: 410, PeakValue: 1, ThresholdLevel: "p95"},   // over −
		{Stat: "HMMLOD", Chrom: "c", Start: 900, End: 950, PeakPos: 920, PeakValue: 7, ThresholdLevel: "p99"},
	}
	q := ConsolidateQTLs(peaks, 50)
	if len(q) != 3 {
		t.Fatalf("QTLs %+v", q)
	}
	plus, minus, none := q[0], q[1], q[2]
	if plus.Direction != "+" || plus.Peaks[slotGstat] == nil || plus.Peaks[slotED] != nil || plus.End != 330 ||
		minus.Direction != "−" || minus.Peaks[slotED] == nil || minus.Peaks[slotGstat] != nil || minus.End != 520 ||
		none.Direction != "none" || none.Peaks[slotHMMLOD] == nil || none.Start != 900 {
		t.Fatalf("plus %+v\nminus %+v\nnone %+v", plus, minus, none)
	}
}

// Detection output is ordered by chromosome then position and identical
// across runs.
func TestDetectionDeterministic(t *testing.T) {
	var sm []SmoothedStats
	for _, c := range []string{"chr9", "chr1", "chr3"} {
		for i := range 60 {
			v := 0.0
			if (i >= 10 && i < 20) || (i >= 40 && i < 50) {
				v = 20
			}
			sm = append(sm, SmoothedStats{CHROM: c, POS: int64(10_000 * (i + 1)), SmGstat: v, CompositeZ: v, HMMLOD: v, SmDeltaSI: v / 40})
		}
	}
	var th Thresholds
	th.TwoBulk.GstatP99, th.TwoBulk.GstatP95 = 10, 5
	th.TwoBulk.DeltaSIP99, th.TwoBulk.DeltaSIP95, th.TwoBulk.DeltaSIMp99, th.TwoBulk.DeltaSIMp95 = 0.3, 0.2, -0.3, -0.2
	th.TwoBulk.ED4P99, th.TwoBulk.ED4P95 = 1, 1
	th.HMM.LODP99, th.HMM.LODP95 = 5, 4
	th.Z.CompositeZP99, th.Z.CompositeZP95 = 10, 5
	_, q0 := DetectIndividualQTLs(sm, th, "2p2b", true, 1000, 5)
	c0 := DetectCompositeZQTLs(sm, th, 1000, 5)
	if len(q0) != 6 || len(c0) != 6 {
		t.Fatalf("individual %d, CompositeZ %d, want 6 each", len(q0), len(c0))
	}
	for i := 1; i < len(q0); i++ {
		a, b := q0[i-1], q0[i]
		if a.Chrom > b.Chrom || (a.Chrom == b.Chrom && a.Start >= b.Start) {
			t.Fatalf("individual QTLs out of order: %+v then %+v", a, b)
		}
	}
	for range 5 {
		_, q := DetectIndividualQTLs(sm, th, "2p2b", true, 1000, 5)
		c := DetectCompositeZQTLs(sm, th, 1000, 5)
		if !reflect.DeepEqual(q, q0) || !reflect.DeepEqual(c, c0) {
			t.Fatal("detection output differs between runs")
		}
	}
}

func TestQTLWritersReturnErrors(t *testing.T) {
	bad := filepath.Join(t.TempDir(), "missing", "x.xlsx")
	if err := WriteIndividualQTLs([]ConsolidatedQTL{{Chrom: "c"}}, "2p2b", bad); err == nil {
		t.Error("individual: expected an error")
	}
	if err := WriteCompositeZQTLs([]PeakIntersection{{Chrom: "c"}}, bad); err == nil {
		t.Error("CompositeZ: expected an error")
	}
	if err := WriteFinalQTLs([]FinalQTL{{Chrom: "c", BRM: &MethodPeak{Z: math.NaN()}}}, bad); err == nil {
		t.Error("final: expected an error")
	}
	good := filepath.Join(t.TempDir(), "ok.xlsx")
	if err := WriteFinalQTLs([]FinalQTL{{Chrom: "c", Consensus: true, BRM: &MethodPeak{Pos: 5, Value: 0.3, Z: 4}}}, good); err != nil {
		t.Fatal(err)
	}
	if err := WriteIndividualQTLs([]ConsolidatedQTL{{Chrom: "c", Threshold: "p95", Direction: "+", Peaks: [nSlots]*Peak{slotGstat: {Pos: 5, Value: 9}}}}, "2p2b", filepath.Join(filepath.Dir(good), "ind.xlsx")); err != nil {
		t.Fatal(err)
	}
	if b, err := os.ReadFile(filepath.Join(filepath.Dir(good), "ind.tsv")); err != nil || !strings.Contains(string(b), "p95 (suggestive)\t+") {
		t.Fatalf("individual TSV: %q, %v", b, err)
	}
	if err := WriteCompositeZQTLs([]PeakIntersection{{Chrom: "c", ThresholdLevel: "p99", PeakPos: 7, PeakValue: 20}}, filepath.Join(filepath.Dir(good), "cz.xlsx")); err != nil {
		t.Fatal(err)
	}
	if b, err := os.ReadFile(filepath.Join(filepath.Dir(good), "cz.tsv")); err != nil || !strings.Contains(string(b), "p99\t7\t20") {
		t.Fatalf("CompositeZ TSV: %q, %v", b, err)
	}
	for _, f := range []string{good, strings.TrimSuffix(good, ".xlsx") + ".tsv"} {
		if fi, err := os.Stat(f); err != nil || fi.Size() == 0 {
			t.Fatalf("%s not written", f)
		}
	}
}
