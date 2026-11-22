package stats

import (
	"math"
	"testing"
)

// With unknown polarity the ΔSI track is |ΔSI|: one-sided thresholds (99%/95%
// of the genome maxima) and no lower tail; signed ΔSI is two-tailed.
func TestNullExtremesAbsDeltaSI(t *testing.T) {
	maxima := make([]float64, 200)
	for i := range maxima {
		maxima[i] = 0.1 + 0.001*float64(i)
	}
	for _, signed := range []bool{true, false} {
		for _, e := range nullExtremes(smoothLayout{highBulk: true, lowBulk: true, signedDSI: signed}) {
			if e.get(&SmoothedStats{SmDeltaSI: 7}) != 7 {
				continue
			}
			if e.twoTailed != signed {
				t.Fatalf("signed=%v: ΔSI two-tailed = %v", signed, e.twoTailed)
			}
			e.max, e.min = maxima, maxima
			var th Thresholds
			e.apply(&th)
			tb := th.TwoBulk
			if signed != !math.IsNaN(tb.DeltaSIMp99) || !(tb.DeltaSIP99 > tb.DeltaSIP95 && tb.DeltaSIP95 > 0) {
				t.Fatalf("signed=%v: ΔSI thresholds %+v", signed, tb)
			}
		}
	}
}

func TestAbsDeltaSITracks(t *testing.T) {
	var names []string
	for _, c := range absDeltaSITracks {
		names = append(names, c.name)
		if c.name == "DeltaSI+" || c.name == "DeltaSI-" || c.name == "HMMLOD" {
			t.Fatalf("track %q in the |ΔSI| set (no signed ΔSI or HMM without polarity)", c.name)
		}
	}
	if len(absDeltaSITracks) != len(twoBulkTracks)-2 || names[1] != "AbsDeltaSI" {
		t.Fatalf("unexpected |ΔSI| track set %v", names)
	}
}
