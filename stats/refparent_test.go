package stats

import (
	"math"
	"testing"
)

// With unknown polarity the ΔSI null is |ΔSI|: the upper thresholds are
// non-negative, match the signed two-sided ones, and there is no lower tail.
func TestSimulateTwoBulkAbsDeltaSI(t *testing.T) {
	p := simParams{highDepth: 50, lowDepth: 50, p0High: 0.5, p0Low: 0.5, rep: 20000,
		highBulkSize: 20, lowBulkSize: 20, population: "F2"}
	signed := simulateTwoBulk(p)
	p.absDeltaSI = true
	abs := simulateTwoBulk(p)

	if !math.IsNaN(abs.DeltaSIMp99) || !math.IsNaN(abs.DeltaSIMp95) {
		t.Fatalf("|ΔSI| lower-tail thresholds = %v/%v, want NaN", abs.DeltaSIMp99, abs.DeltaSIMp95)
	}
	if abs.DeltaSIP95 <= 0 || abs.DeltaSIP99 <= abs.DeltaSIP95 {
		t.Fatalf("|ΔSI| thresholds P99=%v P95=%v", abs.DeltaSIP99, abs.DeltaSIP95)
	}
	// q99(|Δ|) is the two-sided 99% bound of a symmetric Δ: q99.5(Δ).
	if d := math.Abs(abs.DeltaSIP99 - signed.DeltaSIP99); d > 0.03 {
		t.Fatalf("|ΔSI| P99 %v vs signed P99 %v differ by %v", abs.DeltaSIP99, signed.DeltaSIP99, d)
	}
	p.highDepth = 51 // different depth; also check the cache key separates the two
	if twoBulkCacheKey(p) == twoBulkCacheKey(simParams{highDepth: 51, lowDepth: 50, p0High: 0.5, p0Low: 0.5,
		rep: 20000, highBulkSize: 20, lowBulkSize: 20, population: "F2"}) {
		t.Fatal("cache key ignores absDeltaSI")
	}
}

func TestAbsDeltaSIStatConfigs(t *testing.T) {
	var names []string
	for _, c := range absDeltaSIStatConfigs {
		names = append(names, c.name)
		if c.name == "DeltaSI+" || c.name == "DeltaSI-" || c.name == "HMMLOD" {
			t.Fatalf("track %q in the |ΔSI| set (no signed ΔSI or HMM without polarity)", c.name)
		}
	}
	if len(absDeltaSIStatConfigs) != len(twoBulkStatConfigs)-2 || names[1] != "AbsDeltaSI" {
		t.Fatalf("unexpected |ΔSI| track set %v", names)
	}
}

func TestConsolidateAbsDeltaSINotFolded(t *testing.T) {
	// 2b without --ref-parent: no HMM, so k = 2.
	sm := []SmoothedStats{{ZDeltaSI: -2, ZGstat: 0}}
	consolidate(sm, smoothLayout{highBulk: true, lowBulk: true}, nil)
	if want := -2 / math.Sqrt(2); math.Abs(sm[0].CompositeZ-want) > 1e-6 {
		t.Fatalf("unsigned CompositeZ = %v, want %v (Z(|ΔSI|) not folded)", sm[0].CompositeZ, want)
	}
	consolidate(sm, smoothLayout{highBulk: true, lowBulk: true, signedDSI: true, hmm: true}, nil)
	if want := 2 / math.Sqrt(2); math.Abs(sm[0].CompositeZ-want) > 1e-6 {
		t.Fatalf("signed CompositeZ = %v, want %v", sm[0].CompositeZ, want)
	}
}
