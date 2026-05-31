package stats

import (
	"math"
	"testing"
)

func TestBetaBinomialOneBulkUsesObservedNullLikelihood(t *testing.T) {
	got := BetaBinomialOneBulk(8, 2, 0.25)
	wantLogNull := 8*math.Log(0.25) + 2*math.Log(0.75)
	wantLogAlt := logBeta(9, 3) - logBeta(1, 1)
	want := (wantLogAlt - wantLogNull) / math.Log(10)
	if math.Abs(got-want) > 1e-12 {
		t.Fatalf("BetaBinomialOneBulk() = %v, want %v", got, want)
	}
}

func TestExpectedNullAFF3(t *testing.T) {
	if got := ExpectedNullAF("F3", false); got != 0.5 {
		t.Fatalf("ExpectedNullAF(F3) = %v, want 0.5", got)
	}
}
