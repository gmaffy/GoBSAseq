package stats

import "testing"

func TestWilsonCI(t *testing.T) {
	lo, hi := WilsonCI(50, 100, 1.96)
	if lo >= 0.5 || hi <= 0.5 {
		t.Fatalf("expected interval to contain observed proportion: lo=%v hi=%v", lo, hi)
	}
	if lo < 0 || hi > 1 {
		t.Fatalf("interval outside [0,1]: lo=%v hi=%v", lo, hi)
	}
}

func TestBenjaminiHochbergMonotone(t *testing.T) {
	q := BenjaminiHochberg([]float64{0.01, 0.03, 0.2, 0.04})
	if len(q) != 4 {
		t.Fatalf("expected 4 q-values, got %d", len(q))
	}
	for i, v := range q {
		if v < 0 || v > 1 {
			t.Fatalf("q[%d] outside [0,1]: %v", i, v)
		}
	}
	if q[0] > q[2] {
		t.Fatalf("smallest p-value should not have larger q-value than largest p-value: %v", q)
	}
}
