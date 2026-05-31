package utils

import "testing"

func TestSimulateAFF3IsSupported(t *testing.T) {
	got := SimulateAF("F3", 40, 2000)
	if got < 0.35 || got > 0.65 {
		t.Fatalf("SimulateAF(F3) = %v, want near 0.5", got)
	}
}
