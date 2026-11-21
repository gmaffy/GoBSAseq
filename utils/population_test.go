package utils

import (
	"math"
	"testing"
)

// ParsePopulation reproduces the per-property formulas for every supported
// population, and rejects anything else.
func TestParsePopulation(t *testing.T) {
	bcScale := func(n int) float64 { q := math.Pow(0.5, float64(n)); return 2 * (1 - 0.5*q) / (1 - q) }
	bcF2Scale := func(n int) float64 { q := math.Pow(0.5, float64(n-1)); return (4 - q) / (3 - q) }
	for _, c := range []struct {
		name              string
		p0, scale, mapExp float64
	}{
		{"F2", 0.5, 2, 1}, {"f3", 0.5, 4.0 / 3, 1.5}, {" F4 ", 0.5, 8.0 / 7, 1.75}, {"RIL", 0.5, 1, 2},
		{"BC1H", 0.75, bcScale(1), 1}, {"BC1L", 0.25, bcScale(1), 1}, {"BC3H", 0.9375, bcScale(3), 3},
		{"BC1F2H", 0.75, bcF2Scale(1), 2}, {"bc2f2l", 0.125, bcF2Scale(2), 3},
	} {
		p, err := ParsePopulation(c.name)
		if err != nil {
			t.Fatalf("%q: %v", c.name, err)
		}
		if math.Abs(p.P0-c.p0) > 1e-12 || math.Abs(p.VarScale-c.scale) > 1e-12 || p.MapExpansion != c.mapExp {
			t.Errorf("%q: got %+v, want P0 %v VarScale %v MapExpansion %v", c.name, p, c.p0, c.scale, c.mapExp)
		}
		if af, err := ExpectedAF(c.name); err != nil || af != p.P0 {
			t.Errorf("%q: ExpectedAF %v, %v", c.name, af, err)
		}
	}
	for _, bad := range []string{"", "F5", "F_2", "DH", "BC", "BCH", "BC0H", "BC1", "BC1F2", "BC1X", "BC-1H", "BC01H", "BC1F2F2H", "BCxH"} {
		if p, err := ParsePopulation(bad); err == nil {
			t.Errorf("%q: expected an error, got %+v", bad, p)
		}
	}
}
