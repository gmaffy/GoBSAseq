package utils

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Population is a mapping population's null (no-QTL) properties. It is
// parsed once from the --population name (ParsePopulation), so every part of
// the analysis uses the same values.
type Population struct {
	// Name is the canonical (upper-case) name.
	Name string
	// P0 is the expected frequency of the high-parent allele: 0.5 for
	// F2/F3/F4/RIL; for BCn and BCnF2, 1 − 0.5^(n+1) when backcrossed to the
	// high parent (…H) and 0.5^(n+1) when backcrossed to the low parent (…L).
	// Selfing does not change allele frequency, so BCnF2 shares the BCn value.
	P0 float64
	// VarScale is c in the pool frequency variance P0(1 − P0)/(c·n) of n
	// plants (genotype dosage): F2 2, F3 4/3, F4 8/7, RIL 1;
	// BCn 2(1 − q/2)/(1 − q) with q = 0.5^n; BCnF2 (4 − q′)/(3 − q′) with
	// q′ = 0.5^(n−1).
	VarScale float64
	// MapExpansion is the number of meioses' worth of recombination the
	// population's chromosomes carry relative to one F2 meiosis: F2 1, F3 1.5,
	// F4 1.75, RIL 2, BCn n, BCnF2 n + 1.
	MapExpansion float64
}

// ParsePopulation parses a population name (case-insensitive): F2, F3, F4,
// RIL, BCnH, BCnL, BCnF2H, BCnF2L with n ≥ 1. Unknown names are an error.
func ParsePopulation(name string) (Population, error) {
	p := strings.ToUpper(strings.TrimSpace(name))
	switch p {
	case "F2":
		return Population{p, 0.5, 2, 1}, nil
	case "F3":
		return Population{p, 0.5, 4.0 / 3, 1.5}, nil
	case "F4":
		return Population{p, 0.5, 8.0 / 7, 1.75}, nil
	case "RIL":
		return Population{p, 0.5, 1, 2}, nil
	}
	if s, ok := strings.CutPrefix(p, "BC"); ok && len(s) >= 2 {
		high := strings.HasSuffix(s, "H")
		if !high && !strings.HasSuffix(s, "L") {
			return Population{}, unknownPopulation(name)
		}
		s = s[:len(s)-1]
		s, f2 := strings.CutSuffix(s, "F2")
		n, err := strconv.Atoi(s)
		if err != nil || n < 1 || s != strconv.Itoa(n) {
			return Population{}, unknownPopulation(name)
		}
		p0 := 1 - math.Pow(0.5, float64(n+1))
		if !high {
			p0 = 1 - p0
		}
		if f2 {
			q := math.Pow(0.5, float64(n-1))
			return Population{p, p0, (4 - q) / (3 - q), float64(n + 1)}, nil
		}
		q := math.Pow(0.5, float64(n))
		return Population{p, p0, 2 * (1 - 0.5*q) / (1 - q), float64(n)}, nil
	}
	return Population{}, unknownPopulation(name)
}

func unknownPopulation(name string) error {
	return fmt.Errorf("unknown population structure %q (supported: F2, F3, F4, RIL, BCnH, BCnL, BCnF2H, BCnF2L with n ≥ 1)", name)
}

// ExpectedAF returns the population's P0 (ParsePopulation).
func ExpectedAF(name string) (float64, error) {
	pop, err := ParsePopulation(name)
	return pop.P0, err
}
