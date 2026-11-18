package utils

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// ExpectedAF returns the expected frequency of the HIGH-parent allele under H₀
// (no QTL) for a mapping population: 0.5 for F2/F3/F4/RIL; for BCn and BCnF2,
// 1 − 0.5^(n+1) when backcrossed to the high parent (…H) and 0.5^(n+1) when
// backcrossed to the low parent (…L). Selfing does not change allele frequency,
// so BCnF2 shares the BCn value.
func ExpectedAF(popStruc string) (float64, error) {
	p := strings.ToUpper(strings.TrimSpace(popStruc))

	switch p {
	case "F2", "F3", "F4", "RIL":
		return 0.5, nil
	}

	if strings.HasPrefix(p, "BC") {
		suffix := p[2:]
		var n int
		var err error
		var isHigh bool

		if strings.HasSuffix(suffix, "F2H") || strings.HasSuffix(suffix, "F2L") {
			isHigh = strings.HasSuffix(suffix, "F2H")
			n, err = strconv.Atoi(suffix[:len(suffix)-3])
		} else if strings.HasSuffix(suffix, "H") || strings.HasSuffix(suffix, "L") {
			isHigh = strings.HasSuffix(suffix, "H")
			n, err = strconv.Atoi(suffix[:len(suffix)-1])
		}

		if err == nil && n >= 1 {
			p0 := 1.0 - math.Pow(0.5, float64(n+1))
			if isHigh {
				return p0, nil
			}
			return 1.0 - p0, nil
		}
	}

	return 0, fmt.Errorf("unknown population structure: %s", popStruc)
}
