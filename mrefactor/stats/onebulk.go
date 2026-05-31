package stats

import "math"

// ExpectedNullAF returns the Mendelian null allele frequency for unlinked loci.
func ExpectedNullAF(population string, recurrentBC bool) float64 {
	switch population {
	case "F2", "F3":
		return 0.5
	case "RIL":
		return 0.5
	case "BC":
		if recurrentBC {
			return 0.75
		}
		return 0.25
	default:
		return 0.5
	}
}

// GStatisticOneBulk tests deviation from a 50:50 (or BC) split within one bulk.
func GStatisticOneBulk(alt, ref float64, nullAF float64) float64 {
	a := alt + 0.5
	r := ref + 0.5
	total := a + r
	if total == 0 {
		return 0
	}
	expAlt := total * nullAF
	expRef := total * (1 - nullAF)
	g := 0.0
	if a > 0 {
		g += a * math.Log(a/expAlt)
	}
	if r > 0 {
		g += r * math.Log(r/expRef)
	}
	return 2 * g
}

// LodOneBulk is log10 likelihood ratio vs nullAF.
func LodOneBulk(alt, ref float64, nullAF float64) float64 {
	total := alt + ref
	if total == 0 {
		return 0
	}
	const eps = 1e-10
	p := math.Max(eps, math.Min(1-eps, alt/total))
	p0 := nullAF
	return float64(total) * (p*math.Log(p/p0) + (1-p)*math.Log((1-p)/(1-p0))) / math.Log(10)
}

// BetaBinomialOneBulk compares estimated p vs nullAF.
func BetaBinomialOneBulk(alt, ref float64, nullAF float64) float64 {
	total := alt + ref
	if total == 0 {
		return 0
	}
	logNull := float64(total) * (nullAF*math.Log(nullAF) + (1-nullAF)*math.Log(1-nullAF))
	logAlt := logBeta(alt+1, ref+1) - logBeta(1, 1)
	return (logAlt - logNull) / math.Log(10)
}

// AbsDeviationFromNull returns |SI - expected|.
func AbsDeviationFromNull(si, nullAF float64) float64 {
	return math.Abs(si - nullAF)
}

// ED4OneBulk is |SI-null|^4.
func ED4OneBulk(si, nullAF float64) float64 {
	d := si - nullAF
	return d * d * d * d
}
