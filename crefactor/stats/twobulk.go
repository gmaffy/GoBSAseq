package stats

import "math"

// SNPIndex is ALT / (REF + ALT) allele depth fraction.
func SNPIndex(altDepth, refDepth int) float64 {
	total := altDepth + refDepth
	if total == 0 {
		return 0
	}
	return float64(altDepth) / float64(total)
}

// GStatisticTwoBulk is the pooled G-test on a 2×2 allele-count table (Magwene et al.).
func GStatisticTwoBulk(highAlt, highRef, lowAlt, lowRef int) float64 {
	hba := float64(highAlt) + 0.5
	hbr := float64(highRef) + 0.5
	lba := float64(lowAlt) + 0.5
	lbr := float64(lowRef) + 0.5

	highTotal := hba + hbr
	lowTotal := lba + lbr
	total := highTotal + lowTotal
	if total == 0 {
		return 0
	}

	expHighAlt := highTotal * (hba + lba) / total
	expHighRef := highTotal * (hbr + lbr) / total
	expLowAlt := lowTotal * (hba + lba) / total
	expLowRef := lowTotal * (hbr + lbr) / total

	g := 0.0
	if hba > 0 {
		g += hba * math.Log(hba/expHighAlt)
	}
	if hbr > 0 {
		g += hbr * math.Log(hbr/expHighRef)
	}
	if lba > 0 {
		g += lba * math.Log(lba/expLowAlt)
	}
	if lbr > 0 {
		g += lbr * math.Log(lbr/expLowRef)
	}
	return 2 * g
}

// EuclideanDistance4 is (Δ allele frequency)^4 (Magwene-style ED).
func EuclideanDistance4(highSI, lowSI float64) float64 {
	d := highSI - lowSI
	return d * d * d * d
}

func logBeta(a, b float64) float64 {
	la, _ := math.Lgamma(a)
	lb, _ := math.Lgamma(b)
	lab, _ := math.Lgamma(a + b)
	return la + lb - lab
}

// LODTwoBulk compares separate bulk frequencies vs a pooled null.
func LODTwoBulk(ref1, alt1, ref2, alt2 int) float64 {
	n1 := float64(ref1 + alt1)
	n2 := float64(ref2 + alt2)
	total := n1 + n2
	if n1 == 0 || n2 == 0 || total == 0 {
		return 0
	}

	const eps = 1e-10
	clamp := func(p float64) float64 {
		if p <= 0 {
			return eps
		}
		if p >= 1 {
			return 1 - eps
		}
		return p
	}

	p1 := clamp(float64(alt1) / n1)
	p2 := clamp(float64(alt2) / n2)
	p0 := clamp(float64(alt1+alt2) / total)

	logLik := func(k, n, p float64) float64 {
		if n == 0 {
			return 0
		}
		return k*math.Log(p) + (n-k)*math.Log(1-p)
	}

	logL1 := logLik(float64(alt1), n1, p1) + logLik(float64(alt2), n2, p2)
	logL0 := logLik(float64(alt1), n1, p0) + logLik(float64(alt2), n2, p0)
	return (logL1 - logL0) / math.Log(10)
}

// BetaBinomialLogBFTwoBulk supports separate vs shared allele frequencies across bulks.
func BetaBinomialLogBFTwoBulk(highSucc, highFail, lowSucc, lowFail int) float64 {
	alphaH, betaH := 0.5, 0.5
	alphaL, betaL := 0.5, 0.5

	logAlt := logBeta(alphaH+float64(highSucc), betaH+float64(highFail)) - logBeta(alphaH, betaH)
	logAlt += logBeta(alphaL+float64(lowSucc), betaL+float64(lowFail)) - logBeta(alphaL, betaL)

	alpha0, beta0 := 0.5, 0.5
	logNull := logBeta(alpha0+float64(highSucc+lowSucc), beta0+float64(highFail+lowFail)) - logBeta(alpha0, beta0)
	return (logAlt - logNull) / math.Log(10)
}
