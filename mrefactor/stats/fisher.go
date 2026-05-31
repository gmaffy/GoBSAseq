package stats

import "math"

// FisherExact2x2 returns a two-sided Fisher exact p-value for the contingency table:
//
//	| a | b |
//	| c | d |
//
// Implementation follows the hypergeometric sum used in PyBSASeq pre-filtering.
func FisherExact2x2(a, b, c, d int) float64 {
	if a < 0 || b < 0 || c < 0 || d < 0 {
		return 1
	}
	row1 := a + b
	row2 := c + d
	col1 := a + c
	n := row1 + row2
	if n == 0 {
		return 1
	}

	obs := hypergeomPMF(a, row1, col1, n)
	p := obs
	// two-sided: sum probabilities of tables as or more extreme
	minA := maxInt(0, col1-row2)
	maxA := minInt(row1, col1)
	for i := minA; i <= maxA; i++ {
		if i == a {
			continue
		}
		pi := hypergeomPMF(i, row1, col1, n)
		if pi <= obs+1e-15 {
			p += pi
		}
	}
	if p > 1 {
		return 1
	}
	return p
}

func hypergeomPMF(k, rowSum, colSum, n int) float64 {
	return math.Exp(logHypergeomPMF(k, rowSum, colSum, n))
}

func logHypergeomPMF(k, rowSum, colSum, n int) float64 {
	return logChoose(colSum, k) + logChoose(n-colSum, rowSum-k) - logChoose(n, rowSum)
}

func logChoose(n, k int) float64 {
	if k < 0 || k > n {
		return math.Inf(-1)
	}
	if k > n-k {
		k = n - k
	}
	r := 0.0
	for i := 1; i <= k; i++ {
		r += math.Log(float64(n-k+i)) - math.Log(float64(i))
	}
	return r
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
