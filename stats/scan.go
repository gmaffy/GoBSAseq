package stats

import "math"

// crossing is one region of consecutive markers at or beyond a threshold
// line: markers [first, last], ends interpolated where the margin crosses 0,
// and the marker with the largest strength.
type crossing struct {
	first, last int
	start, end  float64
	peak        int
}

// scanCrossings finds the regions of n markers (sorted by position) whose
// margin — how far beyond the line the marker is — is ≥ 0. A region's start
// is linearly interpolated between the last marker below the line and the
// first at or beyond it (the first marker's position at a chromosome start),
// its end likewise; its peak is the marker with the largest strength. A NaN
// margin counts as below the line. It is the one threshold-crossing scanner
// for QTL detection (FindPeakIntersections) and BRM (scanBRMBlocks).
func scanCrossings(n int, pos func(int) float64, margin, strength func(int) float64) []crossing {
	var out []crossing
	var cur crossing
	in := false
	prev := 0.0
	at := func(i int, m0, m1 float64) float64 { // margins m0 < 0 at i−1, m1 ≥ 0 at i (or the reverse)
		x0, x1 := pos(i-1), pos(i)
		return x0 + m0/(m0-m1)*(x1-x0)
	}
	for i := range n {
		m := margin(i)
		beyond := m >= 0 // false for NaN
		switch {
		case beyond && !in:
			in = true
			cur = crossing{first: i, last: i, start: pos(i), peak: i}
			if i > 0 && !math.IsInf(prev, -1) {
				cur.start = at(i, prev, m)
			}
		case beyond:
			cur.last = i
			if strength(i) > strength(cur.peak) {
				cur.peak = i
			}
		case in:
			cur.end = pos(i - 1)
			if !math.IsNaN(m) && !math.IsInf(m, -1) {
				cur.end = at(i, prev, m)
			}
			out = append(out, cur)
			in = false
		}
		if math.IsNaN(m) {
			m = math.Inf(-1)
		}
		prev = m
	}
	if in {
		cur.end = pos(n - 1)
		out = append(out, cur)
	}
	return out
}
