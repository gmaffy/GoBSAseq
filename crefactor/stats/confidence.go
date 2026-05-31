package stats

import (
	"math"
	"sort"
)

// WilsonCI returns a two-sided Wilson score interval for a binomial proportion.
// It is more stable than the Wald interval for low coverage and near fixation.
func WilsonCI(successes, trials int, z float64) (float64, float64) {
	if trials <= 0 {
		return 0, 0
	}
	n := float64(trials)
	p := float64(successes) / n
	z2 := z * z
	den := 1 + z2/n
	center := (p + z2/(2*n)) / den
	margin := z * math.Sqrt((p*(1-p)+z2/(4*n))/n) / den
	lo := math.Max(0, center-margin)
	hi := math.Min(1, center+margin)
	return lo, hi
}

// BenjaminiHochberg converts p-values to monotone BH q-values.
func BenjaminiHochberg(pvalues []float64) []float64 {
	type point struct {
		idx int
		p   float64
	}
	points := make([]point, 0, len(pvalues))
	for i, p := range pvalues {
		if math.IsNaN(p) || p < 0 {
			p = 1
		}
		if p > 1 {
			p = 1
		}
		points = append(points, point{idx: i, p: p})
	}
	sort.Slice(points, func(i, j int) bool { return points[i].p < points[j].p })

	qvalues := make([]float64, len(pvalues))
	minQ := 1.0
	m := float64(len(points))
	for i := len(points) - 1; i >= 0; i-- {
		rank := float64(i + 1)
		q := points[i].p * m / rank
		if q < minQ {
			minQ = q
		}
		if minQ > 1 {
			minQ = 1
		}
		qvalues[points[i].idx] = minQ
	}
	return qvalues
}
