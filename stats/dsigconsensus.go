package stats

import (
	"fmt"
	"math"
	"slices"
	"sort"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// DsiGconsensusZ: one consensus score for the ΔSI and G evidence.
//
// Each component is standardised against the genome-wide null simulation at
// every marker, not against the genome's own spread (the old robust-Z
// CompositeZ divided by the genome's median and MAD, which made its scale
// depend on how much of the genome is null, put its genome-wide line near
// 20–25, and let it veto QTLs that the HMM LOD and BRM both called):
//
//   - ΔSI family (NZDeltaSI): z = signal/σ, BRM's standardised statistic
//     (smoothed signed ΔSI; one-bulk SmAFDev; 2b without --ref-parent the
//     polarity-free √ΔSI²c), divided by its root-mean-square over null
//     genomes. Either direction is evidence: v₁ = NZDeltaSI².
//   - G family (NZGstat): G′ (one-bulk G1) standardised by its null mean and
//     SD at the marker, v₂ = (G′ − E₀[G′])/SD₀[G′].
//
// Each vᵢ becomes a one-sided normal score sᵢ = Φ⁻¹(F₀ᵢ(vᵢ)), F₀ᵢ the pooled
// null distribution of vᵢ (empirical in the body, an exponential tail beyond
// its 99th percentile), so sᵢ ~ N(0, 1) under the null whatever the shape of
// the statistic (G′ is right-skewed; |ΔSI| folded; √ΔSI²c has mass at 0).
// NZGstat is s₂. The two are combined by Stouffer's method corrected for
// their null correlation r:
//
//	DsiGconsensusZ = (s₁ + s₂) / √(2 + 2r)
//
// The per-marker null moments of G, the root-mean-square of z, the two null
// distributions and r are estimated from their own null genomes
// (consensusMomentReps), separate from those that set the thresholds; the
// thresholds are genome-wide quantiles of the maximum DsiGconsensusZ over the
// threshold genomes (extremeQuantile), as for every other statistic.
//
// ΔSI and G′ are two views of the same reads, so r is high (about 0.9 on
// simulated F2 data) and the consensus carries little beyond ΔSI alone; r is
// reported so that is visible.

// ConsensusThresholds are DsiGconsensusZ's genome-wide thresholds and the
// null model that standardises it.
type ConsensusThresholds struct {
	P99, P95 float64
	// R is the null correlation of the two components' normal scores.
	R float64
	// Model standardises the components; nil when not computed.
	Model *ConsensusModel
}

// ConsensusModel is the null standardisation of DsiGconsensusZ's components.
type ConsensusModel struct {
	// DSIScale is the root-mean-square of z = signal/σ over the null genomes.
	DSIScale float64
	// Moments holds the G track's null mean and SD at the markers the null
	// scored, per chromosome (interpolated to other markers); Fallback is used
	// on a chromosome the null did not score.
	Moments  map[string]*GMoments
	Fallback GMoments
	// CDF1, CDF2 are the pooled null distributions of v₁ = NZDeltaSI² and of
	// the standardised G, v₂.
	CDF1, CDF2 NullCDF
	R          float64
	Reps       int // null genomes behind the model
}

// GMoments are a chromosome's G-track null means and SDs at positions Pos
// (sorted).
type GMoments struct {
	Pos      []int64
	Mean, SD []float64
}

// at returns the moments at pos: linear between the scored positions, the
// nearest end outside them.
func (g *GMoments) at(pos int64) (mean, sd float64) {
	n := len(g.Pos)
	if n == 0 {
		return 0, 0
	}
	i := sort.Search(n, func(k int) bool { return g.Pos[k] >= pos })
	switch {
	case i == 0:
		return g.Mean[0], g.SD[0]
	case i == n:
		return g.Mean[n-1], g.SD[n-1]
	case g.Pos[i] == pos:
		return g.Mean[i], g.SD[i]
	}
	t := float64(pos-g.Pos[i-1]) / float64(g.Pos[i]-g.Pos[i-1])
	return g.Mean[i-1] + t*(g.Mean[i]-g.Mean[i-1]), g.SD[i-1] + t*(g.SD[i]-g.SD[i-1])
}

// consensusSpec is a mode's two components at a smoothed marker.
type consensusSpec struct {
	brm  *brmModel
	rho  float64
	g    func(*SmoothedStats) float64
	dsiN string // component names for messages
	gN   string
}

func newConsensusSpec(cfg utils.AnalysisConfig, bsaType string, rho float64) (*consensusSpec, error) {
	brm, err := newBRMModel(cfg, bsaType)
	if err != nil {
		return nil, err
	}
	_, _, both, _ := BulkFlags(bsaType)
	c := &consensusSpec{brm: brm, rho: rho}
	switch {
	case both && cfg.SignedDeltaSI(bsaType):
		c.g, c.dsiN, c.gN = func(s *SmoothedStats) float64 { return s.SmGstat }, "ΔSI", "G′"
	case both:
		c.g, c.dsiN, c.gN = func(s *SmoothedStats) float64 { return s.SmGstat }, "|ΔSI|", "G′"
	default:
		c.g, c.dsiN, c.gN = func(s *SmoothedStats) float64 { return s.SmOneBulkG }, "AFDev", "G1"
	}
	return c, nil
}

// z is the ΔSI-family signal/σ at a marker (signed; ≥ 0 polarity-free).
func (c *consensusSpec) z(s *SmoothedStats) float64 {
	sd := math.Sqrt(c.brm.sigma2(s, c.rho))
	if !(sd > 0) {
		return 0
	}
	return c.brm.signal(s) / sd
}

// scores returns NZDeltaSI, NZGstat and DsiGconsensusZ at a marker whose G
// null moments are (mean, sd).
func (m *ConsensusModel) scores(c *consensusSpec, s *SmoothedStats, mean, sd float64) (nzd, nzg, cons float64) {
	if m.DSIScale > 0 {
		nzd = c.z(s) / m.DSIScale
	}
	s1 := m.CDF1.score(nzd * nzd)
	nzg = m.CDF2.score(standardise(c.g(s), mean, sd))
	return nzd, nzg, (s1 + nzg) / stoufferNorm(m.R)
}

func standardise(x, mean, sd float64) float64 {
	if !(sd > 0) {
		return 0
	}
	return (x - mean) / sd
}

// stoufferNorm is the SD of s₁ + s₂ for unit-variance scores with correlation r.
func stoufferNorm(r float64) float64 { return math.Sqrt(math.Max(2+2*r, 1e-9)) }

// applyDsiGConsensus fills NZDeltaSI, NZGstat and DsiGconsensusZ on the real
// data's smoothed markers from th's null model (zeros when there is none).
func applyDsiGConsensus(cfg utils.AnalysisConfig, bsaType string, sm []SmoothedStats, th Thresholds) error {
	m := th.Consensus.Model
	if m == nil {
		return nil
	}
	spec, err := newConsensusSpec(cfg, bsaType, th.Rho)
	if err != nil {
		return err
	}
	for i := range sm {
		s := &sm[i]
		mom := m.Moments[s.CHROM]
		if mom == nil {
			mom = &m.Fallback
		}
		mean, sd := mom.at(s.POS)
		s.NZDeltaSI, s.NZGstat, s.DsiGconsensusZ = m.scores(spec, s, mean, sd)
	}
	return nil
}

// ── fitting the model on null genomes ───────────────────────────────────────

// consensusMomentReps is the number of null genomes that estimate the
// model: half the threshold genomes, at least 50. They need no HMM, so they
// cost a fraction of a threshold genome.
func consensusMomentReps(reps int) int { return max(50, reps/2) }

// maxConsensusSamples caps the (z, G) pairs kept for the pooled null
// distributions and r.
const maxConsensusSamples = 2_000_000

// consensusFitter accumulates the moment genomes, in replicate order.
type consensusFitter struct {
	spec     *consensusSpec
	chrom    []string // scored markers of a null genome (the same in every one)
	pos      []int64
	gSum     []float64
	gSq      []float64
	z2       float64 // Σ z²
	nz       int
	stride   int // keep every stride-th (genome, marker) pair
	seen     int
	keepK    []int32 // marker index of each kept pair
	keepZ    []float64
	keepG    []float64
	reps     int
	expected int
}

func newConsensusFitter(spec *consensusSpec, reps int) *consensusFitter {
	return &consensusFitter{spec: spec, expected: reps}
}

// add takes one null genome's scored markers.
func (f *consensusFitter) add(sm []SmoothedStats) error {
	if f.reps == 0 {
		f.chrom, f.pos = make([]string, len(sm)), make([]int64, len(sm))
		for k := range sm {
			f.chrom[k], f.pos[k] = sm[k].CHROM, sm[k].POS
		}
		f.gSum, f.gSq = make([]float64, len(sm)), make([]float64, len(sm))
		f.stride = max(1, (f.expected*len(sm)+maxConsensusSamples-1)/maxConsensusSamples)
	} else if len(sm) != len(f.pos) {
		return fmt.Errorf("null genome scored %d markers, the first %d", len(sm), len(f.pos))
	}
	for k := range sm {
		g, z := f.spec.g(&sm[k]), f.spec.z(&sm[k])
		f.gSum[k] += g
		f.gSq[k] += g * g
		f.z2 += z * z
		f.nz++
		if f.seen%f.stride == 0 {
			f.keepK = append(f.keepK, int32(k))
			f.keepZ = append(f.keepZ, z)
			f.keepG = append(f.keepG, g)
		}
		f.seen++
	}
	f.reps++
	return nil
}

// model returns the fitted model and the per-marker G moments in null-scored
// order (index-aligned with the threshold genomes).
func (f *consensusFitter) model() (*ConsensusModel, []float64, []float64, error) {
	if f.reps < 2 || len(f.pos) == 0 {
		return nil, nil, nil, fmt.Errorf("DsiGconsensusZ: too few null genomes (%d) or scored markers (%d)", f.reps, len(f.pos))
	}
	n := float64(f.reps)
	mean, sd := make([]float64, len(f.pos)), make([]float64, len(f.pos))
	for k := range f.pos {
		mean[k] = f.gSum[k] / n
		sd[k] = math.Sqrt(math.Max(f.gSq[k]-n*mean[k]*mean[k], 0) / (n - 1))
	}
	m := &ConsensusModel{Moments: map[string]*GMoments{}, Reps: f.reps}
	if f.nz > 0 {
		m.DSIScale = math.Sqrt(f.z2 / float64(f.nz))
	}
	for k, c := range f.chrom {
		g := m.Moments[c]
		if g == nil {
			g = &GMoments{}
			m.Moments[c] = g
		}
		g.Pos, g.Mean, g.SD = append(g.Pos, f.pos[k]), append(g.Mean, mean[k]), append(g.SD, sd[k])
	}
	m.Fallback = GMoments{Pos: []int64{0}, Mean: []float64{median(mean)}, SD: []float64{median(sd)}}

	v1, v2 := make([]float64, len(f.keepK)), make([]float64, len(f.keepK))
	for i, k := range f.keepK {
		nzd := 0.0
		if m.DSIScale > 0 {
			nzd = f.keepZ[i] / m.DSIScale
		}
		v1[i], v2[i] = nzd*nzd, standardise(f.keepG[i], mean[k], sd[k])
	}
	m.CDF1, m.CDF2 = fitNullCDF(v1), fitNullCDF(v2)
	s1, s2 := make([]float64, len(v1)), make([]float64, len(v2))
	for i := range v1 {
		s1[i], s2[i] = m.CDF1.score(v1[i]), m.CDF2.score(v2[i])
	}
	m.R = pearson(s1, s2)
	return m, mean, sd, nil
}

func median(x []float64) float64 {
	if len(x) == 0 {
		return 0
	}
	s := slices.Clone(x)
	slices.Sort(s)
	if n := len(s); n%2 == 1 {
		return s[n/2]
	} else {
		return (s[n/2-1] + s[n/2]) / 2
	}
}

func pearson(a, b []float64) float64 {
	n := float64(len(a))
	if n < 2 {
		return 0
	}
	var ma, mb float64
	for i := range a {
		ma += a[i]
		mb += b[i]
	}
	ma, mb = ma/n, mb/n
	var sab, saa, sbb float64
	for i := range a {
		da, db := a[i]-ma, b[i]-mb
		sab += da * db
		saa += da * da
		sbb += db * db
	}
	if saa <= 0 || sbb <= 0 {
		return 0
	}
	return sab / math.Sqrt(saa*sbb)
}

// reportConsensusModel prints the model's r, with a note when the two
// components are near-copies.
func reportConsensusModel(spec *consensusSpec, m *ConsensusModel) {
	color.Cyan("DsiGconsensusZ: %s and %s standardised against %d null genomes; null correlation of the components r = %.3f",
		spec.dsiN, spec.gN, m.Reps, m.R)
	if m.R > 0.9 {
		color.Yellow("DsiGconsensusZ: r = %.2f, so %s and %s are near-copies here: the consensus carries little beyond %s alone, and their agreement is not independent confirmation",
			m.R, spec.dsiN, spec.gN, spec.dsiN)
	}
}

// ── pooled null distribution → normal score ─────────────────────────────────

// NullCDF is a pooled null distribution, stored compactly: quantiles at
// probabilities 0.001, 0.002, …, 0.990 (Knots), the sample minimum, and an
// exponential upper tail beyond the 0.99 quantile with mean excess Beta.
type NullCDF struct {
	Knots []float64
	Min   float64
	N     int
	Beta  float64
}

const (
	cdfStep  = 0.001
	cdfKnots = 990  // knot i is the quantile at (i+1)·cdfStep
	cdfTailP = 0.01 // probability beyond the last knot
)

// fitNullCDF summarises sample x (any order) as a NullCDF.
func fitNullCDF(x []float64) NullCDF {
	n := len(x)
	if n == 0 {
		return NullCDF{}
	}
	s := slices.Clone(x)
	slices.Sort(s)
	c := NullCDF{Knots: make([]float64, cdfKnots), Min: s[0], N: n}
	for i := range c.Knots { // type-7 quantile
		h := float64(n-1) * float64(i+1) * cdfStep
		lo := int(h)
		v := s[lo]
		if lo+1 < n {
			v += (h - float64(lo)) * (s[lo+1] - s[lo])
		}
		c.Knots[i] = v
	}
	t, sum, k := c.Knots[cdfKnots-1], 0.0, 0
	for i := n - 1; i >= 0 && s[i] > t; i-- {
		sum += s[i] - t
		k++
	}
	if k > 0 {
		c.Beta = sum / float64(k)
	}
	if !(c.Beta > 0) {
		c.Beta = 1e-12
	}
	return c
}

// score returns Φ⁻¹(F(v)), the normal score of v under the null
// distribution: N(0, 1) for null samples. Ties (e.g. a mass at 0) get the
// midpoint of their probability range.
func (c *NullCDF) score(v float64) float64 {
	if len(c.Knots) == 0 {
		return 0
	}
	last := c.Knots[cdfKnots-1]
	if v > last { // exponential tail: P(V > v) = 0.01·exp(−(v − last)/β)
		return normalUpperLog(-math.Log(cdfTailP) + (v-last)/c.Beta)
	}
	pMin := 0.5 / float64(max(c.N, 1))
	i := sort.SearchFloat64s(c.Knots, v) // first knot ≥ v
	var p float64
	switch {
	case c.Knots[i] == v:
		j := i
		for j+1 < cdfKnots && c.Knots[j+1] == v {
			j++
		}
		p = float64(i+j+2) / 2 * cdfStep
		if i == 0 && v <= c.Min { // a mass at the minimum reaches down to 0
			p = float64(j+1) / 2 * cdfStep
		}
	case i == 0:
		if v <= c.Min || c.Knots[0] <= c.Min {
			p = pMin
		} else {
			p = pMin + (v-c.Min)/(c.Knots[0]-c.Min)*(cdfStep-pMin)
		}
	default:
		lo, hi := c.Knots[i-1], c.Knots[i]
		p = (float64(i) + (v-lo)/(hi-lo)) * cdfStep
	}
	return normalQuantile(math.Min(math.Max(p, pMin), 1-cdfTailP))
}

// normalUpperLog returns z with P(Z > z) = exp(−L), accurate for any L > 0
// (beyond float64's smallest probability it uses the Mills-ratio expansion).
func normalUpperLog(L float64) float64 {
	if L < 700 {
		return -normalQuantile(math.Exp(-L))
	}
	// −ln P(Z > z) ≈ z²/2 + ln z + ln √(2π)
	z := math.Sqrt(2 * L)
	for range 3 {
		z = math.Sqrt(2 * (L - math.Log(z) - 0.5*math.Log(2*math.Pi)))
	}
	return z
}

// normalQuantile is Φ⁻¹(p) (Wichura's AS 241, PPND16: relative accuracy
// about 1e-16, including the far tails).
func normalQuantile(p float64) float64 {
	if p <= 0 {
		return math.Inf(-1)
	}
	if p >= 1 {
		return math.Inf(1)
	}
	q := p - 0.5
	if math.Abs(q) <= 0.425 {
		r := 0.180625 - q*q
		return q * (((((((2509.0809287301226727*r+33430.575583588128105)*r+67265.770927008700853)*r+45921.953931549871457)*r+
			13731.693765509461125)*r+1971.5909503065514427)*r+133.14166789178437745)*r + 3.387132872796366608) /
			(((((((5226.495278852854561*r+28729.085735721942674)*r+39307.89580009271061)*r+21213.794301586595867)*r+
				5394.1960214247511077)*r+687.1870074920579083)*r+42.313330701600911252)*r + 1)
	}
	r := p
	if q > 0 {
		r = 1 - p
	}
	r = math.Sqrt(-math.Log(r))
	var v float64
	if r <= 5 {
		r -= 1.6
		v = (((((((7.7454501427834140764e-4*r+0.0227238449892691845833)*r+0.24178072517745061177)*r+1.27045825245236838258)*r+
			3.64784832476320460504)*r+5.7694972214606914055)*r+4.6303378461565452959)*r + 1.42343711074968357734) /
			(((((((1.05075007164441684324e-9*r+5.475938084995344946e-4)*r+0.0151986665636164571966)*r+0.14810397642748007459)*r+
				0.68976733498510000455)*r+1.6763848301838038494)*r+2.05319162663775882187)*r + 1)
	} else {
		r -= 5
		v = (((((((2.01033439929228813265e-7*r+2.71155556874348757815e-5)*r+0.0012426609473880784386)*r+0.026532189526576123093)*r+
			0.29656057182850489123)*r+1.7848265399172913358)*r+5.4637849111641143699)*r + 6.6579046435011037772) /
			(((((((2.04426310338993978564e-15*r+1.4215117583164458887e-7)*r+1.8463183175100546818e-5)*r+7.868691311456132591e-4)*r+
				0.0148753612908506148525)*r+0.13692988092273580531)*r+0.59983220655588793769)*r + 1)
	}
	if q < 0 {
		return -v
	}
	return v
}

// consensusMaximum returns the largest DsiGconsensusZ over a null genome's
// scored markers, with the G moments index-aligned (mean[k], sd[k]).
func (m *ConsensusModel) consensusMaximum(c *consensusSpec, sm []SmoothedStats, mean, sd []float64) (float64, error) {
	if len(sm) != len(mean) {
		return 0, fmt.Errorf("null genome scored %d markers, the model %d", len(sm), len(mean))
	}
	best := math.Inf(-1)
	for k := range sm {
		_, _, v := m.scores(c, &sm[k], mean[k], sd[k])
		best = math.Max(best, v)
	}
	return best, nil
}
