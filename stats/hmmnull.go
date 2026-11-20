package stats

import (
	"fmt"
	"math"
	"math/rand"
	"sort"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// Genome-wide null simulation for the HMM LOD and CompositeZ.
//
// Each replicate is a no-QTL genome at the real markers:
//   - each bulk's pool frequency evolves along every chromosome by the same
//     linked pool model the HMM uses (prior Binomial(U, p0), switch
//     probabilities from the gap in Morgans);
//   - reads are drawn at every marker's real per-bulk depth, with the
//     sequencing error and (one-bulk) the estimated mapping bias;
//   - per-site statistics are computed by the same rawCalc, then the replicate
//     goes through the same scorer as the real data (HMM LOD, kernel
//     smoothing, robust Z, CompositeZ).
//
// HMM LOD thresholds are the 99% and 95% quantiles of each replicate's
// genome-wide maximum (genome-wide error control). CompositeZ thresholds keep
// their pointwise meaning: the 99% and 95% quantiles of CompositeZ pooled
// over sites of all replicates (now from smoothed, linked, real-depth nulls
// instead of single-SNP draws).

// nullCompositeSample is about how many markers per replicate are smoothed
// and sampled for the CompositeZ null.
const nullCompositeSample = 5000

// genomeNullThresholds runs cfg.NullReps null replicates.
func genomeNullThresholds(cfg utils.AnalysisConfig, bsaType string, raw []BSAstats) (ZThresholds, HMMThresholds, error) {
	reps := cfg.NullReps
	if reps <= 0 {
		reps = 200
	}
	if len(raw) == 0 {
		return ZThresholds{}, HMMThresholds{}, fmt.Errorf("no sites for the null simulation")
	}
	sc, err := newScorer(cfg, bsaType)
	if err != nil {
		return ZThresholds{}, HMMThresholds{}, err
	}
	sim, err := newNullSim(cfg, bsaType, raw)
	if err != nil {
		return ZThresholds{}, HMMThresholds{}, err
	}

	color.Cyan("Genome-wide null: %d simulated no-QTL genomes at the real markers and depths\n", reps)
	bar := progressbar.Default(int64(reps), "Null genomes")
	rng := rand.New(rand.NewSource(nextSeed()))
	var maxima, cz []float64
	buf := make([]BSAstats, len(raw))
	every := max(1, len(raw)/nullCompositeSample)
	order := newMarkerOrder(raw) // replicates share the markers: sort once
	for range reps {
		sim.replicate(rng, buf)
		if sc.hmm != nil { // as RawStats: HMM LOD from the per-site counts, before smoothing
			computeHMMLOD(sc.hmm, buf, order)
		}
		mx := 0.0
		sm, err := sc.scoreOrdered(buf, order, false, every)
		if err != nil {
			return ZThresholds{}, HMMThresholds{}, err
		}
		for i := range buf { // genome-wide maximum over all markers
			mx = math.Max(mx, buf[i].HMMLOD)
		}
		maxima = append(maxima, mx)
		for i := range sm {
			cz = append(cz, sm[i].CompositeZ)
		}
		_ = bar.Add(1)
	}
	_ = bar.Finish()

	q := func(v []float64, p float64) float64 {
		sort.Float64s(v)
		return r6(stat.Quantile(p, stat.Empirical, v, nil))
	}
	zt := ZThresholds{ZP99: 2.576, ZP95: 1.960, CompositeZP99: q(cz, 0.99), CompositeZP95: q(cz, 0.95),
		NumStats: countZStats(bsaType)}
	var ht HMMThresholds
	if sc.lay.hmm {
		ht = HMMThresholds{LODP99: q(maxima, 0.99), LODP95: q(maxima, 0.95)}
		color.Green("HMM LOD genome-wide thresholds: p99 %.2f, p95 %.2f\n", ht.LODP99, ht.LODP95)
	}
	color.Green("CompositeZ null thresholds: p99 %.3f, p95 %.3f\n", zt.CompositeZP99, zt.CompositeZP95)
	return zt, ht, nil
}

// nullSim holds what a replicate needs: the template markers in chromosome
// order, the gaps, the pool models and the per-site calculator.
type nullSim struct {
	raw        []BSAstats
	chroms     [][]int     // indices into raw, by position
	dM         [][]float64 // gap (Morgans) before each marker, per chromosome
	high, low  *poolModel
	calc       rawCalc
	p0         float64
	oneBulk    bool
	sites      []filter.Site
	jHigh, jLo []int
}

func newNullSim(cfg utils.AnalysisConfig, bsaType string, raw []BSAstats) (*nullSim, error) {
	hb, lb, both, one := BulkFlags(bsaType)
	highN, lowN, oneN, err := cfg.BulkSizes(bsaType)
	if err != nil {
		return nil, err
	}
	p0, err := ExpectedAF(cfg.Population)
	if err != nil {
		return nil, err
	}
	s := &nullSim{raw: raw, p0: p0, oneBulk: one, sites: make([]filter.Site, len(raw))}
	mk := func(n int) (*poolModel, error) {
		m, err := newPoolModel(cfg.Population, n)
		return &m, err
	}
	switch {
	case both:
		if s.high, err = mk(highN); err == nil {
			s.low, err = mk(lowN)
		}
	case hb:
		s.high, err = mk(oneN)
	case lb:
		s.low, err = mk(oneN)
	}
	if err != nil {
		return nil, err
	}

	// Per-site calculator as in RawStats; the one-bulk p0 by polarity is read
	// back from the real sites.
	roles := siteRoles{hp: -1, lp: -1, hb: -1, lb: -1}
	if hb {
		roles.hb = 0
	}
	if lb {
		roles.lb = 1
	}
	s.calc = rawCalc{roles: roles, signedDSI: cfg.SignedDeltaSI(bsaType), oneBulk: one, p0: [2]float64{p0, p0}}
	for i := range raw {
		k := 0
		if raw[i].HighIsAlt {
			k = 1
		}
		s.calc.p0[k] = raw[i].OneBulkP0
	}
	if !one {
		s.calc.p0 = [2]float64{p0, p0}
	}

	// Chromosomes and gaps, exactly as the HMM computes them.
	byChrom := map[string][]int{}
	var order []string
	for i := range raw {
		c := raw[i].CHROM
		if _, ok := byChrom[c]; !ok {
			order = append(order, c)
		}
		byChrom[c] = append(byChrom[c], i)
	}
	for _, c := range order {
		idx := byChrom[c]
		sort.SliceStable(idx, func(a, b int) bool { return raw[idx[a]].POS < raw[idx[b]].POS })
		s.chroms = append(s.chroms, idx)
		s.dM = append(s.dM, gapsMorgans(raw, idx, cfg.CMPerMb))
	}
	return s, nil
}

// gapsMorgans returns the gap in Morgans before each marker of a chromosome
// (idx sorted by position): cmPerMb, or 1 Morgan over the marker span if 0.
func gapsMorgans(raw []BSAstats, idx []int, cmPerMb float64) []float64 {
	n := len(idx)
	mPerBp := cmPerMb / 100 / 1e6
	if cmPerMb <= 0 {
		mPerBp = 0
		if span := float64(raw[idx[n-1]].POS - raw[idx[0]].POS); span > 0 {
			mPerBp = 1 / span
		}
	}
	dM := make([]float64, n)
	for x := 1; x < n; x++ {
		dM[x] = float64(raw[idx[x]].POS-raw[idx[x-1]].POS) * mPerBp
	}
	return dM
}

// replicate fills out with one simulated no-QTL genome's per-site statistics.
func (s *nullSim) replicate(rng *rand.Rand, out []BSAstats) {
	for c, idx := range s.chroms {
		var jh, jl int
		for x, i := range idx {
			if s.high != nil {
				jh = s.high.stepState(rng, jh, x, s.dM[c][x])
			}
			if s.low != nil {
				jl = s.low.stepState(rng, jl, x, s.dM[c][x])
			}
			t := &s.raw[i]
			site := &s.sites[i]
			*site = filter.Site{Chrom: t.CHROM, Pos: t.POS, Ref: t.REF, Alt: t.ALT, HighIsAlt: t.HighIsAlt, NSamples: 2}
			shift := 0.0
			if s.oneBulk {
				shift = t.OneBulkP0 - s.p0
			}
			if s.high != nil {
				d := t.HighDepth
				site.Samples[0] = simReads(rng, s.high, jh, d, shift, t.HighIsAlt)
			}
			if s.low != nil {
				d := t.LowDepth
				site.Samples[1] = simReads(rng, s.low, jl, d, shift, t.HighIsAlt)
			}
			out[i] = s.calc.site(site)
		}
	}
}

// stepState draws the hidden state at marker x: from the prior at the first
// marker, otherwise by per-unit switches across the gap.
func (m *poolModel) stepState(rng *rand.Rand, j, x int, dM float64) int {
	if x == 0 {
		return binomSample(rng, m.U, m.pi)
	}
	a, b := m.switchProbs(dM)
	return j - binomSample(rng, j, a) + binomSample(rng, m.U-j, b)
}

// simReads draws a bulk's reads at depth d from pool state j.
func simReads(rng *rand.Rand, m *poolModel, j, d int, shift float64, highIsAlt bool) filter.SampleAD {
	f := m.freq[j]
	p := f*(1-m.err) + (1-f)*m.err + shift
	p = math.Min(math.Max(p, hmmMinProb), 1-hmmMinProb)
	h := int32(binomSample(rng, d, p))
	sa := filter.SampleAD{Ref: h, Alt: int32(d) - h, HasAD: true}
	if highIsAlt {
		sa.Ref, sa.Alt = sa.Alt, sa.Ref
	}
	return sa
}

// binomSample draws Binomial(n, p): inversion for small n·p (the common case
// for switch counts and modest depths), gonum otherwise.
func binomSample(rng *rand.Rand, n int, p float64) int {
	if n <= 0 || p <= 0 {
		return 0
	}
	if p >= 1 {
		return n
	}
	if p > 0.5 {
		return n - binomSample(rng, n, 1-p)
	}
	if float64(n)*p < 30 {
		u := rng.Float64()
		pk := math.Pow(1-p, float64(n))
		cdf, ratio := pk, p/(1-p)
		k := 0
		for u > cdf && k < n {
			pk *= float64(n-k) / float64(k+1) * ratio
			k++
			cdf += pk
		}
		return k
	}
	return int(distuv.Binomial{N: float64(n), P: p, Src: rng}.Rand())
}
