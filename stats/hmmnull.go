package stats

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"
	"sync/atomic"
	"unsafe"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// Genome-wide null simulation: the source of every threshold (Thresholds).
//
// Each replicate is a no-QTL genome at the real markers:
//   - each bulk's pool frequency evolves along every chromosome by the same
//     linked pool model the HMM uses (prior Binomial(U, p0), switch
//     probabilities from the gap in Morgans);
//   - reads are drawn at every marker's real per-bulk depth, with the
//     sequencing error, (one-bulk) the estimated mapping bias, and the
//     read-level overdispersion ρ estimated from the data (estimateRho), so
//     overdispersed data are not scored against a binomial null;
//   - per-site statistics are computed by the same rawCalc and the HMM LOD as
//     in RawStats, then the replicate goes through the same kernel smoothing
//     as the real data.
//
// Two sets of replicates: consensusMomentReps genomes (smoothing only, no
// HMM) fit DsiGconsensusZ's null standardisation (dsigconsensus.go); NullReps
// further genomes set every threshold, DsiGconsensusZ's included.
//
// Every threshold is genome-wide: a quantile of the most extreme value each
// scanned statistic reaches in a replicate (maximum, and minimum for lower
// tails), from a generalized extreme-value fit to the replicates'
// extremes (extremeQuantile). The fit uses every replicate, so the 99% quantile is
// stable at the default 200 replicates, where the empirical one rests on the
// top two.
//
// Replicates run in parallel (nullPool), each from its own seeds, and are
// consumed in replicate order, so the thresholds do not depend on the number
// of cores.

// nullSeed seeds the null simulation, so thresholds are reproducible.
const nullSeed = 20240611

// nullExtreme collects one smoothed track's extremes over the replicates and
// turns them into thresholds.
type nullExtreme struct {
	get       func(*SmoothedStats) float64
	twoTailed bool
	set       func(t *Thresholds, p99, p95, m99, m95 float64)
	max, min  []float64
}

// nullExtremes lists the smoothed tracks that get thresholds in layout lay.
// (DsiGconsensusZ and BRM are handled separately: they need the null model.)
func nullExtremes(lay smoothLayout) []*nullExtreme {
	var out []*nullExtreme
	add := func(two bool, get func(*SmoothedStats) float64, set func(t *Thresholds, p99, p95, m99, m95 float64)) {
		out = append(out, &nullExtreme{get: get, twoTailed: two, set: set})
	}
	if lay.highBulk {
		add(true, func(s *SmoothedStats) float64 { return s.SmHighSI }, func(t *Thresholds, p99, p95, m99, m95 float64) {
			t.TwoBulk.HighSIP99, t.TwoBulk.HighSIP95, t.TwoBulk.HighSIMp99, t.TwoBulk.HighSIMp95 = p99, p95, m99, m95
		})
	}
	if lay.lowBulk {
		add(true, func(s *SmoothedStats) float64 { return s.SmLowSI }, func(t *Thresholds, p99, p95, m99, m95 float64) {
			t.TwoBulk.LowSIP99, t.TwoBulk.LowSIP95, t.TwoBulk.LowSIMp99, t.TwoBulk.LowSIMp95 = p99, p95, m99, m95
		})
	}
	if lay.both() {
		add(lay.signedDSI, func(s *SmoothedStats) float64 { return s.SmDeltaSI }, func(t *Thresholds, p99, p95, m99, m95 float64) {
			t.TwoBulk.DeltaSIP99, t.TwoBulk.DeltaSIP95, t.TwoBulk.DeltaSIMp99, t.TwoBulk.DeltaSIMp95 = p99, p95, m99, m95
		})
		add(false, func(s *SmoothedStats) float64 { return s.SmGstat }, func(t *Thresholds, p99, p95, _, _ float64) {
			t.TwoBulk.GstatP99, t.TwoBulk.GstatP95 = p99, p95
		})
		add(false, func(s *SmoothedStats) float64 { return s.SmED }, func(t *Thresholds, p99, p95, _, _ float64) {
			t.TwoBulk.ED4P99, t.TwoBulk.ED4P95 = p99, p95
		})
	}
	if lay.one() {
		add(true, func(s *SmoothedStats) float64 { return s.SmAFDev }, func(t *Thresholds, p99, p95, m99, m95 float64) {
			t.OneBulk.AFDevP99, t.OneBulk.AFDevP95, t.OneBulk.AFDevMp99, t.OneBulk.AFDevMp95 = p99, p95, m99, m95
		})
		add(false, func(s *SmoothedStats) float64 { return s.SmOneBulkG }, func(t *Thresholds, p99, p95, _, _ float64) {
			t.OneBulk.OneBulkGstatP99, t.OneBulk.OneBulkGstatP95 = p99, p95
		})
	}
	return out
}

// observe records the track's extremes over one replicate's scored markers.
func (e *nullExtreme) observe(sm []SmoothedStats) {
	mx, mn := math.Inf(-1), math.Inf(1)
	for i := range sm {
		v := e.get(&sm[i])
		mx, mn = math.Max(mx, v), math.Min(mn, v)
	}
	e.max, e.min = append(e.max, mx), append(e.min, -mn)
}

// apply sets the track's thresholds: one-sided at the 99%/95% quantiles of
// the maxima; two-tailed at 99.5%/97.5% of the maxima and of the negated
// minima (NaN lower tail for one-sided tracks).
func (e *nullExtreme) apply(t *Thresholds) {
	if !e.twoTailed {
		e.set(t, r6(extremeQuantile(e.max, 0.99)), r6(extremeQuantile(e.max, 0.95)), math.NaN(), math.NaN())
		return
	}
	e.set(t, r6(extremeQuantile(e.max, 0.995)), r6(extremeQuantile(e.max, 0.975)),
		r6(-extremeQuantile(e.min, 0.995)), r6(-extremeQuantile(e.min, 0.975)))
}

// genomeNullThresholds runs the genome-wide null simulation and returns every
// threshold and DsiGconsensusZ's null model.
func genomeNullThresholds(cfg utils.AnalysisConfig, bsaType string, raw []BSAstats) (Thresholds, error) {
	reps := cfg.NullReps
	if reps <= 0 {
		reps = 200
	}
	if len(raw) == 0 {
		return Thresholds{}, fmt.Errorf("no sites for the null simulation")
	}
	sc, err := newScorer(cfg, bsaType)
	if err != nil {
		return Thresholds{}, err
	}
	sim, err := newNullSim(cfg, bsaType, raw)
	if err != nil {
		return Thresholds{}, err
	}
	order := newMarkerOrder(raw) // replicates share the markers: sort once
	rho, nObs := estimateRho(raw, order, sim)
	sim.rho = rho
	if nObs < minRhoObs {
		color.Yellow("Read overdispersion: too few informative markers (%d, need %d); the null uses binomial reads (ρ = 0)", nObs, minRhoObs)
	} else {
		color.Cyan("Read overdispersion ρ = %.4f (beta-binomial, from %d observations); used in the null simulation", rho, nObs)
		if rho >= maxRho {
			color.Yellow("Read overdispersion is at its cap (%.2f): check depths, bulk sizes and the population type", maxRho)
		}
	}
	every := nullEvery(order, raw, sc.bw)
	ext := nullExtremes(sc.lay)
	alpha := cfg.BrmAlpha
	var brm *brmModel // BRM's genome-wide threshold, when --brm-alpha is valid
	if alpha > 0 && alpha < 1 {
		if brm, err = newBRMModel(cfg, bsaType); err != nil {
			return Thresholds{}, err
		}
	}
	spec, err := newConsensusSpec(cfg, bsaType, rho)
	if err != nil {
		return Thresholds{}, err
	}
	pool := newNullPool(sim, sc, order, every)
	color.Cyan("Null simulation: %d replicate(s) at a time, %d core(s) each\n", pool.outer, pool.inner)

	// DsiGconsensusZ's null model, from its own genomes (no HMM needed).
	mReps := consensusMomentReps(reps)
	color.Cyan("DsiGconsensusZ null model: %d simulated no-QTL genomes\n", mReps)
	fit := newConsensusFitter(spec, mReps)
	bar := progressbar.Default(int64(mReps), "Null model")
	if err := pool.run(mReps, nullSeed+1, false, func(_ int, sm []SmoothedStats, _ float64) error {
		_ = bar.Add(1)
		return fit.add(sm)
	}); err != nil {
		return Thresholds{}, err
	}
	_ = bar.Finish()
	model, gMean, gSD, err := fit.model()
	if err != nil {
		return Thresholds{}, err
	}
	reportConsensusModel(spec, model)

	color.Cyan("Genome-wide null: %d simulated no-QTL genomes at the real markers and depths\n", reps)
	bar = progressbar.Default(int64(reps), "Null genomes")
	hmmMax, brmMax, consMax := make([]float64, 0, reps), make([]float64, 0, reps), make([]float64, 0, reps)
	err = pool.run(reps, nullSeed, sc.hmm != nil, func(_ int, sm []SmoothedStats, lod float64) error {
		if len(sm) == 0 {
			return fmt.Errorf("no markers were scored in the null simulation")
		}
		hmmMax = append(hmmMax, lod)
		for _, e := range ext {
			e.observe(sm)
		}
		if brm != nil {
			mz := 0.0
			for i := range sm {
				mz = math.Max(mz, brm.z(&sm[i], rho))
			}
			brmMax = append(brmMax, mz)
		}
		cz, err := model.consensusMaximum(spec, sm, gMean, gSD)
		if err != nil {
			return err
		}
		consMax = append(consMax, cz)
		_ = bar.Add(1)
		return nil
	})
	if err != nil {
		return Thresholds{}, err
	}
	_ = bar.Finish()

	th := Thresholds{Rho: rho, NullReps: reps}
	if brm != nil {
		th.BRMU, th.BRMAlpha = r6(extremeQuantile(brmMax, 1-alpha)), alpha
		color.Green("BRM genome-wide threshold: |z| ≥ %.3f (α = %g)\n", th.BRMU, alpha)
	}
	for _, e := range ext {
		e.apply(&th)
	}
	if sc.lay.hmm {
		th.HMM = HMMThresholds{LODP99: r6(extremeQuantile(hmmMax, 0.99)), LODP95: r6(extremeQuantile(hmmMax, 0.95))}
		color.Green("HMM LOD genome-wide thresholds: p99 %.2f, p95 %.2f\n", th.HMM.LODP99, th.HMM.LODP95)
	}
	th.Consensus = ConsensusThresholds{P99: r6(extremeQuantile(consMax, 0.99)), P95: r6(extremeQuantile(consMax, 0.95)),
		R: r6(model.R), Model: model}
	color.Green("DsiGconsensusZ genome-wide thresholds: p99 %.3f, p95 %.3f\n", th.Consensus.P99, th.Consensus.P95)
	return th, nil
}

// nullMemBudget bounds the memory of the replicates the null runs at once.
const nullMemBudget = 2 << 30

// smoothBytesPerMarker bounds the smoothing state per marker of the
// chromosome being smoothed (positions, weights, values, query results).
const smoothBytesPerMarker = 256

// nullOuterOverride, if > 0, fixes how many replicates run at once (tests).
var nullOuterOverride int

// nullPool runs null replicates in parallel: outer replicates at a time, each
// on inner goroutines (outer × inner ≈ the cores). Whole replicates in
// parallel keep every core busy, where one replicate at a time waited on its
// longest chromosome three times over (simulation, HMM, smoothing); outer is
// bounded by nullMemBudget, since each replicate holds a genome buffer.
type nullPool struct {
	sim          *nullSim
	sc           *scorer
	order        *markerOrder
	every        int
	outer, inner int
}

func newNullPool(sim *nullSim, sc *scorer, order *markerOrder, every int) *nullPool {
	procs := max(runtime.GOMAXPROCS(0), 1)
	maxN := 0
	for _, idx := range order.chroms {
		maxN = max(maxN, len(idx))
	}
	per := int64(len(sim.raw))*int64(unsafe.Sizeof(BSAstats{})) + // replicate buffer
		int64(maxN)*smoothBytesPerMarker + // smoothing state of the largest chromosome
		int64(len(sim.raw)/max(every, 1)+1)*int64(unsafe.Sizeof(SmoothedStats{})) // scored markers
	if sc.hmm != nil {
		per += sc.hmm.chromBytes(maxN)
	}
	outer := int(min(int64(procs), max(1, nullMemBudget/max(per, 1))))
	if nullOuterOverride > 0 {
		outer = nullOuterOverride
	}
	return &nullPool{sim: sim, sc: sc, order: order, every: every, outer: outer, inner: max(1, procs/outer)}
}

// run simulates n replicates from seed and passes each one's scored markers
// (and, with hmm, its maximum HMM LOD) to consume, in replicate order, on the
// calling goroutine. Replicate r's per-chromosome seeds are the r-th block of
// a stream seeded with seed, whatever the scheduling, so the results do not
// depend on the number of workers. The first error stops the run.
func (p *nullPool) run(n int, seed int64, hmm bool, consume func(rep int, sm []SmoothedStats, lod float64) error) error {
	nc := len(p.sim.chroms)
	rng := rand.New(rand.NewSource(seed))
	seeds := make([]int64, n*nc)
	for i := range seeds {
		seeds[i] = rng.Int63()
	}
	type result struct {
		rep int
		sm  []SmoothedStats
		lod float64
		err error
	}
	results := make(chan result, p.outer)
	// A worker takes a token before a replicate number and the collector
	// returns it once that replicate is consumed: at most 2·outer replicates
	// are held at once, and the lowest unconsumed one always holds a token.
	tokens := make(chan struct{}, 2*p.outer)
	var next atomic.Int64
	var stop atomic.Bool
	var wg sync.WaitGroup
	for range p.outer {
		wg.Add(1)
		go func() {
			defer wg.Done()
			w := p.sim.newWorker()
			for !stop.Load() {
				tokens <- struct{}{}
				r := int(next.Add(1)) - 1
				if r >= n {
					<-tokens
					return
				}
				p.sim.replicateWith(w, seeds[r*nc:(r+1)*nc], p.inner)
				res := result{rep: r}
				if hmm {
					res.lod = p.sc.hmm.maxHMMLODWorkers(w.buf, p.order, p.inner)
				}
				res.sm, res.err = p.sc.scoreWorkers(w.buf, p.order, false, p.every, p.inner)
				results <- res
			}
		}()
	}
	go func() { wg.Wait(); close(results) }()

	pending := map[int]result{}
	want := 0
	var firstErr error
	for res := range results {
		pending[res.rep] = res
		for {
			r, ok := pending[want]
			if !ok {
				break
			}
			delete(pending, want)
			want++
			if firstErr == nil {
				if firstErr = r.err; firstErr == nil {
					firstErr = consume(r.rep, r.sm, r.lod)
				}
				if firstErr != nil {
					stop.Store(true)
				}
			}
			<-tokens
		}
	}
	return firstErr
}

// nullEvery is how many markers apart the null replicates are scored: enough
// that the scored markers are on average ≤ σ/10 apart, so the replicate's
// track maxima are not underestimated (the smoothed tracks vary on the
// scale of σ), and no more, since only the maximum is kept.
func nullEvery(order *markerOrder, raw []BSAstats, sigma float64) int {
	span, n := 0.0, 0
	for _, idx := range order.chroms {
		span += float64(raw[idx[len(idx)-1]].POS - raw[idx[0]].POS)
		n += len(idx)
	}
	if n == 0 || span <= 0 {
		return 1
	}
	return max(1, int(sigma/10/(span/float64(n))))
}

// extremeQuantile fits a generalized extreme-value (GEV) distribution to
// sample maxima x by probability-weighted moments (Hosking, Wallis & Wood
// 1985) and returns its p quantile. The GEV is the limit law of maxima; its
// shape lets the tail be heavier than Gumbel's (maxima of strongly linked
// smoothed tracks) or lighter, and it reduces to the Gumbel when
// the shape is 0. PWM estimates are stable for small samples, so the 99%
// quantile is usable at 200 replicates, where the empirical one rests on the
// top two. Falls back to the empirical quantile for small or constant x.
func extremeQuantile(x []float64, p float64) float64 {
	n := len(x)
	s := append([]float64(nil), x...)
	sort.Float64s(s)
	if n < 20 || s[0] == s[n-1] {
		return stat.Quantile(p, stat.Empirical, s, nil)
	}
	var b0, b1, b2 float64
	fn := float64(n)
	for j, v := range s { // j = rank − 1
		fj := float64(j)
		b0 += v
		b1 += fj / (fn - 1) * v
		b2 += fj * (fj - 1) / ((fn - 1) * (fn - 2)) * v
	}
	b0, b1, b2 = b0/fn, b1/fn, b2/fn
	l2 := 2*b1 - b0 // second L-moment
	c := l2/(3*b2-b0) - math.Ln2/math.Log(3)
	k := 7.8590*c + 2.9554*c*c // Hosking's shape (k < 0: heavy tail)
	y := -math.Log(p)          // reduced variate: x_p = ξ + α(1 − y^k)/k
	if math.Abs(k) < 1e-6 {    // Gumbel
		alpha := l2 / math.Ln2
		return b0 - 0.5772156649015329*alpha - alpha*math.Log(y)
	}
	g := math.Gamma(1 + k)
	alpha := l2 * k / (g * (1 - math.Pow(2, -k)))
	xi := b0 + alpha*(g-1)/k
	return xi + alpha*(1-math.Pow(y, k))/k
}

// Overdispersion estimation (estimateRho).
const (
	minRhoDepth = 10    // every bulk depth of an observation
	maxRhoGapM  = 0.001 // neighbours at most 0.1 cM apart share the pool frequency
	minRhoObs   = 500
	maxRho      = 0.25
	chi2Median1 = 0.454936423119572 // median of χ²₁
)

// estimateRho estimates the read-level overdispersion ρ, the beta-binomial
// intraclass correlation: Var(reads | p) = d·p(1−p)·(1 + (d−1)ρ).
//
// With known polarity it uses neighbouring markers ≤ 0.1 cM apart (both
// representative of their position, every bulk depth ≥ 10). They share the
// pool frequency, so the difference of their signals (ΔSI, or SI − P0 in
// one-bulk modes) has no pool or QTL component, only read noise, with
// variance V(ρ) = Σ q·(1 + (d−1)ρ)/d over the bulk observations, q being the
// pool model's E[p(1−p)]. ρ solves median(D²/V(ρ)) = median(χ²₁), which
// ignores outliers. Without polarity (2b, --ref-parent none) neighbours'
// signals may have opposite signs, so single sites are used and V adds the
// pool model's frequency variance; any excess pool variance then shows up as ρ.
//
// It returns ρ in [0, maxRho] and the number of observations; ρ = 0 when
// there are fewer than minRhoObs.
func estimateRho(raw []BSAstats, order *markerOrder, s *nullSim) (float64, int) {
	type bulkMom struct{ q, vp float64 }
	mom := func(m *poolModel) (bm bulkMom) {
		var mean, mean2 float64
		for j, w := range m.prior {
			p := m.freq[j]*(1-m.err) + (1-m.freq[j])*m.err
			bm.q += w * p * (1 - p)
			mean += w * p
			mean2 += w * p * p
		}
		bm.vp = mean2 - mean*mean
		return bm
	}
	var hi, lo bulkMom
	if s.high != nil {
		hi = mom(s.high)
	}
	if s.low != nil {
		lo = mom(s.low)
	}
	// site returns a marker's signal x and its variance terms: read variance
	// a + ρ·b, pool variance vp.
	site := func(i int) (x, a, b, vp float64, ok bool) {
		t := &raw[i]
		add := func(m bulkMom, d int) bool {
			if d < minRhoDepth {
				return false
			}
			a += m.q / float64(d)
			b += m.q * float64(d-1) / float64(d)
			vp += m.vp
			return true
		}
		if s.high != nil && !add(hi, t.HighDepth) {
			return
		}
		if s.low != nil && !add(lo, t.LowDepth) {
			return
		}
		x = t.HighSI - t.LowSI
		if s.oneBulk {
			x = t.OneBulkAFDev
		}
		return x, a, b, vp, true
	}

	var d2, va, vb []float64
	if s.oneBulk || s.calc.signedDSI {
		for c, idx := range s.chroms {
			prev, gap := -1, 0.0
			for k, i := range idx {
				if k > 0 {
					gap += s.dM[c][k]
				}
				if !order.rep[i] {
					continue
				}
				if prev >= 0 && gap <= maxRhoGapM && raw[i].POS != raw[prev].POS {
					x1, a1, b1, _, ok1 := site(i)
					x0, a0, b0, _, ok0 := site(prev)
					if ok0 && ok1 {
						d := x1 - x0
						d2, va, vb = append(d2, d*d), append(va, a0+a1), append(vb, b0+b1)
					}
				}
				prev, gap = i, 0
			}
		}
	} else {
		for i := range raw {
			if !order.rep[i] {
				continue
			}
			if x, a, b, vp, ok := site(i); ok {
				d2, va, vb = append(d2, x*x), append(va, a+vp), append(vb, b)
			}
		}
	}
	n := len(d2)
	if n < minRhoObs {
		return 0, n
	}
	// below(ρ) counts D²/V(ρ) ≤ median(χ²₁); it increases with ρ, and ρ is
	// where it reaches half the observations.
	below := func(rho float64) int {
		c := 0
		for k := range d2 {
			if d2[k] <= chi2Median1*(va[k]+rho*vb[k]) {
				c++
			}
		}
		return c
	}
	half := (n + 1) / 2
	if below(0) >= half {
		return 0, n
	}
	if below(maxRho) < half {
		return maxRho, n
	}
	l, h := 0.0, maxRho
	for range 40 {
		m := (l + h) / 2
		if below(m) >= half {
			h = m
		} else {
			l = m
		}
	}
	return h, n
}

// nullSim holds what a replicate needs: the template markers in chromosome
// order, the gaps, the pool models, the per-site calculator and the
// overdispersion. It is read-only while replicates run; each worker has its
// own simWorker (genome buffer and RNG streams).
type nullSim struct {
	raw       []BSAstats
	chroms    [][]int     // indices into raw, by position
	dM        [][]float64 // gap (Morgans) before each marker, per chromosome
	high, low *poolModel
	calc      rawCalc
	p0        float64
	rho       float64 // read overdispersion (estimateRho)
	oneBulk   bool
	bySize    []int // chromosome indices, most markers first (scheduling)
}

// simWorker is one worker's replicate state: the genome buffer and one RNG
// stream per chromosome, reseeded for every replicate.
type simWorker struct {
	buf  []BSAstats
	srcs []splitMix
	rngs []*rand.Rand // wrap srcs
}

func (s *nullSim) newWorker() *simWorker {
	w := &simWorker{buf: make([]BSAstats, len(s.raw)), srcs: make([]splitMix, len(s.chroms)), rngs: make([]*rand.Rand, len(s.chroms))}
	for c := range w.srcs {
		w.rngs[c] = rand.New(&w.srcs[c])
	}
	return w
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
	s := &nullSim{raw: raw, p0: p0, oneBulk: one}
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
	s.bySize = make([]int, len(s.chroms))
	for c := range s.chroms {
		s.bySize[c] = c
	}
	sort.SliceStable(s.bySize, func(a, b int) bool { return len(s.chroms[s.bySize[a]]) > len(s.chroms[s.bySize[b]]) })
	return s, nil
}

// splitMix is SplitMix64, a small rand.Source64 that is cheap to reseed: one
// per chromosome per replicate, so replicates are reproducible whatever the
// scheduling.
type splitMix struct{ s uint64 }

func (r *splitMix) Uint64() uint64 {
	r.s += 0x9e3779b97f4a7c15
	z := r.s
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}
func (r *splitMix) Int63() int64    { return int64(r.Uint64() >> 1) }
func (r *splitMix) Seed(seed int64) { r.s = uint64(seed) }

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

// replicateWith fills w.buf with one simulated no-QTL genome's per-site
// statistics. Chromosome c is drawn from its own stream seeded with seeds[c],
// so the result does not depend on how many goroutines (inner) share the
// chromosomes.
func (s *nullSim) replicateWith(w *simWorker, seeds []int64, inner int) {
	for c := range w.srcs {
		w.srcs[c].Seed(seeds[c])
	}
	if inner <= 1 {
		for _, c := range s.bySize {
			s.replicateChrom(c, w.buf, w.rngs[c])
		}
		return
	}
	var next atomic.Int64
	var wg sync.WaitGroup
	for range min(inner, len(s.chroms)) {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for {
				k := int(next.Add(1)) - 1
				if k >= len(s.bySize) {
					return
				}
				c := s.bySize[k]
				s.replicateChrom(c, w.buf, w.rngs[c])
			}
		}()
	}
	wg.Wait()
}

func (s *nullSim) replicateChrom(c int, out []BSAstats, rng *rand.Rand) {
	var jh, jl int
	for x, i := range s.chroms[c] {
		if s.high != nil {
			jh = s.high.stepState(rng, jh, x, s.dM[c][x])
		}
		if s.low != nil {
			jl = s.low.stepState(rng, jl, x, s.dM[c][x])
		}
		t := &s.raw[i]
		site := &filter.Site{Chrom: t.CHROM, Pos: t.POS, Ref: t.REF, Alt: t.ALT, HighIsAlt: t.HighIsAlt, NSamples: 2}
		shift := 0.0
		if s.oneBulk {
			shift = t.OneBulkP0 - s.p0
		}
		if s.high != nil {
			site.Samples[0] = simReads(rng, s.high, jh, t.HighDepth, shift, s.rho, t.HighIsAlt)
		}
		if s.low != nil {
			site.Samples[1] = simReads(rng, s.low, jl, t.LowDepth, shift, s.rho, t.HighIsAlt)
		}
		out[i] = s.calc.site(site)
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

// simReads draws a bulk's reads at depth d from pool state j. With rho > 0 the
// read probability gets a per-site perturbation of variance ρ·p(1−p) (normal,
// clamped), which gives the beta-binomial read variance
// d·p(1−p)·(1 + (d−1)ρ) without the cost of beta draws.
func simReads(rng *rand.Rand, m *poolModel, j, d int, shift, rho float64, highIsAlt bool) filter.SampleAD {
	f := m.freq[j]
	p := f*(1-m.err) + (1-f)*m.err + shift
	p = math.Min(math.Max(p, hmmMinProb), 1-hmmMinProb)
	if rho > 0 {
		p += math.Sqrt(rho*p*(1-p)) * rng.NormFloat64()
		p = math.Min(math.Max(p, hmmMinProb), 1-hmmMinProb)
	}
	h := int32(binomSample(rng, d, p))
	sa := filter.SampleAD{Ref: h, Alt: int32(d) - h, HasAD: true}
	if highIsAlt {
		sa.Ref, sa.Alt = sa.Alt, sa.Ref
	}
	return sa
}

// binomSample draws Binomial(n, p). For n·p < binomInvMax (read counts at
// usual depths, switch counts) it sums inversion draws over chunks with
// n·p < 30 each: exact, no allocation, cost ∝ n·p. Deeper sites use gonum.
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
	np := float64(n) * p
	if np >= binomInvMax {
		return int(distuv.Binomial{N: float64(n), P: p, Src: rng}.Rand())
	}
	chunk := n
	if np >= 30 {
		chunk = int(30 / p)
	}
	k := 0
	for n > 0 {
		m := min(n, chunk)
		k += binomInversion(rng, m, p)
		n -= m
	}
	return k
}

// binomInvMax bounds n·p for the chunked inversion in binomSample.
const binomInvMax = 300

// binomInversion draws Binomial(n, p) by inversion; n·p < 30, p ≤ 0.5.
func binomInversion(rng *rand.Rand, n int, p float64) int {
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
