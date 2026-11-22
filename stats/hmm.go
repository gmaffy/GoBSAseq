package stats

import (
	"fmt"
	"math"
	"runtime"
	"sync"

	"github.com/gmaffy/GoBSAseq/utils"
)

// HMM-based LOD (after MULTIPOOL, Edwards & Gifford 2012).
//
// Each bulk is a pool of plants. At every SNP a hidden number j of "units"
// carry the high-parent allele, giving a pool frequency f = j/U. Reads are
// binomial draws from f. Between neighbouring SNPs j changes only through
// crossovers, so the hidden frequencies form a Markov chain along the
// chromosome and every linked SNP informs the same pool frequency.
//
// Pool model. A bulk of n plants is modelled as U = round(c·n) exchangeable
// units, c being the population variance scale (utils.Population.VarScale).
// Each unit carries the high allele with probability π = the population's
// expected high-allele frequency, so the prior of j is Binomial(U, π). This
// matches the mean p0 and the pool-sampling variance p0(1−p0)/(c·n) of every
// supported population exactly (F2: U = 2n chromosomes; RIL: U = n plants of
// two identical copies), and the distribution shape exactly for F2 and RIL.
//
// Linkage. A unit switches from high to low between SNPs with probability
// a = 2(1−π)·r and from low to high with b = 2π·r, where r is the per-meiosis
// Haldane recombination fraction of the gap times the population's map
// expansion (F2 1, F3 1.5, F4 1.75, RIL 2, BCn n, BCnF2 n+1). This keeps π
// stationary and gives the pool frequency the correct linkage
// autocorrelation, 1 − a − b = 1 − 2r (exact for F2).
//
// Reads. P(high-allele read) = f(1−e) + (1−f)e + shift, with sequencing error
// e and, in one-bulk modes, shift = the estimated reference-mapping bias at
// the site (OneBulkP0 − p0). Reads are binomial (no extra overdispersion).
//
// LOD at SNP x, from one forward–backward pass per bulk:
//
//	L(j) = P(all reads on the chromosome | j_x = j) ∝ α_x(j)·β_x(j) / prior(j)
//
//	one-bulk: LOD = log10 max_j L(j) − log10 Σ_j prior(j)·L(j)
//	          (pool frequency free vs. no QTL: j_x from its prior around p0)
//	two-bulk: LOD = log10 max_j L_H(j) + log10 max_j L_L(j)
//	          − log10 max_f Σ_j Bin(j; U_H, f)·L_H(j) · Σ_j Bin(j; U_L, f)·L_L(j)
//	          (bulks free vs. both bulks drawn around one shared frequency f,
//	          which is itself free, so segregation distortion that shifts
//	          both bulks equally is not scored as a QTL)
//
// The two-bulk denominator is maximised over a grid of f values with a
// parabolic refinement; real data and the null simulation use the same code,
// so thresholds stay consistent.
//
// Markers within hmmBinMorgans of the first marker of a bin form one HMM step
// (reads combined); every marker gets its bin's LOD.
//
// The HMM needs the high-parent allele at every SNP, so it is not computed in
// 2b mode without --ref-parent (a phase-free joint two-bulk model was tried
// and is kept in oldcode/hmmjoint.go.backup; it was ~20× slower).

const (
	hmmSeqError = 0.01   // per-read probability of reading the other allele
	hmmMinProb  = 1e-6   // clamp for read probabilities
	hmmTiny     = 1e-300 // floor for log of likelihood sums
	hmmFGrid    = 41     // grid points for the shared frequency f (two-bulk)
	// hmmBinMorgans: markers closer than this (0.05 cM) share one HMM step;
	// their reads are combined. A crossover inside a bin is rare (≈ U·5e-4 per
	// bin), so this barely changes the LOD but makes dense marker sets cheap.
	// Against 0.01 cM bins on 1 kb-spaced simulated F2 data: the HMM ran 3–3.4×
	// faster, the genome maximum LOD and the peak LOD of every simulated QTL were
	// unchanged (to 5e-4), and 1.5-LOD support intervals moved by ≤ 10 kb for a
	// strong QTL and by ≤ 0.8 Mb on the flat 3.5–17 Mb plateaus of moderate ones.
	hmmBinMorgansDefault = 5e-4
)

// hmmBinMorgans is the bin width in use (a variable so tests can compare
// against unbinned runs).
var hmmBinMorgans = hmmBinMorgansDefault

// poolModel is one bulk's hidden-pool model.
type poolModel struct {
	U     int
	pi    float64   // stationary high-allele probability per unit (p0)
	exp   float64   // map expansion factor
	prior []float64 // Binomial(U, pi), j = 0..U
	freq  []float64 // j/U
	err   float64   // sequencing error rate in the read model
}

func newPoolModel(pop string, bulkSize int) (poolModel, error) {
	if bulkSize <= 0 {
		return poolModel{}, fmt.Errorf("bulk size must be > 0, got %d", bulkSize)
	}
	p, err := utils.ParsePopulation(pop)
	if err != nil {
		return poolModel{}, err
	}
	U := max(1, int(math.Round(p.VarScale*float64(bulkSize))))
	m := poolModel{U: U, pi: p.P0, exp: p.MapExpansion, prior: binomPMF(U, p.P0), freq: make([]float64, U+1), err: hmmSeqError}
	for j := range m.freq {
		m.freq[j] = float64(j) / float64(U)
	}
	return m, nil
}

// binomPMF returns Binomial(n, p) probabilities for 0..n.
func binomPMF(n int, p float64) []float64 {
	out := make([]float64, n+1)
	if p <= 0 {
		out[0] = 1
		return out
	}
	if p >= 1 {
		out[n] = 1
		return out
	}
	lp, lq := math.Log(p), math.Log1p(-p)
	for k := 0; k <= n; k++ {
		lc, _ := math.Lgamma(float64(n + 1))
		lk, _ := math.Lgamma(float64(k + 1))
		lnk, _ := math.Lgamma(float64(n - k + 1))
		out[k] = math.Exp(lc - lk - lnk + float64(k)*lp + float64(n-k)*lq)
	}
	return out
}

// switchProbs returns the per-unit switch probabilities (high→low a,
// low→high b) across a gap of dM Morgans.
func (m *poolModel) switchProbs(dM float64) (a, b float64) {
	if dM <= 0 {
		return 0, 0
	}
	r := 0.5 * (1 - math.Exp(-2*m.exp*dM))
	return math.Min(2*(1-m.pi)*r, 1), math.Min(2*m.pi*r, 1)
}

// stepTables holds one gap's per-unit switch probabilities and scratch for
// the binomial loss and gain distributions: powA[j] = (1−a)^j and
// powB[k] = (1−b)^k start them.
type stepTables struct {
	a, b, ra, rb     float64
	powA, powB       []float64
	lossBuf, gainBuf []float64
}

func (m *poolModel) tables(t *stepTables, a, b float64) {
	if len(t.powA) != m.U+1 {
		t.powA, t.powB = make([]float64, m.U+1), make([]float64, m.U+1)
		t.lossBuf, t.gainBuf = make([]float64, m.U+1), make([]float64, m.U+1)
	} else if t.a == a && t.b == b { // same gap as the previous step (evenly spaced bins)
		return
	}
	t.a, t.b, t.ra, t.rb = a, b, a/(1-a), b/(1-b)
	t.powA[0], t.powB[0] = 1, 1
	for k := 1; k <= m.U; k++ {
		t.powA[k] = t.powA[k-1] * (1 - a)
		t.powB[k] = t.powB[k-1] * (1 - b)
	}
}

// Truncation of the binomial loss/gain distributions: stop past the mean once
// terms fall below binomTail. Not raised to save time: on linked simulated
// data 1e-10 changed the LOD by < 3e-7 for a 1.0–1.2× speed-up, but where
// neighbouring markers disagree strongly (outlier sites; unlinked test data)
// the likelihood rests on these tiny many-switch transitions, and 1e-10
// moved the LOD by up to 0.23.
const binomTail = 1e-15

// binomInto fills buf[0:k] with P(Binomial(n, p) = i), i < k, truncated past
// the mean once terms fall below binomTail, and returns k. pow0 = (1−p)^n,
// ratio = p/(1−p).
func binomInto(buf []float64, n int, p, ratio, pow0 float64) int {
	if p == 0 {
		buf[0] = 1
		return 1
	}
	pk, mean := pow0, float64(n)*p
	for i := 0; i <= n; i++ {
		buf[i] = pk
		pk *= float64(n-i) / float64(i+1) * ratio
		if float64(i) > mean && pk < binomTail {
			return i + 1
		}
	}
	return n + 1
}

// forwardStep sets dst = Tᵀ src: the distribution of j at the next marker
// given src at this one. Exact per-unit model: from state j, l ~ Bin(j, a)
// high units switch to low and, independently, g ~ Bin(U − j, b) low units
// switch to high; j' = j − l + g (as the null simulation draws it). The
// scratch argument is unused (kept for the callers' shape).
func (m *poolModel) forwardStep(dst, src, _ []float64, t *stepTables) {
	if t.a == 0 && t.b == 0 {
		copy(dst, src)
		return
	}
	clear(dst)
	for j, v := range src {
		if v == 0 {
			continue
		}
		nl := binomInto(t.lossBuf, j, t.a, t.ra, t.powA[j])
		ng := binomInto(t.gainBuf, m.U-j, t.b, t.rb, t.powB[m.U-j])
		gain := t.gainBuf[:ng]
		for l := range nl {
			base := v * t.lossBuf[l]
			row := dst[j-l : j-l+ng]
			for g, pg := range gain {
				row[g] += base * pg
			}
		}
	}
}

// backwardStep sets dst(j) = Σ_j' T(j→j') src(j') for the same transition.
func (m *poolModel) backwardStep(dst, src, _ []float64, t *stepTables) {
	if t.a == 0 && t.b == 0 {
		copy(dst, src)
		return
	}
	for j := range dst {
		nl := binomInto(t.lossBuf, j, t.a, t.ra, t.powA[j])
		ng := binomInto(t.gainBuf, m.U-j, t.b, t.rb, t.powB[m.U-j])
		gain := t.gainBuf[:ng]
		s := 0.0
		for l := range nl {
			row := src[j-l : j-l+ng]
			inner := 0.0
			for g, pg := range gain {
				inner += pg * row[g]
			}
			s += t.lossBuf[l] * inner
		}
		dst[j] = s
	}
}

// logProbs holds log P(high-allele read | f_j) and log P(other read | f_j)
// for one read-probability shift. A run has very few distinct shifts (0 in
// two-bulk modes; theory ± mapping bias in one-bulk modes), so they are
// computed once per shift instead of once per marker.
type logProbs struct {
	shift        float64
	logP, log1mP []float64
}

// logCache returns the logProbs for shift, computing them on first use. Not
// safe for concurrent use: each forward–backward run has its own cache.
type logCache []logProbs

func (m *poolModel) logsFor(c *logCache, shift float64) *logProbs {
	for i := range *c {
		if (*c)[i].shift == shift {
			return &(*c)[i]
		}
	}
	lp := logProbs{shift: shift, logP: make([]float64, m.U+1), log1mP: make([]float64, m.U+1)}
	for j, f := range m.freq {
		p := f*(1-m.err) + (1-f)*m.err + shift
		p = math.Min(math.Max(p, hmmMinProb), 1-hmmMinProb)
		lp.logP[j], lp.log1mP[j] = math.Log(p), math.Log1p(-p)
	}
	*c = append(*c, lp)
	return &(*c)[len(*c)-1]
}

// binEmission fills e(j) ∝ P(reads of markers lo..hi−1 | f_j), scaled to
// max 1. cache holds the per-shift log probabilities.
func (m *poolModel) binEmission(e []float64, data hmmBulkData, lo, hi int, cache *logCache) {
	clear(e)
	for x := lo; x < hi; x++ {
		h, d := data.h[x], data.d[x]
		if d == 0 {
			continue
		}
		shift := 0.0
		if data.shift != nil {
			shift = data.shift[x]
		}
		lp := m.logsFor(cache, shift)
		fh, fl := float64(h), float64(d-h)
		for j := range e {
			e[j] += fh*lp.logP[j] + fl*lp.log1mP[j]
		}
	}
	mx := math.Inf(-1)
	for _, v := range e {
		mx = math.Max(mx, v)
	}
	for j := range e {
		e[j] = math.Exp(e[j] - mx)
	}
}

// hmmBulkData is one bulk's reads along one chromosome.
type hmmBulkData struct {
	h, d  []int     // high-allele reads and depth per marker
	shift []float64 // read-probability shift per marker (mapping bias), or nil
}

// chromLikelihoods runs forward–backward for one bulk along one chromosome
// and returns, per bin b, L_b(j) ∝ P(all reads | j_b = j), normalised to
// max 1. Bin b covers markers binStart[b]..binStart[b+1]−1 (binStart has a
// final sentinel = number of markers); dM[b] is the gap in Morgans from bin
// b−1 to bin b.
func (m *poolModel) chromLikelihoods(data hmmBulkData, binStart []int, dM []float64) [][]float32 {
	nb, S := len(binStart)-1, m.U+1
	alpha := make([][]float32, nb)
	emis := make([][]float32, nb)
	cur, prev, e, scratch := make([]float64, S), make([]float64, S), make([]float64, S), make([]float64, S)
	var tab stepTables
	var cache logCache
	// Per-bin switch probabilities, computed once in the forward pass and
	// reused by the backward pass.
	ab := make([]float64, 2*nb)

	for x := range nb {
		if x == 0 {
			copy(cur, m.prior)
		} else {
			a, b := m.switchProbs(dM[x])
			ab[2*x], ab[2*x+1] = a, b
			m.tables(&tab, a, b)
			m.forwardStep(cur, prev, scratch, &tab)
		}
		m.binEmission(e, data, binStart[x], binStart[x+1], &cache)
		er := make([]float32, S)
		sum := 0.0
		for j := range cur {
			er[j] = float32(e[j])
			cur[j] *= e[j]
			sum += cur[j]
		}
		emis[x] = er
		if sum <= 0 {
			sum = hmmTiny
		}
		row := make([]float32, S)
		for j := range cur {
			cur[j] /= sum
			row[j] = float32(cur[j])
		}
		alpha[x] = row
		prev, cur = cur, prev
	}

	// Backward pass; overwrite alpha with L = α·β/prior, normalised to max 1.
	beta, next := make([]float64, S), make([]float64, S)
	for j := range beta {
		beta[j] = 1
	}
	for x := nb - 1; x >= 0; x-- {
		if x < nb-1 {
			for j := range next {
				next[j] *= float64(emis[x+1][j])
			}
			m.tables(&tab, ab[2*(x+1)], ab[2*(x+1)+1])
			m.backwardStep(beta, next, scratch, &tab)
			mx := 0.0
			for _, v := range beta {
				mx = math.Max(mx, v)
			}
			if mx > 0 {
				for j := range beta {
					beta[j] /= mx
				}
			}
		}
		row := alpha[x]
		mx := 0.0
		for j := range row {
			v := 0.0
			if m.prior[j] > 0 {
				v = float64(row[j]) * beta[j] / m.prior[j]
			}
			scratch[j] = v
			mx = math.Max(mx, v)
		}
		for j := range row {
			if mx > 0 {
				row[j] = float32(scratch[j] / mx)
			} else {
				row[j] = 0
			}
		}
		copy(next, beta)
	}
	return alpha
}

// markerBins returns one bin per marker (for callers without binning).
func markerBins(n int) []int {
	b := make([]int, n+1)
	for i := range b {
		b[i] = i
	}
	return b
}

// lodOneBulk: log10 max_j L(j) − log10 Σ prior(j)·L(j); L has max 1.
func (m *poolModel) lodOneBulk(L []float32) float64 {
	s := 0.0
	for j, v := range L {
		s += m.prior[j] * float64(v)
	}
	return -math.Log10(math.Max(s, hmmTiny))
}

// twoBulkLOD holds the f-grid tables for the two-bulk denominator.
type twoBulkLOD struct {
	fGrid      []float64
	tabH, tabL [][]float64 // Bin(j; U, f) per grid point
}

func newTwoBulkLOD(mh, ml *poolModel) *twoBulkLOD {
	t := &twoBulkLOD{}
	for i := range hmmFGrid {
		f := (float64(i) + 0.5) / hmmFGrid
		t.fGrid = append(t.fGrid, f)
		t.tabH = append(t.tabH, binomPMF(mh.U, f))
		t.tabL = append(t.tabL, binomPMF(ml.U, f))
	}
	return t
}

// lod: both L vectors have max 1, so the numerator is 0 and the LOD is minus
// the log10 of the best shared-frequency likelihood.
func (t *twoBulkLOD) lod(LH, LL []float32) float64 {
	ll := make([]float64, len(t.fGrid))
	best := 0
	for i := range t.fGrid {
		sh, sl := 0.0, 0.0
		for j, v := range LH {
			sh += t.tabH[i][j] * float64(v)
		}
		for j, v := range LL {
			sl += t.tabL[i][j] * float64(v)
		}
		ll[i] = math.Log(math.Max(sh, hmmTiny)) + math.Log(math.Max(sl, hmmTiny))
		if ll[i] > ll[best] {
			best = i
		}
	}
	mx := ll[best]
	if best > 0 && best < len(ll)-1 { // parabolic refinement of the maximum
		y0, y1, y2 := ll[best-1], ll[best], ll[best+1]
		if den := y0 - 2*y1 + y2; den < 0 {
			mx = y1 - (y0-y2)*(y0-y2)/(8*den)
		}
	}
	return math.Max(0, -mx/math.Ln10)
}

// lodUpperStride picks the coarse subset of the f grid (every 8th point from
// the 4th) that lodUpper evaluates.
const lodUpperStride = 8

// lodUpper is an upper bound on lod(LH, LL) at about 1/8 of its cost: the
// shared-frequency log likelihood maximised over a subset of the f grid is at
// most its maximum over the whole grid, which is at most the parabolic
// refinement, so the bound is ≥ the exact LOD. The null simulation only needs
// each genome's maximum LOD, so it skips the exact LOD wherever the bound is
// already below the running maximum.
func (t *twoBulkLOD) lodUpper(LH, LL []float32) float64 {
	best := math.Inf(-1)
	for i := lodUpperStride / 2; i < len(t.fGrid); i += lodUpperStride {
		sh, sl := 0.0, 0.0
		for j, v := range LH {
			sh += t.tabH[i][j] * float64(v)
		}
		for j, v := range LL {
			sl += t.tabL[i][j] * float64(v)
		}
		best = math.Max(best, math.Log(math.Max(sh, hmmTiny))+math.Log(math.Max(sl, hmmTiny)))
	}
	return math.Max(0, -best/math.Ln10)
}

// hmmSetup is everything the HMM needs for one run of a mode.
type hmmSetup struct {
	high, low *poolModel // nil if the mode lacks that bulk
	two       *twoBulkLOD
	p0        float64 // theoretical p0 (one-bulk shift reference)
	cmPerMb   float64 // 0 = 1 Morgan per chromosome
}

func newHMMSetup(cfg utils.AnalysisConfig, bsaType string) (*hmmSetup, error) {
	hb, lb, both, _ := BulkFlags(bsaType)
	high, low, one, err := cfg.BulkSizes(bsaType)
	if err != nil {
		return nil, err
	}
	p0, err := ExpectedAF(cfg.Population)
	if err != nil {
		return nil, err
	}
	s := &hmmSetup{p0: p0, cmPerMb: cfg.CMPerMb}
	mk := func(n int) (*poolModel, error) {
		m, err := newPoolModel(cfg.Population, n)
		return &m, err
	}
	switch {
	case both:
		if s.high, err = mk(high); err != nil {
			return nil, err
		}
		if s.low, err = mk(low); err != nil {
			return nil, err
		}
		s.two = newTwoBulkLOD(s.high, s.low)
	case hb:
		s.high, err = mk(one)
	case lb:
		s.low, err = mk(one)
	}
	return s, err
}

// HMMAvailable reports whether the HMM LOD is computed for a mode: it needs
// the high-parent allele at every SNP, which 2b without --ref-parent lacks.
func HMMAvailable(cfg utils.AnalysisConfig, bsaType string) bool {
	return cfg.SignedDeltaSI(bsaType)
}

// computeHMMLOD fills raw[i].HMMLOD for every site, one chromosome per
// worker, using the shared marker order (chromosomes sorted by position).
func computeHMMLOD(s *hmmSetup, raw []BSAstats, order *markerOrder) {
	computeHMMLODWorkers(s, raw, order, s.workers(order))
}

// hmmMemBudget bounds the HMM's per-chromosome state held at once (the
// forward probabilities and emissions of every bin, float32).
const hmmMemBudget = 2 << 30

// workers is how many chromosomes the HMM can run at once: every core, unless
// the per-chromosome state of the largest chromosome would exceed
// hmmMemBudget.
func (s *hmmSetup) workers(order *markerOrder) int {
	maxN := 0
	for _, idx := range order.chroms {
		maxN = max(maxN, len(idx))
	}
	return max(1, min(runtime.GOMAXPROCS(0), int(hmmMemBudget/max(int64(1), s.chromBytes(maxN)))))
}

// chromBytes bounds the HMM state for a chromosome of n markers (at most n
// bins): alpha and emissions per bulk, float32.
func (s *hmmSetup) chromBytes(n int) int64 {
	u := 0
	for _, m := range []*poolModel{s.high, s.low} {
		if m != nil {
			u += m.U + 1
		}
	}
	return int64(n) * int64(u) * 2 * 4
}

// computeHMMLODWorkers is computeHMMLOD with a fixed number of workers (1 runs
// the chromosomes in order on the calling goroutine).
func computeHMMLODWorkers(s *hmmSetup, raw []BSAstats, order *markerOrder, workers int) {
	if workers <= 1 {
		for _, idx := range order.chroms {
			s.chromLOD(raw, idx, order.rep)
		}
		return
	}
	jobs := make(chan []int)
	var wg sync.WaitGroup
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for idx := range jobs {
				s.chromLOD(raw, idx, order.rep)
			}
		}()
	}
	for _, idx := range order.chroms {
		jobs <- idx
	}
	close(jobs)
	wg.Wait()
}

// maxHMMLOD returns the largest HMM LOD over raw's sites without storing
// per-site values (the null simulation needs only each genome's maximum).
// Chromosomes run in order on the calling goroutine. It equals the maximum
// of computeHMMLOD's values exactly: two-bulk bins whose lodUpper bound is
// below the running maximum cannot hold the maximum and are skipped.
func (s *hmmSetup) maxHMMLOD(raw []BSAstats, order *markerOrder) float64 {
	best := 0.0
	for _, idx := range order.chroms {
		best = s.chromMaxLOD(raw, idx, order.rep, best)
	}
	return best
}

// maxHMMLODWorkers is maxHMMLOD with the chromosomes shared among workers
// goroutines (the maximum, and so the result, does not depend on the split).
func (s *hmmSetup) maxHMMLODWorkers(raw []BSAstats, order *markerOrder, workers int) float64 {
	if workers <= 1 {
		return s.maxHMMLOD(raw, order)
	}
	jobs := make(chan []int)
	var mu sync.Mutex
	var wg sync.WaitGroup
	best := 0.0
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			local := 0.0
			for idx := range jobs {
				local = s.chromMaxLOD(raw, idx, order.rep, local)
			}
			mu.Lock()
			best = math.Max(best, local)
			mu.Unlock()
		}()
	}
	for _, idx := range order.chroms {
		jobs <- idx
	}
	close(jobs)
	wg.Wait()
	return best
}

// chromLOD computes the HMM LOD for one chromosome's sites (indices into raw,
// sorted by position). Only representative views (rep) contribute reads: the
// views of a split multi-allelic site share their reads, so multiplying them
// would count those reads more than once.
func (s *hmmSetup) chromLOD(raw []BSAstats, idx []int, rep []bool) {
	binStart, lod := s.chromBinLODs(raw, idx, rep)
	for b := 0; b+1 < len(binStart); b++ {
		v := lod(b)
		for x := binStart[b]; x < binStart[b+1]; x++ {
			raw[idx[x]].HMMLOD = v
		}
	}
}

// chromMaxLOD returns max(floor, the chromosome's largest bin LOD).
func (s *hmmSetup) chromMaxLOD(raw []BSAstats, idx []int, rep []bool, floor float64) float64 {
	binStart, lod, upper := s.chromBins(raw, idx, rep)
	best := floor
	for b := 0; b+1 < len(binStart); b++ {
		if upper != nil && upper(b) <= best {
			continue
		}
		best = math.Max(best, lod(b))
	}
	return best
}

// chromBinLODs runs the HMM on one chromosome and returns its bins (binStart,
// with a final sentinel) and the LOD of each bin.
func (s *hmmSetup) chromBinLODs(raw []BSAstats, idx []int, rep []bool) ([]int, func(int) float64) {
	binStart, lod, _ := s.chromBins(raw, idx, rep)
	return binStart, lod
}

// chromBins runs forward–backward on one chromosome and returns its bins, the
// LOD of a bin, and (two-bulk) a cheap upper bound on it (nil otherwise).
func (s *hmmSetup) chromBins(raw []BSAstats, idx []int, rep []bool) ([]int, func(int) float64, func(int) float64) {
	n := len(idx)
	gaps := gapsMorgans(raw, idx, s.cmPerMb)

	// Bins of markers within hmmBinMorgans of the bin's first marker.
	binStart := []int{0}
	binDM := []float64{0}
	acc := 0.0
	for x := 1; x < n; x++ {
		acc += gaps[x]
		if acc >= hmmBinMorgans {
			binStart = append(binStart, x)
			binDM = append(binDM, acc)
			acc = 0
		}
	}
	binStart = append(binStart, n)

	bulk := func(high bool) hmmBulkData {
		d := hmmBulkData{h: make([]int, n), d: make([]int, n)}
		oneBulk := s.two == nil
		if oneBulk {
			d.shift = make([]float64, n)
		}
		for x, i := range idx {
			r := &raw[i]
			if !rep[i] { // non-representative view: no reads (d = 0)
				continue
			}
			if high {
				d.h[x], d.d[x] = r.HighBulkH, r.HighBulkH+r.HighBulkL
			} else {
				d.h[x], d.d[x] = r.LowBulkH, r.LowBulkH+r.LowBulkL
			}
			if oneBulk {
				d.shift[x] = r.OneBulkP0 - s.p0
			}
		}
		return d
	}

	switch {
	case s.two != nil:
		LH := s.high.chromLikelihoods(bulk(true), binStart, binDM)
		LL := s.low.chromLikelihoods(bulk(false), binStart, binDM)
		return binStart, func(b int) float64 { return s.two.lod(LH[b], LL[b]) },
			func(b int) float64 { return s.two.lodUpper(LH[b], LL[b]) }
	case s.high != nil:
		L := s.high.chromLikelihoods(bulk(true), binStart, binDM)
		return binStart, func(b int) float64 { return s.high.lodOneBulk(L[b]) }, nil
	case s.low != nil:
		L := s.low.chromLikelihoods(bulk(false), binStart, binDM)
		return binStart, func(b int) float64 { return s.low.lodOneBulk(L[b]) }, nil
	}
	return binStart, func(int) float64 { return 0 }, nil
}
