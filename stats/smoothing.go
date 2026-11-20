package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// minEffectiveMarkers is the smallest effective number of markers,
// n_eff = (Σw)²/Σw² over the kernel weights, a marker needs to be smoothed.
// Unlike a plain count within ±3σ, it ignores markers whose Gaussian weight is
// tiny, so an isolated marker (whose smoothed value would be essentially its
// own) is dropped instead of producing a spike.
const minEffectiveMarkers = 5.0

const kernelCutoffSigmas = 3.0

type SmoothedStats struct {
	CHROM string
	POS   int64
	REF   string
	ALT   string
	Depth int // effective depth (see BSAstats.Depth)
	// Per-bulk depths of the site, for two-bulk thresholds and BRM.
	HighDepth int
	LowDepth  int
	// NEff is the effective number of markers behind the smoothed values,
	// (Σw)²/Σw² over the kernel weights.
	NEff float64

	// ── smoothed raw values ─────────────────────────────────────────────
	// Two-bulk
	SmHighSI  float64
	SmLowSI   float64
	SmDeltaSI float64
	SmGstat   float64
	SmED      float64
	// One-bulk
	OneBulkP0  float64 // the site's own one-bulk p0 (not smoothed; differs by polarity)
	SmAFDev    float64
	SmOneBulkG float64
	// HMMLOD is the pool-HMM LOD at the site (hmm.go): two-bulk, bulks differ;
	// one-bulk, the bulk deviates from p0. Already linkage-aware, so it is
	// carried through unsmoothed. 0 where the HMM is not computed (2b without
	// --ref-parent).
	HMMLOD float64

	// ── robust Z-scores of smoothed values ──────────────────────────────
	ZHighSI   float64
	ZLowSI    float64
	ZDeltaSI  float64
	ZGstat    float64
	ZED       float64
	ZAFDev    float64
	ZOneBulkG float64
	ZHMMLOD   float64

	// ── consolidated scores ──────────────────────────────────────────────
	// CompositeZ = (|ZDeltaSI| + ZGstat + ZHMMLOD) / √k in two-bulk modes (see
	// consolidate). It is not bounded below: ZGstat and ZHMMLOD are negative
	// where the value is under the genome-wide median.
	CompositeZ float64
	// max|Z| across all valid statistics (unsigned)
	MaxAbsZ float64
}

func SmoothAndNormalise(cfg utils.AnalysisConfig, bsaType string, rawStats []BSAstats) ([]SmoothedStats, error) {
	if len(rawStats) == 0 {
		return nil, fmt.Errorf("no raw stats supplied to SmoothAndNormalise")
	}
	if cfg.WindowSize <= 0 {
		return nil, fmt.Errorf("WindowSize must be > 0, got %d", cfg.WindowSize)
	}

	color.Cyan("\n============================ Smoothing & Normalising (%s) ============================\n\n", bsaType)
	color.Cyan("Gaussian kernel σ = %.0f bp  (WindowSize=%d / 2)\n\n", float64(cfg.WindowSize)/2, cfg.WindowSize)
	if HMMAvailable(cfg, bsaType) {
		if cfg.CMPerMb > 0 {
			color.Cyan("HMM LOD: recombination %.3g cM/Mb\n", cfg.CMPerMb)
		} else {
			color.Cyan("HMM LOD: recombination assumed 1 Morgan per chromosome (set --cm-per-mb to change)\n")
		}
	} else {
		color.Yellow("HMM LOD not computed: 2b without --ref-parent has no allele polarity (CompositeZ uses |ΔSI| and G only)\n")
	}

	sc, err := newScorer(cfg, bsaType)
	if err != nil {
		return nil, err
	}
	smoothed, err := sc.score(rawStats, true, 1)
	if err != nil {
		return nil, err
	}

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.smoothed_and_normalised.tsv", bsaType))
	if err := writeSmoothedTSV(outPath, smoothed, sc.lay); err != nil {
		return nil, fmt.Errorf("writing %s: %w", outPath, err)
	}

	color.Green("Smoothed & normalised stats written to %s (%d variants)\n", outPath, len(smoothed))
	return smoothed, nil
}

// scorer runs the in-memory scoring pipeline (HMM LOD → kernel smoothing →
// robust Z → CompositeZ). The real data and every null-simulation replicate
// go through it, so thresholds are calibrated on exactly what is scanned.
type scorer struct {
	bw     float64
	step   float64 // smoothing grid spacing (gridStep)
	weight func(*BSAstats) float64
	hmm    *hmmSetup // nil when the HMM LOD is not computed
	lay    smoothLayout
}

func newScorer(cfg utils.AnalysisConfig, bsaType string) (*scorer, error) {
	hb, lb, _, _ := BulkFlags(bsaType)
	weight, err := markerWeight(cfg, bsaType)
	if err != nil {
		return nil, err
	}
	bw := float64(cfg.WindowSize) / 2
	sc := &scorer{bw: bw, step: gridStep(int64(cfg.StepSize), bw), weight: weight,
		lay: smoothLayout{highBulk: hb, lowBulk: lb, signedDSI: cfg.SignedDeltaSI(bsaType), hmm: HMMAvailable(cfg, bsaType)}}
	if sc.lay.hmm {
		if sc.hmm, err = newHMMSetup(cfg, bsaType); err != nil {
			return nil, err
		}
	}
	return sc, nil
}

// score fills raw[i].HMMLOD and returns the smoothed, normalised and
// consolidated markers. report prints the sparse-marker summary. every > 1
// smooths only every n-th marker of each chromosome (all markers still feed
// the kernel); the null simulation uses it, since it only samples CompositeZ.
func (sc *scorer) score(raw []BSAstats, report bool, every int) ([]SmoothedStats, error) {
	return sc.scoreOrdered(raw, newMarkerOrder(raw), report, every)
}

// scoreOrdered is score with a precomputed marker order (the null simulation
// builds it once and reuses it for every replicate, whose markers are the
// same). With every > 1 the output is compact and in chromosome order, not
// input order.
func (sc *scorer) scoreOrdered(raw []BSAstats, order *markerOrder, report bool, every int) ([]SmoothedStats, error) {
	l := sc.lay
	both, one := l.highBulk && l.lowBulk, l.highBulk != l.lowBulk
	if sc.hmm != nil {
		computeHMMLOD(sc.hmm, raw, order)
	}
	smoothed, err := gaussianSmoothOpts(raw, smoothOpts{bw: sc.bw, step: sc.step, weight: sc.weight,
		hasHighBulk: l.highBulk, hasLowBulk: l.lowBulk, signedDSI: l.signedDSI, report: report, every: every, order: order})
	if err != nil {
		return nil, err
	}
	normalise(smoothed, l.highBulk, l.lowBulk, both, one, l.hmm)
	consolidate(smoothed, l.highBulk, l.lowBulk, both, one, l.signedDSI, l.hmm)
	return smoothed, nil
}

// ─── smoothing ─────────────────────────────────────────────────────────────────

// markerOrder groups marker indices by chromosome (chromosomes sorted by
// name), each group sorted by position (stable, so duplicate positions keep
// input order). It depends only on CHROM/POS, so it is built once per marker
// set and shared by the HMM and the smoother.
type markerOrder struct {
	chroms [][]int
}

func newMarkerOrder(raw []BSAstats) *markerOrder {
	byChrom := map[string][]int{}
	for i := range raw {
		byChrom[raw[i].CHROM] = append(byChrom[raw[i].CHROM], i)
	}
	names := make([]string, 0, len(byChrom))
	for name := range byChrom {
		names = append(names, name)
	}
	sort.Strings(names)
	o := &markerOrder{chroms: make([][]int, 0, len(names))}
	for _, name := range names {
		idx := byChrom[name]
		sort.SliceStable(idx, func(a, b int) bool { return raw[idx[a]].POS < raw[idx[b]].POS })
		o.chroms = append(o.chroms, idx)
	}
	return o
}

// markerWeight returns each marker's inverse-variance weight for smoothing.
//
// A bulk's allele frequency estimate has two variance components (per unit
// p(1−p)): pool sampling of the plants, 1/(c·n), and read sampling, 1/d, with
// n the bulk size, c the population variance scale and d the depth. ΔSI adds
// both bulks. The weight is 1 / the total, so once depth is well above the
// pool size, pool sampling dominates and deeper markers gain little extra
// weight. (Plain depth weights gave a 200x marker 5× the weight of a 40x one
// with 20-plant F2 bulks, though it is only ~1.7× as precise, and deep
// markers are where collapsed repeats concentrate.) The same weight is used
// for every statistic: it is the precision of the underlying allele
// frequencies.
func markerWeight(cfg utils.AnalysisConfig, bsaType string) (func(*BSAstats) float64, error) {
	high, low, one, err := cfg.BulkSizes(bsaType)
	if err != nil {
		return nil, err
	}
	c := PopulationVarianceScale(cfg.Population)
	if _, _, both, _ := BulkFlags(bsaType); both {
		pool := 1/(c*float64(high)) + 1/(c*float64(low))
		return func(s *BSAstats) float64 {
			dH, dL := s.HighDepth, s.LowDepth
			if dH <= 0 || dL <= 0 {
				dH, dL = s.Depth, s.Depth
			}
			return 1 / (pool + 1/math.Max(float64(dH), 1) + 1/math.Max(float64(dL), 1))
		}, nil
	}
	pool := 1 / (c * float64(one))
	return func(s *BSAstats) float64 { return 1 / (pool + 1/math.Max(float64(s.Depth), 1)) }, nil
}

// smoothKernel is how a statistic is smoothed.
type smoothKernel uint8

const (
	// kGauss: Gaussian × marker weights (SI, ΔSI, AFDev): the weighted local
	// average in the interior, and a local linear fit where a chromosome end
	// cuts the kernel off (see edgeBlend), which removes the average's bias
	// toward the interior there.
	kGauss smoothKernel = iota
	// kTricube: tricube kernel average over ±h, Magwene's G′ (G, G1).
	kTricube
	// kLoess: local quadratic regression with tricube weights over ±h, as in
	// MMAPPR's loess of ED⁴.
	kLoess
)

// smoothField pairs a raw statistic with the smoothed field it fills; lo/hi
// clamp the fitted value to the statistic's range (local regressions can
// overshoot at chromosome ends).
type smoothField struct {
	get    func(*BSAstats) float64
	set    func(*SmoothedStats, float64)
	kernel smoothKernel
	lo, hi float64
}

// smoothFields lists the statistics smoothed for a mode's layout, grouped by
// kernel: [Gaussian | tricube | loess].
func smoothFields(hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk, signedDSI bool) []smoothField {
	inf := math.Inf(1)
	var g, t, l []smoothField
	if hasHighBulk {
		g = append(g, smoothField{func(s *BSAstats) float64 { return s.HighSI }, func(s *SmoothedStats, v float64) { s.SmHighSI = v }, kGauss, 0, 1})
	}
	if hasLowBulk {
		g = append(g, smoothField{func(s *BSAstats) float64 { return s.LowSI }, func(s *SmoothedStats, v float64) { s.SmLowSI = v }, kGauss, 0, 1})
	}
	if hasBothBulks {
		dsiLo := -1.0
		if !signedDSI {
			dsiLo = 0 // |ΔSI|
		}
		g = append(g, smoothField{func(s *BSAstats) float64 { return s.DeltaSI }, func(s *SmoothedStats, v float64) { s.SmDeltaSI = v }, kGauss, dsiLo, 1})
		t = append(t, smoothField{func(s *BSAstats) float64 { return s.Gstat }, func(s *SmoothedStats, v float64) { s.SmGstat = v }, kTricube, 0, inf})
		l = append(l, smoothField{func(s *BSAstats) float64 { return s.ED }, func(s *SmoothedStats, v float64) { s.SmED = v }, kLoess, 0, inf})
	}
	if hasOneBulk {
		g = append(g, smoothField{func(s *BSAstats) float64 { return s.OneBulkAFDev }, func(s *SmoothedStats, v float64) { s.SmAFDev = v }, kGauss, -1, 1})
		t = append(t, smoothField{func(s *BSAstats) float64 { return s.OneBulkGstat }, func(s *SmoothedStats, v float64) { s.SmOneBulkG = v }, kTricube, 0, inf})
	}
	return append(append(g, t...), l...)
}

// smoothOpts configures gaussianSmooth.
type smoothOpts struct {
	bw                                 float64                 // Gaussian σ in bp (WindowSize/2); tricube half-width h = 2·bw
	step                               float64                 // grid spacing for dense chromosomes; 0 = always exact at markers
	weight                             func(*BSAstats) float64 // marker weight (markerWeight)
	hasHighBulk, hasLowBulk, signedDSI bool
	report                             bool         // print the sparse-marker summary
	every                              int          // > 1: output only every n-th marker of each chromosome (compact, chromosome order)
	order                              *markerOrder // precomputed marker order; nil = build it
}

// gridStep is the smoothing grid spacing: --step-size, capped at σ/20 so
// linear interpolation between grid points stays accurate (the smoothed
// tracks vary on the scale of σ; interpolation error falls with the square
// of the spacing, ≈ 0.25% of a track's range at σ/20 on pure-noise input).
func gridStep(stepSize int64, bw float64) float64 {
	s := bw / 20
	if stepSize > 0 && float64(stepSize) < s {
		s = float64(stepSize)
	}
	return s
}

// smoothChunk bounds the number of query points a worker smooths at a time.
// Chromosomes are processed one at a time, so chunks are sized to give every
// core several jobs even on short chromosomes or coarse grids.
const smoothChunk = 256

// gaussianSmooth smooths every statistic with default options (exact at every
// marker, signed ΔSI); see gaussianSmoothOpts.
func gaussianSmooth(raw []BSAstats, bw float64, weight func(*BSAstats) float64, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk, report bool) ([]SmoothedStats, error) {
	return gaussianSmoothOpts(raw, smoothOpts{bw: bw, weight: weight, hasHighBulk: hasHighBulk, hasLowBulk: hasLowBulk,
		signedDSI: true, report: report, every: 1})
}

// gaussianSmoothOpts smooths every statistic per chromosome, with the kernel
// set by its smoothField. At a query position q, with Δ = pos_m − q:
//
//	kGauss:   weights w_m = weight(m)·exp(−½(Δ/bw)²), |Δ| ≤ 3·bw; the
//	          weighted mean, blended near chromosome ends into β0 of the
//	          weighted linear fit x ≈ β0 + β1·Δ (edgeBlend)
//	kTricube: Σ t_m·x_m / Σ t_m,  t_m = weight(m)·(1 − |Δ/h|³)³, |Δ| < h
//	          (Magwene's G′)
//	kLoess:   β0 of the weighted fit x ≈ β0 + β1·u + β2·u², u = Δ/h,
//	          weights t_m (local quadratic, as R's loess)
//
// h = 2·bw = WindowSize (QTLseqr's windowSize). Fitted values are clamped to
// each statistic's range. Markers whose effective number of markers under the
// Gaussian weights, (Σw)²/Σw², is below minEffectiveMarkers are dropped.
//
// Query points. On a chromosome with more markers than grid points (spacing
// gridStep), the kernels are evaluated exactly on a regular grid from the
// first to the last marker and linearly interpolated to each marker (values
// and NEff), so the cost grows with chromosome length, not marker density.
// Otherwise they are evaluated exactly at every marker.
//
// Memory. Chromosomes are processed one at a time: a chromosome's positions,
// weights, values and query results are allocated, used and released before
// the next. Within a chromosome, query points are smoothed in parallel chunks;
// all kernels share one pass over the neighbours, and each sum runs over the
// neighbours in position order, so results do not depend on scheduling.
//
// The output keeps the input (raw) order. With every > 1 (null simulation) it
// is compact (only the output markers, no len(raw) array) and in chromosome
// order.
func gaussianSmoothOpts(raw []BSAstats, o smoothOpts) ([]SmoothedStats, error) {
	bw := o.bw
	cutoff := kernelCutoffSigmas * bw
	h := 2 * bw // tricube half-width = WindowSize (inside the ±3σ window)
	invBw := 1 / bw
	hasBothBulks, hasOneBulk := o.hasHighBulk && o.hasLowBulk, o.hasHighBulk != o.hasLowBulk
	fields := smoothFields(o.hasHighBulk, o.hasLowBulk, hasBothBulks, hasOneBulk, o.signedDSI)
	nf := len(fields)
	ng, nt := 0, 0
	for _, fd := range fields {
		switch fd.kernel {
		case kGauss:
			ng++
		case kTricube:
			nt++
		}
	}
	nl := nf - ng - nt
	useTri := nt+nl > 0

	order := o.order
	if order == nil {
		order = newMarkerOrder(raw)
	}
	every := max(o.every, 1)
	var out []SmoothedStats
	var keep []bool
	if every > 1 {
		nOut := 0
		for _, idx := range order.chroms {
			nOut += (len(idx) + every - 1) / every
		}
		out = make([]SmoothedStats, 0, nOut)
	} else {
		out = make([]SmoothedStats, len(raw))
		keep = make([]bool, len(raw))
	}
	var drops []chromDrop

	workers := max(runtime.GOMAXPROCS(0), 1)
	// One persistent worker pool; chromosomes are fed one at a time.
	type job struct {
		q0, q1 int
		c      *chromWork
		done   *sync.WaitGroup
	}
	jobs := make(chan job, 64)
	var pool sync.WaitGroup
	for range workers {
		pool.Add(1)
		go func() {
			defer pool.Done()
			accG := make([]float64, 2*ng) // Σw·x, Σw·d·x per Gaussian field
			accT := make([]float64, nt)   // Σt·x per tricube field
			accL := make([]float64, 3*nl) // Σt·x, Σt·u·x, Σt·u²·x per loess field
			for jb := range jobs {
				c := jb.c
				pos, n := c.npos, len(c.npos) // kernel neighbours (markers, or aggregated grid intervals)
				first, last := c.pos[0], c.pos[len(c.pos)-1]
				edgeFit := last-first >= 2*cutoff // local linear near ends only on chromosomes ≥ 6σ
				q0 := c.q[jb.q0]
				lo := sort.Search(n, func(m int) bool { return pos[m] >= q0-cutoff })
				hi := sort.Search(n, func(m int) bool { return pos[m] > q0+cutoff })
				tlo := sort.Search(n, func(m int) bool { return pos[m] > q0-h })
				thi := sort.Search(n, func(m int) bool { return pos[m] >= q0+h })
				for qi := jb.q0; qi < jb.q1; qi++ {
					pK := c.q[qi]
					for lo < n && pos[lo] < pK-cutoff {
						lo++
					}
					for hi < n && pos[hi] <= pK+cutoff {
						hi++
					}
					for tlo < n && pos[tlo] <= pK-h { // tricube window: |Δ| < h
						tlo++
					}
					for thi < n && pos[thi] < pK+h {
						thi++
					}

					clear(accG)
					clear(accT)
					clear(accL)
					var g0, g1, g2, w2 float64     // Gaussian moments Σw·dⁱ, Σw²
					var s0, s1, s2, s3, s4 float64 // tricube moments Σt·uⁱ
					for m := lo; m < hi; m++ {
						d := (pos[m] - pK) * invBw
						k := math.Exp(-0.5 * d * d)
						w := c.nwt[m] * k // marker weight × Gaussian
						wd := w * d
						g0 += w
						g1 += wd
						g2 += wd * d
						if c.nwt2 == nil {
							w2 += w * w
						} else { // aggregated: Σ(wt·K)² = K²·Σwt²
							w2 += k * k * c.nwt2[m]
						}
						row := c.nvals[m*nf : (m+1)*nf : (m+1)*nf]
						for f := range ng {
							accG[2*f] += w * row[f]
							accG[2*f+1] += wd * row[f]
						}
						if m < tlo || m >= thi || !useTri {
							continue
						}
						u := 0.5 * d // Δ/h, since h = 2·bw
						a := math.Abs(u)
						q := 1 - a*a*a
						tw := c.nwt[m] * q * q * q // marker weight × tricube
						s0 += tw
						for f := range nt {
							accT[f] += tw * row[ng+f]
						}
						if nl > 0 {
							u2 := u * u
							s1 += tw * u
							s2 += tw * u2
							s3 += tw * u2 * u
							s4 += tw * u2 * u2
							for f := range nl {
								y := tw * row[ng+nt+f]
								accL[3*f] += y
								accL[3*f+1] += y * u
								accL[3*f+2] += y * u2
							}
						}
					}
					res := c.res[qi*nf : (qi+1)*nf]
					c.neff[qi] = 0
					if w2 > 0 {
						c.neff[qi] = g0 * g0 / w2
					}
					nG, nT := hi-lo, thi-tlo
					alpha := 0.0
					if edgeFit {
						alpha = edgeBlend(math.Min(pK-first, last-pK), bw)
					}
					for f := range nf {
						fd := &fields[f]
						var v float64
						switch {
						case f < ng:
							mean := 0.0
							if g0 > 0 {
								mean = accG[2*f] / g0
							}
							v = mean
							if alpha > 0 {
								v += alpha * (localLinear(g0, g1, g2, accG[2*f], accG[2*f+1], nG) - mean)
							}
						case f < ng+nt:
							if s0 > 0 {
								v = accT[f-ng] / s0
							}
						default:
							l := f - ng - nt
							v = loessFit(s0, s1, s2, s3, s4, accL[3*l], accL[3*l+1], accL[3*l+2], nT)
						}
						res[f] = math.Min(math.Max(v, fd.lo), fd.hi)
					}
				}
				jb.done.Done()
			}
		}()
	}

	for _, idx := range order.chroms {
		c := newChromWork(raw, idx, fields, o)
		var done sync.WaitGroup
		chunk := min(max(len(c.q)/(4*workers), 16), smoothChunk)
		for q0 := 0; q0 < len(c.q); q0 += chunk {
			done.Add(1)
			jobs <- job{q0, min(q0+chunk, len(c.q)), c, &done}
		}
		done.Wait()
		drops = append(drops, c.emit(raw, &out, keep, fields, o))
		// c is released here
	}
	close(jobs)
	pool.Wait()

	if keep != nil { // keep smoothed markers in input order
		valid := out[:0]
		for i := range out {
			if keep[i] {
				valid = append(valid, out[i])
			}
		}
		out = valid
	}
	if o.report {
		reportSparseDrops(drops, len(raw), len(raw)-len(out))
	}
	return out, nil
}

// chromWork is one chromosome's smoothing state: its markers (sorted by
// position), the query points and their results. It lives only while that
// chromosome is processed.
type chromWork struct {
	name string
	idx  []int     // indices into raw, by position
	pos  []float64 // marker positions
	// Kernel neighbours: the markers themselves, or in grid mode one
	// aggregated point per dense grid interval (at the weighted centroid,
	// carrying Σwt, Σwt² and the weighted mean values) plus the individual
	// markers of sparse intervals.
	npos  []float64 // neighbour positions
	nwt   []float64 // neighbour weights (Σ marker weights)
	nwt2  []float64 // Σ marker weight² per neighbour; nil when neighbours are markers
	nvals []float64 // neighbour values (weighted means), row-major [neighbour][field]
	q     []float64 // query positions, sorted (markers, or grid points plus sparse markers)
	grid  bool      // grid mode: dense markers are interpolated between grid points
	gq    []int     // grid mode: q-index of each grid point
	mq    []int     // grid mode: q-index of a marker evaluated exactly, or −1
	res   []float64 // results, row-major [query][field]
	neff  []float64 // NEff per query
}

// sparseIntervalMarkers: in grid mode, markers in a grid interval holding at
// most this many markers are evaluated exactly. In sparse stretches the
// smoothed curve changes quickly as single markers enter or leave the kernel,
// so interpolation would be inaccurate there; few markers means little cost.
const sparseIntervalMarkers = 4

func newChromWork(raw []BSAstats, idx []int, fields []smoothField, o smoothOpts) *chromWork {
	n, nf := len(idx), len(fields)
	c := &chromWork{name: raw[idx[0]].CHROM, idx: idx, pos: make([]float64, n)}
	wt, vals := make([]float64, n), make([]float64, n*nf)
	for k, ri := range idx {
		s := &raw[ri]
		c.pos[k] = float64(s.POS)
		wt[k] = o.weight(s)
		for f, fd := range fields {
			vals[k*nf+f] = fd.get(s)
		}
	}
	c.q, c.npos, c.nwt, c.nvals = c.pos, c.pos, wt, vals

	// Markers sharing a position (views of one split multi-allelic site) come
	// from the same reads and plants, so each position is one kernel
	// neighbour ("unit"): weight = the views' mean weight, value = their
	// weighted mean, Σw² = weight² (so NEff counts it once).
	u := newPosUnits(c.pos, wt, vals, nf)
	if u.n() < n {
		c.npos, c.nwt, c.nwt2, c.nvals = u.pos, u.w, u.w2, u.means()
	}
	if o.step > 0 {
		first, last := c.pos[0], c.pos[n-1]
		// A grid needs a positive span: markers all at one position (e.g. a
		// multi-allelic site split on a tiny contig) stay on the exact path.
		if nGrid := int((last-first)/o.step) + 2; last > first && nGrid < n {
			c.buildGrid(first, last, o.step)
			c.aggregate(u, nf, o.bw/aggregateBinSigmas)
		}
	}
	c.res = make([]float64, len(c.q)*nf)
	c.neff = make([]float64, len(c.q))
	return c
}

// posUnits are a chromosome's distinct marker positions: per position, the
// mean of its markers' weights (w), w² (w2), the weighted value sums
// (sums, row-major [unit][field]) and the index of its first marker.
type posUnits struct {
	pos, w, w2, sums []float64
	first            []int
	nf               int
}

func newPosUnits(pos, wt, vals []float64, nf int) *posUnits {
	u := &posUnits{nf: nf}
	for k := 0; k < len(pos); {
		j := k
		for j < len(pos) && pos[j] == pos[k] {
			j++
		}
		dup := float64(j - k)
		w := 0.0
		for m := k; m < j; m++ {
			w += wt[m] / dup
		}
		u.pos = append(u.pos, pos[k])
		u.w = append(u.w, w)
		u.w2 = append(u.w2, w*w)
		u.first = append(u.first, k)
		for f := range nf {
			s := 0.0
			for m := k; m < j; m++ {
				s += wt[m] / dup * vals[m*nf+f]
			}
			u.sums = append(u.sums, s)
		}
		k = j
	}
	return u
}

func (u *posUnits) n() int { return len(u.pos) }

// means returns the units' weighted-mean values.
func (u *posUnits) means() []float64 {
	m := make([]float64, len(u.sums))
	for i := range u.pos {
		for f := range u.nf {
			if u.w[i] > 0 {
				m[i*u.nf+f] = u.sums[i*u.nf+f] / u.w[i]
			}
		}
	}
	return m
}

// aggregateBinSigmas: in grid mode, dense stretches are aggregated into bins
// of width σ/aggregateBinSigmas (10 kb at the default σ = 1 Mb). A ±3σ window
// then holds ~600 points instead of every marker (~40k at 150 bp spacing).
const aggregateBinSigmas = 100

// aggregate builds the grid-mode kernel neighbours from the position units:
// units of dense grid intervals (whose markers are interpolated, c.mq < 0)
// are grouped into bins of width ≤ width, each becoming one point at its
// weighted centroid with weight Σw, Σw² and weighted-mean values; units
// evaluated exactly (sparse intervals) stay individual. Placing the point at
// the centroid cancels the first-order kernel error, Σw·(pos − centroid) = 0;
// what remains is second order in the bin width (σ/100).
func (c *chromWork) aggregate(u *posUnits, nf int, width float64) {
	c.npos, c.nwt, c.nwt2, c.nvals = nil, nil, nil, nil
	sums := make([]float64, nf)
	add := func(p, w, w2 float64) {
		c.npos = append(c.npos, p)
		c.nwt = append(c.nwt, w)
		c.nwt2 = append(c.nwt2, w2)
		for f := range nf {
			v := 0.0
			if w > 0 {
				v = sums[f] / w
			}
			c.nvals = append(c.nvals, v)
		}
	}
	g := 0 // current grid interval: grid[g] ≤ pos < grid[g+1]
	for k := 0; k < u.n(); {
		if c.mq[u.first[k]] >= 0 { // exactly evaluated (sparse interval): keep the unit
			copy(sums, u.sums[k*nf:(k+1)*nf])
			add(u.pos[k], u.w[k], u.w2[k])
			k++
			continue
		}
		for g+2 < len(c.gq) && c.q[c.gq[g+1]] <= u.pos[k] {
			g++
		}
		hi := c.q[c.gq[g+1]]
		if g+2 == len(c.gq) {
			hi = math.Inf(1) // last interval includes the last position
		}
		hi = math.Min(hi, u.pos[k]+width) // bin: within the interval and ≤ width
		var w, w2, wp float64
		clear(sums)
		j := k
		for ; j < u.n() && c.mq[u.first[j]] < 0 && u.pos[j] < hi; j++ {
			w += u.w[j]
			w2 += u.w2[j]
			wp += u.w[j] * u.pos[j]
			for f := range nf {
				sums[f] += u.sums[j*nf+f]
			}
		}
		centroid := u.pos[k]
		if w > 0 {
			centroid = wp / w
		}
		add(centroid, w, w2)
		k = j
	}
}

// buildGrid sets the query points: grid points every step from first to last,
// merged in order with the markers that lie in sparse grid intervals.
func (c *chromWork) buildGrid(first, last, step float64) {
	var grid []float64
	for p := first; p < last; p += step {
		grid = append(grid, p)
	}
	grid = append(grid, last)
	// Markers per grid interval [grid[g], grid[g+1]).
	count := make([]int, len(grid))
	g := 0
	for _, p := range c.pos {
		for g+2 < len(grid) && grid[g+1] <= p {
			g++
		}
		count[g]++
	}
	c.grid = true
	c.mq = make([]int, len(c.pos))
	c.q = make([]float64, 0, len(grid))
	c.gq = make([]int, len(grid))
	g, k := 0, 0
	for g < len(grid) || k < len(c.pos) {
		// Emit grid points before markers at the same position, so markers in
		// an interval follow its left grid point.
		if k == len(c.pos) || (g < len(grid) && grid[g] <= c.pos[k]) {
			c.gq[g] = len(c.q)
			c.q = append(c.q, grid[g])
			g++
			continue
		}
		gi := max(g-1, 0) // interval of marker k
		if count[gi] <= sparseIntervalMarkers {
			c.mq[k] = len(c.q)
			c.q = append(c.q, c.pos[k])
		} else {
			c.mq[k] = -1
		}
		k++
	}
}

// chromDrop records how many markers a chromosome lost to the sparse rule.
type chromDrop struct {
	name        string
	lost, total int
}

// emit writes the chromosome's smoothed markers (interpolating from the grid
// when one was used): into (*out)[ri], marking keep[ri], or, when keep is nil
// (compact output), appended to *out.
func (c *chromWork) emit(raw []BSAstats, out *[]SmoothedStats, keep []bool, fields []smoothField, o smoothOpts) chromDrop {
	nf := len(fields)
	every := max(o.every, 1)
	hasOneBulk := o.hasHighBulk != o.hasLowBulk
	d := chromDrop{name: c.name, total: len(c.idx)}
	g := 0 // grid interval containing the marker: grid[g] ≤ pos < grid[g+1]
	for k, ri := range c.idx {
		if every > 1 && k%every != 0 {
			continue
		}
		qk := k // query index holding this marker's exact result
		if c.grid {
			qk = c.mq[k]
		}
		var nEff float64
		val := func(f int) float64 { return c.res[qk*nf+f] }
		if qk >= 0 {
			nEff = c.neff[qk]
		} else { // interpolate between the grid points around the marker
			p := c.pos[k]
			for g+2 < len(c.gq) && c.q[c.gq[g+1]] <= p {
				g++
			}
			a, b := c.gq[g], c.gq[g+1]
			t := 0.0
			if span := c.q[b] - c.q[a]; span > 0 {
				t = math.Min(math.Max((p-c.q[a])/span, 0), 1)
			}
			nEff = c.neff[a] + t*(c.neff[b]-c.neff[a])
			val = func(f int) float64 {
				x, y := c.res[a*nf+f], c.res[b*nf+f]
				return x + t*(y-x)
			}
		}
		if nEff < minEffectiveMarkers {
			d.lost++
			continue
		}
		r := &raw[ri]
		var s *SmoothedStats
		if keep != nil {
			s = &(*out)[ri]
			keep[ri] = true
		} else {
			*out = append(*out, SmoothedStats{})
			s = &(*out)[len(*out)-1]
		}
		s.CHROM, s.POS, s.REF, s.ALT = r.CHROM, r.POS, r.REF, r.ALT
		s.Depth, s.HighDepth, s.LowDepth = r.Depth, r.HighDepth, r.LowDepth
		s.NEff = nEff
		s.HMMLOD = r.HMMLOD
		if hasOneBulk {
			s.OneBulkP0 = r.OneBulkP0
		}
		for f := range fields {
			fields[f].set(s, val(f))
		}
	}
	if every > 1 { // markers not output are not "lost"
		d.total = (len(c.idx) + every - 1) / every
	}
	return d
}

// edgeBlend is the weight of the local linear fit for Gaussian tracks (on
// chromosomes spanning at least two full kernel widths, 6σ: on shorter
// contigs every marker is near an end, and a local linear fit there would
// only add variance, so they keep the local average) at a
// query dist bp from the nearer chromosome end: 1 within 2σ (the kernel is
// cut off and the local average is biased toward the interior), 0 beyond the
// 3σ cutoff (the kernel is whole; the local average has lower variance, which
// matters in sparse stretches where a local linear fit extrapolates a
// noise-driven slope), linear in between so the track has no step.
func edgeBlend(dist, bw float64) float64 {
	return math.Min(math.Max((kernelCutoffSigmas*bw-dist)/bw, 0), 1)
}

// localLinear returns β0 of the weighted linear fit y ≈ β0 + β1·d from the
// moments g_i = Σw·dⁱ and r_i = Σw·dⁱ·y; the weighted mean when the fit is
// ill-conditioned (one point, or all at one position).
func localLinear(g0, g1, g2, r0, r1 float64, n int) float64 {
	if g0 <= 0 {
		return 0
	}
	if n >= 2 {
		if den := g0*g2 - g1*g1; den > 1e-12*g0*g2 {
			return (g2*r0 - g1*r1) / den
		}
	}
	return r0 / g0
}

// loessFit returns β0 of the weighted least-squares fit y ≈ β0 + β1·u + β2·u²
// from the moments s_i = Σt·uⁱ and r_i = Σt·uⁱ·y, i.e. the fitted value at the
// query (u = 0). With fewer than 3 points or a near-singular system it falls
// back to the local linear fit, then to the weighted mean.
func loessFit(s0, s1, s2, s3, s4, r0, r1, r2 float64, n int) float64 {
	if s0 <= 0 {
		return 0
	}
	if n >= 3 {
		// Cramer's rule on [s0 s1 s2; s1 s2 s3; s2 s3 s4]·β = [r0 r1 r2].
		det := s0*(s2*s4-s3*s3) - s1*(s1*s4-s3*s2) + s2*(s1*s3-s2*s2)
		if math.Abs(det) > 1e-12*s0*s0*s0 {
			det0 := r0*(s2*s4-s3*s3) - s1*(r1*s4-s3*r2) + s2*(r1*s3-s2*r2)
			return det0 / det
		}
	}
	if n >= 2 {
		if den := s0*s2 - s1*s1; math.Abs(den) > 1e-12*s0*s0 {
			return (s2*r0 - s1*r1) / den
		}
	}
	return r0 / s0
}

// reportSparseDrops prints how many markers each chromosome lost to the
// sparse-kernel rule, so a whole contig cannot disappear unnoticed.
func reportSparseDrops(drops []chromDrop, total, dropped int) {
	color.Cyan("Smoothed %d variants (%d removed — effective number of markers in the kernel below %.0f)\n",
		total-dropped, dropped, minEffectiveMarkers)
	const maxLines = 10
	shown, more := 0, 0
	for _, d := range drops {
		if d.lost == 0 {
			continue
		}
		if shown == maxLines {
			more++
			continue
		}
		shown++
		if d.lost == d.total {
			color.Yellow("  %s: all %d markers removed (too sparse)", d.name, d.lost)
		} else {
			color.Yellow("  %s: %d of %d markers removed (too sparse)", d.name, d.lost, d.total)
		}
	}
	if more > 0 {
		color.Yellow("  … and %d more chromosomes/contigs with removed markers", more)
	}
}

// ─── normalisation ─────────────────────────────────────────────────────────────

// zField pairs a smoothed statistic with its Z-score field.
type zField struct {
	get func(*SmoothedStats) float64
	set func(*SmoothedStats, float64)
}

// normalise computes global robust Z-scores for every smoothed statistic. The
// Z-scores are stored in place.
func normalise(sm []SmoothedStats, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk, hmm bool) {
	var fields []zField
	// SI Z-scores only in one-bulk modes: in two-bulk modes ΔSI carries the
	// contrast and the SI Z-scores fed nothing (no CompositeZ, no detection).
	if hasHighBulk && !hasBothBulks {
		fields = append(fields, zField{func(s *SmoothedStats) float64 { return s.SmHighSI }, func(s *SmoothedStats, v float64) { s.ZHighSI = v }})
	}
	if hasLowBulk && !hasBothBulks {
		fields = append(fields, zField{func(s *SmoothedStats) float64 { return s.SmLowSI }, func(s *SmoothedStats, v float64) { s.ZLowSI = v }})
	}
	if hasBothBulks {
		fields = append(fields,
			zField{func(s *SmoothedStats) float64 { return s.SmDeltaSI }, func(s *SmoothedStats, v float64) { s.ZDeltaSI = v }},
			zField{func(s *SmoothedStats) float64 { return s.SmGstat }, func(s *SmoothedStats, v float64) { s.ZGstat = v }},
			zField{func(s *SmoothedStats) float64 { return s.SmED }, func(s *SmoothedStats, v float64) { s.ZED = v }},
		)
	}
	if hasOneBulk {
		fields = append(fields,
			zField{func(s *SmoothedStats) float64 { return s.SmAFDev }, func(s *SmoothedStats, v float64) { s.ZAFDev = v }},
			zField{func(s *SmoothedStats) float64 { return s.SmOneBulkG }, func(s *SmoothedStats, v float64) { s.ZOneBulkG = v }},
		)
	}
	if hmm {
		fields = append(fields, zField{func(s *SmoothedStats) float64 { return s.HMMLOD }, func(s *SmoothedStats, v float64) { s.ZHMMLOD = v }})
	}

	// At most normaliseWorkers fields at a time, each with one reused buffer
	// of len(sm) floats (the median and MAD sort it in place).
	sem := make(chan struct{}, normaliseWorkers)
	var wg sync.WaitGroup
	for _, f := range fields {
		wg.Add(1)
		sem <- struct{}{}
		go func() {
			defer func() { <-sem; wg.Done() }()
			robustZField(sm, f, make([]float64, len(sm)))
		}()
	}
	wg.Wait()
}

// normaliseWorkers bounds how many fields are normalised at once (memory:
// one len(sm) buffer each).
const normaliseWorkers = 3

// robustZField sets f's robust Z-score, (x − median) / (1.4826 × MAD), on
// every marker, using buf (len(sm)) as the only scratch. The 1.4826 factor
// makes MAD a consistent estimator of σ under normality (Rousseeuw & Croux
// 1993). When MAD = 0 it falls back to the classical Z (population SD), and
// to zeros for a constant track.
func robustZField(sm []SmoothedStats, f zField, buf []float64) {
	const consistencyFactor = 1.4826
	n := len(sm)
	if n == 0 {
		return
	}
	for i := range sm {
		buf[i] = f.get(&sm[i])
	}
	sort.Float64s(buf)
	med := sortedMedian(buf)
	for i := range sm {
		buf[i] = math.Abs(f.get(&sm[i]) - med)
	}
	sort.Float64s(buf)
	mad := sortedMedian(buf) * consistencyFactor
	if mad > 0 {
		for i := range sm {
			f.set(&sm[i], (f.get(&sm[i])-med)/mad)
		}
		return
	}
	// MAD = 0 fallback: classical Z-score (population SD); zeros if constant.
	mean := 0.0
	for i := range sm {
		mean += f.get(&sm[i])
	}
	mean /= float64(n)
	variance := 0.0
	for i := range sm {
		d := f.get(&sm[i]) - mean
		variance += d * d
	}
	sd := math.Sqrt(variance / float64(n))
	for i := range sm {
		z := 0.0
		if sd > 0 {
			z = (f.get(&sm[i]) - mean) / sd
		}
		f.set(&sm[i], z)
	}
}

// sortedMedian is the median of an already sorted slice.
func sortedMedian(s []float64) float64 {
	n := len(s)
	if n%2 == 1 {
		return s[n/2]
	}
	return (s[n/2-1] + s[n/2]) / 2
}

// ─── consolidation ─────────────────────────────────────────────────────────────

// consolidate fills CompositeZ (Stouffer's Z) and MaxAbsZ for every variant.
//
// CompositeZ combines k Z-scores (countZStats gives k):
//
//	Two-bulk: |ZDeltaSI| + ZGstat + ZHMMLOD   (Z(|ΔSI|) unfolded when ΔSI is unsigned)
//	One-bulk: |ZAFDev|  + ZOneBulkG + ZHMMLOD
//
// ZHMMLOD is left out (k = 2) where the HMM is not computed (2b without
// --ref-parent). ΔSI and AFDev are folded so that both tails count as signal.
// ZGstat and ZHMMLOD are robust Z-scores of non-negative statistics, so they
// are negative below the genome-wide median; CompositeZ is not bounded below.
// Direction is kept in the individual tracks (ZDeltaSI / ZAFDev). ZHighSI,
// ZLowSI and ZED are excluded (ΔSI covers them).
//
// MaxAbsZ is max |Z| over all available Z-scores, a diagnostic column.
//
// Its thresholds come from the null simulation (hmmnull.go), which runs the
// same scorer, so the combination and its null always match.
func consolidate(sm []SmoothedStats, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk, signedDSI, hmm bool) {
	for i := range sm {
		s := &sm[i]

		sum, k := 0.0, 0
		if hasBothBulks {
			// Without polarity SmDeltaSI is already smoothed |ΔSI|; its Z is
			// used as is.
			zd := s.ZDeltaSI
			if signedDSI {
				zd = math.Abs(zd)
			}
			sum += zd
			sum += s.ZGstat
			k += 2
		}
		if hasOneBulk {
			sum += math.Abs(s.ZAFDev)
			sum += s.ZOneBulkG
			k += 2
		}
		if hmm && (hasBothBulks || hasOneBulk) {
			sum += s.ZHMMLOD
			k++
		}
		if k == 0 {
			continue
		}
		s.CompositeZ = sum / math.Sqrt(float64(k))

		maxAbs := 0.0
		upd := func(z float64) {
			if a := math.Abs(z); a > maxAbs {
				maxAbs = a
			}
		}
		// SI Z-scores exist only in one-bulk modes (see normalise).
		if hasHighBulk && !hasBothBulks {
			upd(s.ZHighSI)
		}
		if hasLowBulk && !hasBothBulks {
			upd(s.ZLowSI)
		}
		if hasBothBulks {
			upd(s.ZDeltaSI)
			upd(s.ZGstat)
			upd(s.ZED)
		}
		if hasOneBulk {
			upd(s.ZAFDev)
			upd(s.ZOneBulkG)
		}
		if hmm {
			upd(s.ZHMMLOD)
		}
		s.MaxAbsZ = maxAbs
	}
}

// ─── TSV output ────────────────────────────────────────────────────────────────

// smoothLayout is a mode's track layout: which bulks, whether ΔSI is signed,
// and whether the HMM LOD is computed.
type smoothLayout struct {
	highBulk, lowBulk bool
	signedDSI         bool
	hmm               bool
}

func smoothedHeader(lay smoothLayout) []byte {
	both, one := lay.highBulk && lay.lowBulk, lay.highBulk != lay.lowBulk
	header := []string{"CHROM", "POS", "REF", "ALT", "Depth"}
	if both {
		header = append(header, "HighDepth", "LowDepth")
	}
	header = append(header, "NEff")
	if lay.highBulk {
		header = append(header, "Sm_HighSI")
		if !both {
			header = append(header, "Z_HighSI")
		}
	}
	if lay.lowBulk {
		header = append(header, "Sm_LowSI")
		if !both {
			header = append(header, "Z_LowSI")
		}
	}
	if both {
		header = append(header,
			"Sm_"+DeltaSIName(lay.signedDSI), "Z_"+DeltaSIName(lay.signedDSI),
			"Sm_Gstat", "Z_Gstat",
			"Sm_ED4", "Z_ED4",
		)
	}
	if one {
		header = append(header,
			"P0",
			"Sm_AFDev", "Z_AFDev",
			"Sm_OneBulkGstat", "Z_OneBulkGstat",
		)
	}
	if lay.hmm {
		header = append(header, "HMMLOD", "Z_HMMLOD")
	}
	header = append(header, "CompositeZ", "MaxAbsZ")
	return []byte(strings.Join(header, "\t") + "\n")
}

// appendSmoothedRow appends one marker's TSV row (floats to 6 decimals).
func appendSmoothedRow(buf []byte, s *SmoothedStats, lay smoothLayout) []byte {
	both, one := lay.highBulk && lay.lowBulk, lay.highBulk != lay.lowBulk
	f6 := func(b []byte, vs ...float64) []byte {
		for _, v := range vs {
			b = strconv.AppendFloat(append(b, '\t'), v, 'f', 6, 64)
		}
		return b
	}
	itoa := func(b []byte, v int) []byte { return strconv.AppendInt(append(b, '\t'), int64(v), 10) }

	buf = append(buf, s.CHROM...)
	buf = strconv.AppendInt(append(buf, '\t'), s.POS, 10)
	buf = append(append(buf, '\t'), s.REF...)
	buf = append(append(buf, '\t'), s.ALT...)
	buf = itoa(buf, s.Depth)
	if both {
		buf = itoa(itoa(buf, s.HighDepth), s.LowDepth)
	}
	buf = strconv.AppendFloat(append(buf, '\t'), s.NEff, 'f', 2, 64)
	if lay.highBulk {
		buf = f6(buf, s.SmHighSI)
		if !both {
			buf = f6(buf, s.ZHighSI)
		}
	}
	if lay.lowBulk {
		buf = f6(buf, s.SmLowSI)
		if !both {
			buf = f6(buf, s.ZLowSI)
		}
	}
	if both {
		buf = f6(buf, s.SmDeltaSI, s.ZDeltaSI, s.SmGstat, s.ZGstat, s.SmED, s.ZED)
	}
	if one {
		buf = f6(buf, s.OneBulkP0, s.SmAFDev, s.ZAFDev, s.SmOneBulkG, s.ZOneBulkG)
	}
	if lay.hmm {
		buf = f6(buf, s.HMMLOD, s.ZHMMLOD)
	}
	buf = f6(buf, s.CompositeZ, s.MaxAbsZ)
	return append(buf, '\n')
}

// writeSmoothedTSV writes the smoothed TSV. Every write, flush and close error
// is returned, so a full disk cannot leave a silently truncated file.
func writeSmoothedTSV(filename string, sm []SmoothedStats, lay smoothLayout) (err error) {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	w := bufio.NewWriterSize(f, 1<<20)
	if _, err := w.Write(smoothedHeader(lay)); err != nil {
		return err
	}
	buf := make([]byte, 0, 1024)
	for i := range sm {
		buf = appendSmoothedRow(buf[:0], &sm[i], lay)
		if _, err := w.Write(buf); err != nil {
			return err
		}
	}
	return w.Flush()
}

// ─── helpers ───────────────────────────────────────────────────────────────────

// BulkFlags reports which bulks a BSA mode has. Modes are matched exactly
// (not by substring); an unknown mode has no bulks.
func BulkFlags(bsaType string) (hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) {
	switch bsaType {
	case "2p2b", "hp2b", "lp2b", "2b":
		hasHighBulk, hasLowBulk = true, true
	case "2phb", "hphb", "lphb":
		hasHighBulk = true
	case "2plb", "hplb", "lplb":
		hasLowBulk = true
	}
	hasBothBulks = hasHighBulk && hasLowBulk
	hasOneBulk = hasHighBulk != hasLowBulk
	return
}
