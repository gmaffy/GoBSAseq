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
	// ReadVar0 + ρ·ReadVar1 is the read-sampling variance of the smoothed
	// Gaussian tracks (ΔSI, AFDev, ΔSI²c) under the null, for read
	// overdispersion ρ: Σ l_m²·r_m over the effective kernel weights l_m
	// (including the local linear fit near chromosome ends) and each marker's
	// beta-binomial read variance r_m = r0_m + ρ·r1_m (readVariance). BRM adds
	// the pool variance. Zero when not computed.
	ReadVar0, ReadVar1 float64

	// ── smoothed raw values ─────────────────────────────────────────────
	// Two-bulk
	SmHighSI  float64
	SmLowSI   float64
	SmDeltaSI float64
	SmGstat   float64
	SmED      float64
	// SmDeltaSI2 (2b without --ref-parent only): smoothed bias-corrected ΔSI²
	// (deltaSI2c), whose square root estimates |ΔAF| without polarity; the
	// polarity-free BRM tests it.
	SmDeltaSI2 float64
	// One-bulk
	OneBulkP0  float64 // the site's own one-bulk p0 (not smoothed; differs by polarity)
	SmAFDev    float64
	SmOneBulkG float64
	// HMMLOD is the pool-HMM LOD at the site, computed in raw stats (hmm.go):
	// two-bulk, bulks differ; one-bulk, the bulk deviates from p0. A separate
	// track in LOD units: carried through unsmoothed, not normalised, not in
	// CompositeZ. 0 where the HMM is not computed (2b without --ref-parent).
	HMMLOD float64

	// ── robust Z-scores of smoothed values ──────────────────────────────
	ZHighSI   float64
	ZLowSI    float64
	ZDeltaSI  float64
	ZGstat    float64
	ZED       float64
	ZAFDev    float64
	ZOneBulkG float64

	// ── consolidated scores ──────────────────────────────────────────────
	// CompositeZ = (|ZDeltaSI| + ZGstat) / √2 in two-bulk modes (see
	// consolidate). It is not bounded below: ZGstat is negative where G′ is
	// under the genome-wide median.
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

// scorer runs the in-memory scoring pipeline (kernel smoothing → robust Z →
// CompositeZ) on raw stats whose HMMLOD is already filled. The real data and every null-simulation replicate
// go through it, so thresholds are calibrated on exactly what is scanned.
type scorer struct {
	bw      float64
	step    float64 // smoothing grid spacing (gridStep)
	weight  func(*BSAstats) float64
	readVar func(*BSAstats) (float64, float64) // readVariance
	hmm     *hmmSetup                          // for the null simulation; nil when the HMM LOD is not computed
	lay     smoothLayout
}

func newScorer(cfg utils.AnalysisConfig, bsaType string) (*scorer, error) {
	hb, lb, _, _ := BulkFlags(bsaType)
	weight, err := markerWeight(cfg, bsaType)
	if err != nil {
		return nil, err
	}
	readVar, err := readVariance(cfg, bsaType)
	if err != nil {
		return nil, err
	}
	bw := float64(cfg.WindowSize) / 2
	sc := &scorer{bw: bw, step: gridStep(int64(cfg.StepSize), bw), weight: weight, readVar: readVar,
		lay: smoothLayout{highBulk: hb, lowBulk: lb, signedDSI: cfg.SignedDeltaSI(bsaType), hmm: HMMAvailable(cfg, bsaType)}}
	if sc.lay.hmm {
		if sc.hmm, err = newHMMSetup(cfg, bsaType); err != nil {
			return nil, err
		}
	}
	return sc, nil
}

// score returns the smoothed, normalised and consolidated markers (HMMLOD is
// carried from raw, where RawStats or the null simulation computed it). report prints the sparse-marker summary. every > 1
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
	smoothed, err := gaussianSmoothOpts(raw, smoothOpts{bw: sc.bw, step: sc.step, weight: sc.weight,
		lay: sc.lay, report: report, every: every, order: order, readVar: sc.readVar})
	if err != nil {
		return nil, err
	}
	compScales := normalise(smoothed, sc.lay)
	consolidate(smoothed, sc.lay, compScales)
	return smoothed, nil
}

// ─── smoothing ─────────────────────────────────────────────────────────────────

// markerOrder groups marker indices by chromosome (chromosomes sorted by
// name), each group sorted by position (stable, so duplicate positions keep
// input order), and picks each position's representative view. It depends
// only on CHROM, POS and depths, so it is built once per marker set and
// shared by the HMM and the smoother (and by every null replicate).
type markerOrder struct {
	chroms [][]int
	// rep[i]: raw[i] represents its position. Views sharing a position (a
	// split multi-allelic site) are different allele contrasts from the same
	// reads; only the deepest (most informative reads, HighDepth + LowDepth;
	// ties: first) enters the kernels and the HMM, so a position counts once
	// and contrasts are not mixed. Every view still gets the position's output.
	rep []bool
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
	o := &markerOrder{chroms: make([][]int, 0, len(names)), rep: make([]bool, len(raw))}
	for _, name := range names {
		idx := byChrom[name]
		sort.SliceStable(idx, func(a, b int) bool { return raw[idx[a]].POS < raw[idx[b]].POS })
		o.chroms = append(o.chroms, idx)
		for k := 0; k < len(idx); {
			best := k
			j := k + 1
			for ; j < len(idx) && raw[idx[j]].POS == raw[idx[k]].POS; j++ {
				if viewDepth(&raw[idx[j]]) > viewDepth(&raw[idx[best]]) {
					best = j
				}
			}
			o.rep[idx[best]] = true
			k = j
		}
	}
	return o
}

// viewDepth is a view's informative depth (reads on its two alleles, over
// the bulks present).
func viewDepth(s *BSAstats) int { return s.HighDepth + s.LowDepth }

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
	pop, err := utils.ParsePopulation(cfg.Population)
	if err != nil {
		return nil, err
	}
	c := pop.VarScale
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

// readVariance returns a marker's null read-sampling variance terms for the
// Gaussian BRM signal (two-bulk ΔSI, one-bulk AFDev): r0 + ρ·r1, summed over
// the bulks, with per bulk at depth d and null frequency p
//
//	r = q·(1 + (d−1)ρ)/d,  q = E[f(1−f)] = p(1−p) − Var(pool),  Var(pool) = p(1−p)/(c·n)
//
// (beta-binomial reads around the pool frequency f; c the population variance
// scale, n the bulk size). p is p0 in two-bulk modes and the site's
// bias-corrected P0 in one-bulk modes.
func readVariance(cfg utils.AnalysisConfig, bsaType string) (func(*BSAstats) (float64, float64), error) {
	high, low, one, err := cfg.BulkSizes(bsaType)
	if err != nil {
		return nil, err
	}
	pop, err := utils.ParsePopulation(cfg.Population)
	if err != nil {
		return nil, err
	}
	p0, c := pop.P0, pop.VarScale
	bulk := func(p float64, n, depth int) (r0, r1 float64) {
		d := math.Max(float64(depth), 1)
		q := p * (1 - p) * (1 - 1/(c*float64(n)))
		return q / d, q * (d - 1) / d
	}
	if _, _, both, _ := BulkFlags(bsaType); both {
		return func(s *BSAstats) (float64, float64) {
			dH, dL := s.HighDepth, s.LowDepth
			if dH <= 0 || dL <= 0 {
				dH, dL = s.Depth, s.Depth
			}
			a0, a1 := bulk(p0, high, dH)
			b0, b1 := bulk(p0, low, dL)
			return a0 + b0, a1 + b1
		}, nil
	}
	return func(s *BSAstats) (float64, float64) {
		p := s.OneBulkP0
		if !(p > 0 && p < 1) {
			p = p0
		}
		return bulk(p, one, s.Depth)
	}, nil
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

// ─── track layout ──────────────────────────────────────────────────────────────

// smoothLayout is a mode's track layout: which bulks, whether ΔSI is signed,
// and whether the HMM LOD is computed. Its tracks() table is the single
// source for what is smoothed, normalised, combined into CompositeZ and
// written to the TSV.
type smoothLayout struct {
	highBulk, lowBulk bool
	signedDSI         bool
	hmm               bool
}

func (l smoothLayout) both() bool { return l.highBulk && l.lowBulk }
func (l smoothLayout) one() bool  { return l.highBulk != l.lowBulk }

// compositeRole is a track's part in CompositeZ.
type compositeRole uint8

const (
	notInComposite  compositeRole = iota
	compositePlain                // + Z
	compositeFolded               // + |Z| (signed tracks: either direction is signal)
)

// track is one output track of a layout.
type track struct {
	smCol, zCol string                        // TSV columns ("" = none)
	sm          func(*SmoothedStats) *float64 // the smoothed (or carried) value
	z           func(*SmoothedStats) *float64 // its robust Z-score; nil = not normalised
	raw         func(*BSAstats) float64       // smoothed from this; nil = carried from the site (set in emit)
	kernel      smoothKernel
	lo, hi      float64 // clamp for the fitted value
	composite   compositeRole
	// logZ: the robust Z is taken of log(x + c), c = 10% of the track's
	// genome-wide median. For the non-negative, right-skewed tracks (G′, G1,
	// ED⁴, HMM LOD) the median/MAD Z assumes a symmetry they lack; on
	// simulated no-QTL genomes the log brings G′ and ED⁴ from Bowley skew
	// +0.36/+0.49 to +0.14/+0.17 and their largest null Z from 7.8/17.9 to
	// 2.7 (as Magwene's log G′). The offset keeps a clamped zero ~2.3 log
	// units below the median instead of stretching the lower tail.
	logZ bool
}

// tracks lists a layout's tracks in TSV column order.
//
//	two-bulk: HighSI, LowSI (smoothed, not normalised: ΔSI carries the
//	          contrast), ΔSI (|ΔSI| without polarity), G′, ED⁴
//	one-bulk: SI (smoothed and normalised), P0 (carried), AFDev, G1
//	both:     HMM LOD (carried: already linkage-aware) when computed
//
// CompositeZ = Σ of the composite tracks' Z (|Z| when folded) / √k.
func (l smoothLayout) tracks() []track {
	inf := math.Inf(1)
	var t []track
	si := func(name string, raw func(*BSAstats) float64, sm, z func(*SmoothedStats) *float64) {
		tr := track{smCol: "Sm_" + name, sm: sm, raw: raw, kernel: kGauss, lo: 0, hi: 1}
		if !l.both() { // SI Z-scores only in one-bulk modes
			tr.zCol, tr.z = "Z_"+name, z
		}
		t = append(t, tr)
	}
	if l.highBulk {
		si("HighSI", func(s *BSAstats) float64 { return s.HighSI },
			func(s *SmoothedStats) *float64 { return &s.SmHighSI }, func(s *SmoothedStats) *float64 { return &s.ZHighSI })
	}
	if l.lowBulk {
		si("LowSI", func(s *BSAstats) float64 { return s.LowSI },
			func(s *SmoothedStats) *float64 { return &s.SmLowSI }, func(s *SmoothedStats) *float64 { return &s.ZLowSI })
	}
	if l.both() {
		name, lo, role := DeltaSIName(l.signedDSI), -1.0, compositeFolded
		if !l.signedDSI { // |ΔSI|: already folded
			lo, role = 0, compositePlain
		}
		t = append(t,
			track{"Sm_" + name, "Z_" + name, func(s *SmoothedStats) *float64 { return &s.SmDeltaSI }, func(s *SmoothedStats) *float64 { return &s.ZDeltaSI },
				func(s *BSAstats) float64 { return s.DeltaSI }, kGauss, lo, 1, role, false},
			track{"Sm_Gstat", "Z_Gstat", func(s *SmoothedStats) *float64 { return &s.SmGstat }, func(s *SmoothedStats) *float64 { return &s.ZGstat },
				func(s *BSAstats) float64 { return s.Gstat }, kTricube, 0, inf, compositePlain, true},
			track{"Sm_ED4", "Z_ED4", func(s *SmoothedStats) *float64 { return &s.SmED }, func(s *SmoothedStats) *float64 { return &s.ZED },
				func(s *BSAstats) float64 { return s.ED }, kLoess, 0, inf, notInComposite, true},
		)
		if !l.signedDSI { // polarity-free BRM input
			t = append(t, track{smCol: "Sm_DeltaSI2c", sm: func(s *SmoothedStats) *float64 { return &s.SmDeltaSI2 },
				raw: deltaSI2c, kernel: kGauss, lo: -1, hi: 1})
		}
	}
	if l.one() {
		t = append(t,
			track{smCol: "P0", sm: func(s *SmoothedStats) *float64 { return &s.OneBulkP0 }},
			track{"Sm_AFDev", "Z_AFDev", func(s *SmoothedStats) *float64 { return &s.SmAFDev }, func(s *SmoothedStats) *float64 { return &s.ZAFDev },
				func(s *BSAstats) float64 { return s.OneBulkAFDev }, kGauss, -1, 1, compositeFolded, false},
			track{"Sm_OneBulkGstat", "Z_OneBulkGstat", func(s *SmoothedStats) *float64 { return &s.SmOneBulkG }, func(s *SmoothedStats) *float64 { return &s.ZOneBulkG },
				func(s *BSAstats) float64 { return s.OneBulkGstat }, kTricube, 0, inf, compositePlain, true},
		)
	}
	if l.hmm && (l.both() || l.one()) { // carried from raw stats; its own track only
		t = append(t, track{smCol: "HMMLOD", sm: func(s *SmoothedStats) *float64 { return &s.HMMLOD }})
	}
	return t
}

// oneSided reports whether only high values of a track are evidence: the
// non-negative statistics (log-Z tracks) and |ΔSI|.
func (tr track) oneSided() bool { return tr.logZ || tr.smCol == "Sm_AbsDeltaSI" }

// deltaSI2c is a site's bias-corrected squared ΔSI: (SI_H − SI_L)² minus the
// unbiased read-sampling variance SI(1−SI)/(d−1) of each bulk. ΔSI² does not
// depend on which allele is REF, and E[deltaSI2c] = (pool_H − pool_L)², so
// its smoothed value estimates the squared pool frequency difference without
// polarity (2b without --ref-parent).
func deltaSI2c(s *BSAstats) float64 {
	d := s.HighSI - s.LowSI
	v := d * d
	if s.HighDepth > 1 {
		v -= s.HighSI * (1 - s.HighSI) / float64(s.HighDepth-1)
	}
	if s.LowDepth > 1 {
		v -= s.LowSI * (1 - s.LowSI) / float64(s.LowDepth-1)
	}
	return v
}

// compositeCount is k, the number of Z-scores in CompositeZ.
func (l smoothLayout) compositeCount() int {
	k := 0
	for _, tr := range l.tracks() {
		if tr.composite != notInComposite {
			k++
		}
	}
	return k
}

// smoothField pairs a raw statistic with the smoothed field it fills; lo/hi
// clamp the fitted value to the statistic's range (local regressions can
// overshoot at chromosome ends).
type smoothField struct {
	get    func(*BSAstats) float64
	set    func(*SmoothedStats, float64)
	kernel smoothKernel
	lo, hi float64
}

// smoothFields lists the layout's smoothed tracks, grouped by kernel:
// [Gaussian | tricube | loess].
func smoothFields(l smoothLayout) []smoothField {
	var byKernel [3][]smoothField
	for _, tr := range l.tracks() {
		if tr.raw == nil {
			continue
		}
		sm := tr.sm
		byKernel[tr.kernel] = append(byKernel[tr.kernel], smoothField{tr.raw, func(s *SmoothedStats, v float64) { *sm(s) = v }, tr.kernel, tr.lo, tr.hi})
	}
	return append(append(byKernel[kGauss], byKernel[kTricube]...), byKernel[kLoess]...)
}

// ─── smoothing ─────────────────────────────────────────────────────────────────

// smoothOpts configures gaussianSmooth.
type smoothOpts struct {
	bw     float64                 // Gaussian σ in bp (WindowSize/2); tricube half-width h = 2·bw
	step   float64                 // grid spacing for dense chromosomes; 0 = always exact at markers
	weight func(*BSAstats) float64 // marker weight (markerWeight)
	lay    smoothLayout
	report bool         // print the sparse-marker summary
	every  int          // > 1: output only every n-th marker of each chromosome (compact, chromosome order)
	order  *markerOrder // precomputed marker order; nil = build it
	// readVar, if set, gives a marker's null read variance terms (r0, r1:
	// variance r0 + ρ·r1); the Gaussian pass then fills ReadVar0/ReadVar1.
	readVar func(*BSAstats) (r0, r1 float64)
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
	return gaussianSmoothOpts(raw, smoothOpts{bw: bw, weight: weight,
		lay: smoothLayout{highBulk: hasHighBulk, lowBulk: hasLowBulk, signedDSI: true}, report: report, every: 1})
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
// the next. Within a chromosome, query points are smoothed in parallel chunks
// (kernelSpec.eval); all kernels share one pass over the neighbours, and each
// sum runs over the neighbours in position order, so results do not depend on
// scheduling.
//
// The output keeps the input (raw) order. With every > 1 (null simulation) it
// is compact (only the output markers, no len(raw) array) and in chromosome
// order.
func gaussianSmoothOpts(raw []BSAstats, o smoothOpts) ([]SmoothedStats, error) {
	k := newKernelSpec(o.bw, smoothFields(o.lay))
	order := o.order
	if order == nil {
		order = newMarkerOrder(raw)
	}
	o.order = order // newChromWork reads the representative views from it
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

	// One persistent worker pool; chromosomes are fed one at a time.
	type job struct {
		q0, q1 int
		c      *chromWork
		done   *sync.WaitGroup
	}
	workers := max(runtime.GOMAXPROCS(0), 1)
	jobs := make(chan job, 64)
	var pool sync.WaitGroup
	for range workers {
		pool.Add(1)
		go func() {
			defer pool.Done()
			acc := k.newAcc()
			for jb := range jobs {
				c := jb.c
				edgeFit := c.span() >= 2*k.cutoff // local linear near ends only on chromosomes ≥ 6σ
				w := k.startWindow(c.npos, c.q[jb.q0])
				for qi := jb.q0; qi < jb.q1; qi++ {
					w.advance(c.npos, c.q[qi], k.cutoff, k.h)
					k.eval(c, qi, w, acc, edgeFit)
				}
				jb.done.Done()
			}
		}()
	}

	for _, idx := range order.chroms {
		c := newChromWork(raw, idx, k.fields, o)
		var done sync.WaitGroup
		chunk := min(max(len(c.q)/(4*workers), 16), smoothChunk)
		for q0 := 0; q0 < len(c.q); q0 += chunk {
			done.Add(1)
			jobs <- job{q0, min(q0+chunk, len(c.q)), c, &done}
		}
		done.Wait()
		drops = append(drops, c.emit(raw, &out, keep, k.fields, o))
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

// kernelSpec is a smoothing pass's fixed parameters: bandwidths and the
// fields, grouped by kernel (ng Gaussian, nt tricube, nl loess).
type kernelSpec struct {
	bw, invBw, h, cutoff float64
	fields               []smoothField
	nf, ng, nt, nl       int
}

func newKernelSpec(bw float64, fields []smoothField) *kernelSpec {
	k := &kernelSpec{bw: bw, invBw: 1 / bw, h: 2 * bw, cutoff: kernelCutoffSigmas * bw, fields: fields, nf: len(fields)}
	for _, fd := range fields {
		switch fd.kernel {
		case kGauss:
			k.ng++
		case kTricube:
			k.nt++
		}
	}
	k.nl = k.nf - k.ng - k.nt
	return k
}

// kernelAcc is one worker's per-query accumulators.
type kernelAcc struct {
	g, t, l []float64 // Σw·x, Σw·d·x per Gaussian field; Σt·x per tricube field; Σt·x, Σt·u·x, Σt·u²·x per loess field
}

func (k *kernelSpec) newAcc() *kernelAcc {
	return &kernelAcc{g: make([]float64, 2*k.ng), t: make([]float64, k.nt), l: make([]float64, 3*k.nl)}
}

// window is the neighbour index ranges of a query: [lo, hi) within the ±3σ
// Gaussian window, [tlo, thi) within the open ±h tricube window.
type window struct{ lo, hi, tlo, thi int }

// startWindow positions a window for query q by binary search.
func (k *kernelSpec) startWindow(pos []float64, q float64) window {
	n := len(pos)
	return window{
		lo:  sort.Search(n, func(m int) bool { return pos[m] >= q-k.cutoff }),
		hi:  sort.Search(n, func(m int) bool { return pos[m] > q+k.cutoff }),
		tlo: sort.Search(n, func(m int) bool { return pos[m] > q-k.h }),
		thi: sort.Search(n, func(m int) bool { return pos[m] >= q+k.h }),
	}
}

// advance moves the window forward to query q (queries are sorted).
func (w *window) advance(pos []float64, q, cutoff, h float64) {
	n := len(pos)
	for w.lo < n && pos[w.lo] < q-cutoff {
		w.lo++
	}
	for w.hi < n && pos[w.hi] <= q+cutoff {
		w.hi++
	}
	for w.tlo < n && pos[w.tlo] <= q-h { // tricube window: |Δ| < h
		w.tlo++
	}
	for w.thi < n && pos[w.thi] < q+h {
		w.thi++
	}
}

// eval smooths every field at query qi of c over window w, writing c.res and
// c.neff. edgeFit enables the local linear fit near chromosome ends.
func (k *kernelSpec) eval(c *chromWork, qi int, w window, acc *kernelAcc, edgeFit bool) {
	pos, pK, nf, ng, nt, nl := c.npos, c.q[qi], k.nf, k.ng, k.nt, k.nl
	useTri := nt+nl > 0
	clear(acc.g)
	clear(acc.t)
	clear(acc.l)
	var g0, g1, g2, w2 float64 // Gaussian moments Σw·dⁱ, Σw²
	var va, vb [3]float64      // read variance moments Σw²·r·dⁱ for r0 (va) and r1 (vb)
	rvOn := c.nrv != nil
	var s0, s1, s2, s3, s4 float64 // tricube moments Σt·uⁱ
	for m := w.lo; m < w.hi; m++ {
		d := (pos[m] - pK) * k.invBw
		kg := math.Exp(-0.5 * d * d)
		wm := c.nwt[m] * kg // marker weight × Gaussian
		wd := wm * d
		g0 += wm
		g1 += wd
		g2 += wd * d
		if c.nwt2 == nil {
			w2 += wm * wm
		} else { // aggregated: Σ(wt·K)² = K²·Σwt²
			w2 += kg * kg * c.nwt2[m]
		}
		if rvOn { // Σ(wt·K)²·r = K²·Σwt²·r
			k2 := kg * kg
			a0, b0 := k2*c.nrv[2*m], k2*c.nrv[2*m+1]
			va[0] += a0
			va[1] += a0 * d
			va[2] += a0 * d * d
			vb[0] += b0
			vb[1] += b0 * d
			vb[2] += b0 * d * d
		}
		row := c.nvals[m*nf : (m+1)*nf : (m+1)*nf]
		for f := range ng {
			acc.g[2*f] += wm * row[f]
			acc.g[2*f+1] += wd * row[f]
		}
		if m < w.tlo || m >= w.thi || !useTri {
			continue
		}
		u := 0.5 * d // Δ/h, since h = 2·bw
		a := math.Abs(u)
		q := 1 - a*a*a
		tw := c.nwt[m] * q * q * q // marker weight × tricube
		s0 += tw
		for f := range nt {
			acc.t[f] += tw * row[ng+f]
		}
		if nl > 0 {
			u2 := u * u
			s1 += tw * u
			s2 += tw * u2
			s3 += tw * u2 * u
			s4 += tw * u2 * u2
			for f := range nl {
				y := tw * row[ng+nt+f]
				acc.l[3*f] += y
				acc.l[3*f+1] += y * u
				acc.l[3*f+2] += y * u2
			}
		}
	}
	c.neff[qi] = 0
	if w2 > 0 {
		c.neff[qi] = g0 * g0 / w2
	}
	alpha := 0.0
	if edgeFit {
		alpha = edgeBlend(math.Min(pK-c.pos[0], c.pos[len(c.pos)-1]-pK), k.bw)
	}
	if rvOn {
		// The Gaussian estimate is Σ l_m·x_m with effective weights
		// l_m = w_m·(a + b·d_m): a = 1/g0, b = 0 for the local average;
		// blended with the local linear fit β0 (weights w·(g2 − g1·d)/D),
		// a = (1−α)/g0 + α·g2/D, b = −α·g1/D. Its variance is
		// Σ w²·r·(a + b·d)² = a²·M0 + 2ab·M1 + b²·M2.
		a, b := 0.0, 0.0
		if g0 > 0 {
			a = 1 / g0
			if den := g0*g2 - g1*g1; alpha > 0 && w.hi-w.lo >= 2 && den > 1e-12*g0*g2 {
				a = (1-alpha)/g0 + alpha*g2/den
				b = -alpha * g1 / den
			}
		}
		c.rv[2*qi] = a*a*va[0] + 2*a*b*va[1] + b*b*va[2]
		c.rv[2*qi+1] = a*a*vb[0] + 2*a*b*vb[1] + b*b*vb[2]
	}
	res := c.res[qi*nf : (qi+1)*nf]
	for f := range nf {
		var v float64
		switch {
		case f < ng:
			mean := 0.0
			if g0 > 0 {
				mean = acc.g[2*f] / g0
			}
			v = mean
			if alpha > 0 {
				v += alpha * (localLinear(g0, g1, g2, acc.g[2*f], acc.g[2*f+1], w.hi-w.lo) - mean)
			}
		case f < ng+nt:
			if s0 > 0 {
				v = acc.t[f-ng] / s0
			}
		default:
			l := f - ng - nt
			v = loessFit(s0, s1, s2, s3, s4, acc.l[3*l], acc.l[3*l+1], acc.l[3*l+2], w.thi-w.tlo)
		}
		res[f] = math.Min(math.Max(v, k.fields[f].lo), k.fields[f].hi)
	}
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
	nrv   []float64 // Σ marker weight²·(r0, r1) per neighbour, row-major [neighbour][2]; nil if not computed
	q     []float64 // query positions, sorted (markers, or grid points plus sparse markers)
	grid  bool      // grid mode: dense markers are interpolated between grid points
	gq    []int     // grid mode: q-index of each grid point
	mq    []int     // grid mode: q-index of a marker evaluated exactly, or −1
	res   []float64 // results, row-major [query][field]
	neff  []float64 // NEff per query
	rv    []float64 // read variance terms (ReadVar0, ReadVar1) per query, row-major [query][2]; nil if not computed
}

// sparseIntervalMarkers: in grid mode, markers in a grid interval holding at
// most this many markers are evaluated exactly. In sparse stretches the
// smoothed curve changes quickly as single markers enter or leave the kernel,
// so interpolation would be inaccurate there; few markers means little cost.
const sparseIntervalMarkers = 4

// span is the distance from the chromosome's first to its last marker.
func (c *chromWork) span() float64 { return c.pos[len(c.pos)-1] - c.pos[0] }

func newChromWork(raw []BSAstats, idx []int, fields []smoothField, o smoothOpts) *chromWork {
	n, nf := len(idx), len(fields)
	c := &chromWork{name: raw[idx[0]].CHROM, idx: idx, pos: make([]float64, n)}
	wt, vals := make([]float64, n), make([]float64, n*nf)
	var rv []float64 // per marker wt²·(r0, r1)
	if o.readVar != nil {
		rv = make([]float64, 2*n)
	}
	for k, ri := range idx {
		s := &raw[ri]
		c.pos[k] = float64(s.POS)
		wt[k] = o.weight(s)
		for f, fd := range fields {
			vals[k*nf+f] = fd.get(s)
		}
		if rv != nil {
			r0, r1 := o.readVar(s)
			rv[2*k], rv[2*k+1] = wt[k]*wt[k]*r0, wt[k]*wt[k]*r1
		}
	}
	c.q, c.npos, c.nwt, c.nvals, c.nrv = c.pos, c.pos, wt, vals, rv

	// Grid mode for dense chromosomes. A grid needs a positive span: markers
	// all at one position (e.g. a multi-allelic site split on a tiny contig)
	// stay on the exact path.
	first, last := c.pos[0], c.pos[n-1]
	grid := o.step > 0 && last > first && int((last-first)/o.step)+2 < n
	dup := false
	for k := 1; k < n && !dup; k++ {
		dup = c.pos[k] == c.pos[k-1]
	}
	// Markers sharing a position (views of one split multi-allelic site) come
	// from the same reads and plants, so each position is one kernel
	// neighbour ("unit"): the representative (deepest) view's weight and
	// values (markerOrder.rep). Units are only built when needed: duplicates
	// present, or grid mode (aggregation).
	if dup || grid {
		isRep := make([]bool, n)
		for k, ri := range idx {
			isRep[k] = o.order.rep[ri]
		}
		u := newPosUnits(c.pos, wt, vals, rv, nf, isRep)
		if dup {
			c.npos, c.nwt, c.nwt2, c.nvals, c.nrv = u.pos, u.w, u.w2, u.means(), u.rv
		}
		if grid {
			c.buildGrid(first, last, o.step)
			c.aggregate(u, nf, o.bw/aggregateBinSigmas)
		}
	}
	c.res = make([]float64, len(c.q)*nf)
	c.neff = make([]float64, len(c.q))
	if c.nrv != nil {
		c.rv = make([]float64, 2*len(c.q))
	}
	return c
}

// posUnits are a chromosome's distinct marker positions: per position, the
// representative view's weight (w), w² (w2), weighted values (sums, w·value,
// row-major [unit][field]) and the index of the position's first marker.
type posUnits struct {
	pos, w, w2, sums []float64
	rv               []float64 // representative's wt²·(r0, r1), [unit][2]; nil if not computed
	first            []int
	nf               int
}

func newPosUnits(pos, wt, vals, rv []float64, nf int, isRep []bool) *posUnits {
	u := &posUnits{nf: nf}
	for k := 0; k < len(pos); {
		j, r := k, k
		for j < len(pos) && pos[j] == pos[k] {
			if isRep[j] {
				r = j
			}
			j++
		}
		w := wt[r]
		u.pos = append(u.pos, pos[k])
		u.w = append(u.w, w)
		u.w2 = append(u.w2, w*w)
		u.first = append(u.first, k)
		for f := range nf {
			u.sums = append(u.sums, w*vals[r*nf+f])
		}
		if rv != nil {
			u.rv = append(u.rv, rv[2*r], rv[2*r+1])
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
	c.npos, c.nwt, c.nwt2, c.nvals, c.nrv = nil, nil, nil, nil, nil
	sums := make([]float64, nf)
	var rv0, rv1 float64 // Σ wt²·(r0, r1) of the point
	add := func(p, w, w2 float64) {
		c.npos = append(c.npos, p)
		c.nwt = append(c.nwt, w)
		c.nwt2 = append(c.nwt2, w2)
		if u.rv != nil {
			c.nrv = append(c.nrv, rv0, rv1)
		}
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
			if u.rv != nil {
				rv0, rv1 = u.rv[2*k], u.rv[2*k+1]
			}
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
		rv0, rv1 = 0, 0
		j := k
		for ; j < u.n() && c.mq[u.first[j]] < 0 && u.pos[j] < hi; j++ {
			w += u.w[j]
			w2 += u.w2[j]
			if u.rv != nil {
				rv0 += u.rv[2*j]
				rv1 += u.rv[2*j+1]
			}
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
	hasOneBulk := o.lay.one()
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
		var nEff, rv0, rv1 float64
		val := func(f int) float64 { return c.res[qk*nf+f] }
		if qk >= 0 {
			nEff = c.neff[qk]
			if c.rv != nil {
				rv0, rv1 = c.rv[2*qk], c.rv[2*qk+1]
			}
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
			if c.rv != nil {
				rv0 = c.rv[2*a] + t*(c.rv[2*b]-c.rv[2*a])
				rv1 = c.rv[2*a+1] + t*(c.rv[2*b+1]-c.rv[2*a+1])
			}
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
		s.ReadVar0, s.ReadVar1 = rv0, rv1
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

// zField pairs a smoothed statistic with its Z-score field; logZ: take the Z
// of log(x + c), c = 10% of the genome-wide median (see track.logZ).
type zField struct {
	get  func(*SmoothedStats) float64
	set  func(*SmoothedStats, float64)
	logZ bool
}

// zScale is a robust Z transform: z = (x − center) / scale (0 if scale = 0).
type zScale struct {
	center, scale float64
	set           bool
}

func (z zScale) of(x float64) float64 {
	if z.scale > 0 {
		return (x - z.center) / z.scale
	}
	return 0
}

// normalise computes the genome-wide robust Z-score of every track that has
// one (smoothLayout.tracks), in place. It returns, per track index, the
// untransformed robust Z scale of the log-Z tracks that enter CompositeZ:
// the Z columns of those tracks are log-based (symmetric display), while
// CompositeZ uses the untransformed robust Z (see consolidate).
func normalise(sm []SmoothedStats, lay smoothLayout) []zScale {
	tracks := lay.tracks()
	raw := make([]zScale, len(tracks))
	type job struct {
		f   zField
		dst *zScale // raw-scale job: store the scale instead of setting Z
	}
	var jobs []job
	for i, tr := range tracks {
		if tr.z == nil {
			continue
		}
		get, set := tr.sm, tr.z
		jobs = append(jobs, job{f: zField{func(s *SmoothedStats) float64 { return *get(s) }, func(s *SmoothedStats, v float64) { *set(s) = v }, tr.logZ}})
		if tr.logZ && tr.composite != notInComposite {
			jobs = append(jobs, job{f: zField{get: func(s *SmoothedStats) float64 { return *get(s) }}, dst: &raw[i]})
		}
	}

	// At most normaliseWorkers jobs at a time, each with one reused buffer
	// of len(sm) floats (the median and MAD sort it in place).
	sem := make(chan struct{}, normaliseWorkers)
	var wg sync.WaitGroup
	for _, jb := range jobs {
		wg.Add(1)
		sem <- struct{}{}
		go func() {
			defer func() { <-sem; wg.Done() }()
			buf := make([]float64, len(sm))
			if jb.dst != nil {
				*jb.dst = robustZScale(sm, jb.f.get, buf)
				jb.dst.set = true
				return
			}
			robustZField(sm, jb.f, buf)
		}()
	}
	wg.Wait()
	return raw
}

// normaliseWorkers bounds how many fields are normalised at once (memory:
// one len(sm) buffer each).
const normaliseWorkers = 3

// robustZScale returns the robust Z transform of get over sm: center = median,
// scale = 1.4826 × MAD (a consistent estimator of σ under normality; Rousseeuw
// & Croux 1993); when MAD = 0, the classical mean and population SD. buf
// (len(sm)) is the only scratch.
func robustZScale(sm []SmoothedStats, get func(*SmoothedStats) float64, buf []float64) zScale {
	const consistencyFactor = 1.4826
	n := len(sm)
	if n == 0 {
		return zScale{}
	}
	for i := range sm {
		buf[i] = get(&sm[i])
	}
	sort.Float64s(buf)
	med := sortedMedian(buf)
	for i := range sm {
		buf[i] = math.Abs(get(&sm[i]) - med)
	}
	sort.Float64s(buf)
	if mad := sortedMedian(buf) * consistencyFactor; mad > 0 {
		return zScale{center: med, scale: mad}
	}
	mean := 0.0
	for i := range sm {
		mean += get(&sm[i])
	}
	mean /= float64(n)
	variance := 0.0
	for i := range sm {
		d := get(&sm[i]) - mean
		variance += d * d
	}
	return zScale{center: mean, scale: math.Sqrt(variance / float64(n))}
}

// robustZField sets f's robust Z-score (robustZScale) on every marker, using
// buf (len(sm)) as the only scratch; zeros for a constant track. With f.logZ
// the Z is of log(x + c), c = 10% of the median (see track.logZ).
func robustZField(sm []SmoothedStats, f zField, buf []float64) {
	if len(sm) == 0 {
		return
	}
	get := f.get
	if f.logZ {
		for i := range sm {
			buf[i] = f.get(&sm[i])
		}
		sort.Float64s(buf)
		c := 0.1 * sortedMedian(buf)
		if !(c > 0) {
			c = 1e-12
		}
		get = func(s *SmoothedStats) float64 { return math.Log(math.Max(f.get(s), 0) + c) }
	}
	z := robustZScale(sm, get, buf)
	for i := range sm {
		f.set(&sm[i], z.of(get(&sm[i])))
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

// consolidate fills CompositeZ (Stouffer's Z) and MaxAbsZ for every variant,
// from the layout's tracks:
//
//	Two-bulk: (|ZDeltaSI| + ZGstat) / √2   (Z(|ΔSI|) unfolded when ΔSI is unsigned)
//	One-bulk: (|ZAFDev|  + ZOneBulkG) / √2
//
// The HMM LOD is not part of it: it is a separate track in LOD units with its
// own genome-wide thresholds. ΔSI and AFDev are folded so that both tails
// count as signal. ZGstat is a robust Z of a non-negative statistic, so it is
// negative below the genome-wide median; CompositeZ is not bounded below.
// Direction is kept in the individual tracks (ZDeltaSI / ZAFDev).
//
// MaxAbsZ is the largest deviation over all normalised tracks, a diagnostic
// column: |Z| for signed tracks (ΔSI, AFDev, one-bulk SI), and the positive
// part of Z for one-sided tracks (G′, G1, ED⁴, |ΔSI|), where a low value is
// not evidence.
//
// For the log-Z tracks in CompositeZ (G′, G1) it uses their untransformed
// robust Z (raw, from normalise), not the log-based Z column: on overdispersed
// data the log-based components raised CompositeZ's pointwise false-positive
// rate (1.4% vs 1.0% at p99). raw = nil uses the Z columns as they are.
//
// Its thresholds come from the null simulation (hmmnull.go), which runs the
// same scorer, so the combination and its null always match.
func consolidate(sm []SmoothedStats, lay smoothLayout, raw []zScale) {
	tracks := lay.tracks()
	k := lay.compositeCount()
	for i := range sm {
		s := &sm[i]
		if k == 0 {
			continue
		}
		sum, maxAbs := 0.0, 0.0
		for ti, tr := range tracks {
			if tr.z == nil {
				continue
			}
			z := *tr.z(s)
			zc := z // CompositeZ component: untransformed Z for log-Z tracks
			if ti < len(raw) && raw[ti].set {
				zc = raw[ti].of(*tr.sm(s))
			}
			switch tr.composite {
			case compositePlain:
				sum += zc
			case compositeFolded:
				sum += math.Abs(zc)
			}
			a := math.Abs(z)
			if tr.oneSided() { // only high values are evidence
				a = math.Max(z, 0)
			}
			if a > maxAbs {
				maxAbs = a
			}
		}
		s.CompositeZ = sum / math.Sqrt(float64(k))
		s.MaxAbsZ = maxAbs
	}
}

// ─── TSV output ────────────────────────────────────────────────────────────────

func smoothedHeader(lay smoothLayout) []byte {
	header := []string{"CHROM", "POS", "REF", "ALT", "Depth"}
	if lay.both() {
		header = append(header, "HighDepth", "LowDepth")
	}
	header = append(header, "NEff")
	for _, tr := range lay.tracks() {
		for _, col := range []string{tr.smCol, tr.zCol} {
			if col != "" {
				header = append(header, col)
			}
		}
	}
	header = append(header, "CompositeZ", "MaxAbsZ")
	return []byte(strings.Join(header, "\t") + "\n")
}

// appendSmoothedRow appends one marker's TSV row (floats to 6 decimals,
// NEff to 2). tracks is lay.tracks().
func appendSmoothedRow(buf []byte, s *SmoothedStats, lay smoothLayout, tracks []track) []byte {
	f6 := func(b []byte, v float64) []byte { return strconv.AppendFloat(append(b, '\t'), v, 'f', 6, 64) }
	itoa := func(b []byte, v int) []byte { return strconv.AppendInt(append(b, '\t'), int64(v), 10) }

	buf = append(buf, s.CHROM...)
	buf = strconv.AppendInt(append(buf, '\t'), s.POS, 10)
	buf = append(append(buf, '\t'), s.REF...)
	buf = append(append(buf, '\t'), s.ALT...)
	buf = itoa(buf, s.Depth)
	if lay.both() {
		buf = itoa(itoa(buf, s.HighDepth), s.LowDepth)
	}
	buf = strconv.AppendFloat(append(buf, '\t'), s.NEff, 'f', 2, 64)
	for _, tr := range tracks {
		if tr.smCol != "" {
			buf = f6(buf, *tr.sm(s))
		}
		if tr.zCol != "" {
			buf = f6(buf, *tr.z(s))
		}
	}
	buf = f6(f6(buf, s.CompositeZ), s.MaxAbsZ)
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
	tracks := lay.tracks()
	for i := range sm {
		buf = appendSmoothedRow(buf[:0], &sm[i], lay, tracks)
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
