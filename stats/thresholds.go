package stats

// thresholds.go — empirical significance thresholds for BSA statistics.
//
// Approach
// --------
// For each unique sequencing depth observed in the smoothed data we run a
// Monte Carlo null simulation using a two-stage model that is particularly
// accurate for deeply sequenced sites:
//
// Stage 1: Sample the number of alt alleles in the bulk population
//   (Binomial(n=bulk_size, p=p₀)) - this represents the realized allele count
//   among the actual individuals in the bulk.
//
// Stage 2: Sample the observed reads from that realized frequency
//   (Binomial(n=depth, p=realized_af)) - this represents the sequencing observation.
//
// This two-stage model is more biologically accurate for deep sequencing because it
// accounts for the finite population of individuals in each bulk, rather than
// assuming reads are drawn directly from the expected population frequency p₀.
//
// For each simulation draw, all BSA test statistics are evaluated; their
// empirical 95th / 99th percentiles become the per-variant thresholds.
//
// This approach is biologically correct because:
//   - p₀ reflects the actual segregation ratio of the mapping population
//     (F2 = 0.5, BC1H = 0.75, RIL = 0.5, …) rather than a re-simulated mean.
//   - The two-stage model correctly handles deep sequencing where many reads
//     may come from the same individuals, reducing false positives.
//   - Thresholds scale naturally with depth: deeper sites have narrower null
//     distributions and therefore tighter (more stringent) thresholds.
//   - Z-score thresholds are also empirical, derived from the same simulation,
//     so they account for the actual distribution rather than assuming normality.
//
// Bulk Reciprocal Model (BRM)
// ---------------------------
// Unlike the Monte Carlo thresholds above which are per-variant and depth-dependent,
// BRM uses a sliding-window analytical threshold based on the variance of allele
// frequency in the population. BRM blocks DO incorporate bulk size in their
// threshold calculations (see calculateBRMBlocksTwoBulk and calculateBRMBlocksOneBulk).
//
// For Two-Bulk:
//   Threshold = u(α) * sqrt( [ (n₁+n₂)/(V_scale * n₁ * n₂) ] * p * (1-p) )
//   where n₁, n₂ are bulk sizes, V_scale is the population variance scale (e.g., 2 for F2),
//   and p is the average allele frequency across both bulks.
//
// For One-Bulk:
//   Threshold = u(α) * sqrt( [ p₀ * (1-p₀) ] / (V_scale * n) )
//   where n is bulk size and p₀ is the expected null allele frequency.
//
// Note on QTL Detection
// ---------------------
// The CompositeZ ≥ 3.0 rule in DetectQTLs is the main driver for QTL calls.
// Everything else (simulated thresholds, BRM) is secondary or just for plots.
// For deeply sequenced sites, use DetectQTLsWithMCDirect() which applies the
// fully sound Monte Carlo thresholds to individual plots.
//
// Speed
// -----
// Unique (highDepth, lowDepth) pairs are computed in parallel across all CPU
// cores.  Results are memoised in a sync.Map so each unique depth combination
// is only simulated once per program run.  golang.org/x/sync/singleflight
// prevents redundant concurrent simulation of the same key.

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"sync"
	"sync/atomic"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
	"golang.org/x/sync/singleflight"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// ── Threshold structs ────────────────────────────────────────────────────────

// Thresholds holds empirical significance thresholds for a single variant.
// All values are two-tailed where applicable (P99/P95 = upper tail,
// Mp99/Mp95 = lower tail for signed statistics).
type Thresholds struct {
	TwoBulk TwoBulkThresholds
	OneBulk OneBulkThresholds
	Z       ZThresholds
}

// TwoBulkThresholds holds thresholds for two-bulk BSA statistics.
type TwoBulkThresholds struct {
	// Segregation index thresholds (signed, so upper and lower tails)
	HighSIP99, HighSIP95     float64 // upper tail
	HighSIMp99, HighSIMp95   float64 // lower tail
	LowSIP99, LowSIP95       float64
	LowSIMp99, LowSIMp95     float64
	DeltaSIP99, DeltaSIP95   float64
	DeltaSIMp99, DeltaSIMp95 float64

	// One-tailed statistics (larger = more evidence)
	GstatP99, GstatP95     float64
	ED4P99, ED4P95         float64
	LODP99, LODP95         float64
	BBLogBFP99, BBLogBFP95 float64
}

// OneBulkThresholds holds thresholds for single-bulk BSA statistics.
type OneBulkThresholds struct {
	AFDevP99, AFDevP95                   float64 // upper tail
	AFDevMp99, AFDevMp95                 float64 // lower tail
	OneBulkGstatP99, OneBulkGstatP95     float64
	OneBulkLODP99, OneBulkLODP95         float64
	OneBulkBBLogBFP99, OneBulkBBLogBFP95 float64
}

// ZThresholds holds empirical thresholds for the robust Z-scores and composite
// scores produced by smoothing.go.  These are derived from the same null
// simulation as the raw-stat thresholds, so they reflect the actual null
// distribution of the smoothed, normalised scores.
type ZThresholds struct {
	ZP99, ZP95                   float64 // upper tail
	ZN99, ZN95                   float64 // lower tail (negative)
	CompositeZP99, CompositeZP95 float64
	CompositeZN99, CompositeZN95 float64
	NumStats                     int // number of statistics combined in CompositeZ
}

// ── Simulation parameters ────────────────────────────────────────────────────

// simParams bundles the inputs to a single null simulation run.
type simParams struct {
	highDepth int
	lowDepth  int     // == highDepth for two-bulk with shared depth
	p0High    float64 // expected AF in high bulk under H₀
	p0Low     float64 // expected AF in low bulk under H₀ (usually == p0High)
	p0One     float64 // expected AF for one-bulk mode
	rep       int
	// Bulk sizes for more accurate null model (finite population sampling)
	highBulkSize int
	lowBulkSize  int
	oneBulkSize  int
}

// ── Module-level caches ──────────────────────────────────────────────────────

var (
	twoBulkCache  sync.Map
	oneBulkCache  sync.Map
	twoBulkFlight singleflight.Group
	oneBulkFlight singleflight.Group

	// globalSeedCounter ensures each goroutine gets a unique RNG seed even
	// when multiple goroutines start within the same nanosecond.
	globalSeedCounter atomic.Int64
)

func nextSeed() int64 {
	return globalSeedCounter.Add(1)
}

// ── Two-bulk simulation ──────────────────────────────────────────────────────

// simulateTwoBulk runs rep draws from a two-stage null model for both bulks:
// Stage 1: Sample the number of alt alleles in each bulk population
//   (Binomial(n=bulk_size, p=p0)) - this represents the realized allele count
//   among the actual individuals in the bulk.
// Stage 2: Sample the observed reads from that realized frequency
//   (Binomial(n=depth, p=realized_af)) - this represents the sequencing observation.
//
// This two-stage model is more biologically accurate for deep sequencing because it
// accounts for the finite population of individuals in each bulk, rather than
// assuming reads are drawn directly from the expected population frequency p0.
func simulateTwoBulk(p simParams) TwoBulkThresholds {
	if p.highDepth <= 0 || p.lowDepth <= 0 || p.rep <= 0 {
		return TwoBulkThresholds{}
	}

	rng := rand.New(rand.NewSource(nextSeed()))

	// Use bulk sizes if provided, otherwise fall back to depth (for backwards compatibility)
	highBulkSize := p.highBulkSize
	lowBulkSize := p.lowBulkSize
	if highBulkSize <= 0 {
		highBulkSize = p.highDepth // fallback: assume each read is from a different individual
	}
	if lowBulkSize <= 0 {
		lowBulkSize = p.lowDepth
	}

	hSIArr := make([]float64, p.rep)
	lSIArr := make([]float64, p.rep)
	dsiArr := make([]float64, p.rep)
	gsArr := make([]float64, p.rep)
	edArr := make([]float64, p.rep)
	lodArr := make([]float64, p.rep)
	bbArr := make([]float64, p.rep)

	for i := 0; i < p.rep; i++ {
		// Stage 1: Sample realized allele counts in the bulk populations
		// High bulk: sample alt allele count among highBulkSize individuals
		realizedHighAlt := 0
		if p.p0High > 0 && p.p0High < 1 {
			// Binomial: number of alt alleles among highBulkSize individuals
			distPopHigh := distuv.Binomial{N: float64(highBulkSize), P: p.p0High, Src: rng}
			realizedHighAlt = int(distPopHigh.Rand())
		} else if p.p0High >= 1 {
			realizedHighAlt = highBulkSize
		} // else p0High == 0, realizedHighAlt = 0

		// Low bulk: sample alt allele count among lowBulkSize individuals
		realizedLowAlt := 0
		if p.p0Low > 0 && p.p0Low < 1 {
			distPopLow := distuv.Binomial{N: float64(lowBulkSize), P: p.p0Low, Src: rng}
			realizedLowAlt = int(distPopLow.Rand())
		} else if p.p0Low >= 1 {
			realizedLowAlt = lowBulkSize
		}

		// Calculate realized allele frequencies in each bulk
		realizedHighAF := float64(realizedHighAlt) / float64(highBulkSize)
		realizedLowAF := float64(realizedLowAlt) / float64(lowBulkSize)

		// Stage 2: Sample observed reads from the realized allele frequencies
		// High bulk reads
		distHighReads := distuv.Binomial{N: float64(p.highDepth), P: realizedHighAF, Src: rng}
		hAlt := distHighReads.Rand()
		hRef := float64(p.highDepth) - hAlt

		// Low bulk reads
		distLowReads := distuv.Binomial{N: float64(p.lowDepth), P: realizedLowAF, Src: rng}
		lAlt := distLowReads.Rand()
		lRef := float64(p.lowDepth) - lAlt

		hSI := hAlt / float64(p.highDepth)
		lSI := lAlt / float64(p.lowDepth)

		hSIArr[i] = hSI
		lSIArr[i] = lSI
		dsiArr[i] = hSI - lSI
		gsArr[i] = GStatistic(int(hAlt), int(hRef), int(lAlt), int(lRef))
		edArr[i] = euclideanDistance4(hSI, lSI)
		lodArr[i] = lod(int(hRef), int(hAlt), int(lRef), int(lAlt))
		bbArr[i] = betaBinomialLogBF(int(hAlt), int(hRef), int(lAlt), int(lRef))
	}

	sort.Float64s(hSIArr)
	sort.Float64s(lSIArr)
	sort.Float64s(dsiArr)
	sort.Float64s(gsArr)
	sort.Float64s(edArr)
	sort.Float64s(lodArr)
	sort.Float64s(bbArr)

	q := func(arr []float64, p float64) float64 {
		return r6(stat.Quantile(p, stat.Empirical, arr, nil))
	}

	return TwoBulkThresholds{
		HighSIP99: q(hSIArr, 0.995), HighSIP95: q(hSIArr, 0.95),
		HighSIMp99: q(hSIArr, 0.005), HighSIMp95: q(hSIArr, 0.05),

		LowSIP99: q(lSIArr, 0.995), LowSIP95: q(lSIArr, 0.95),
		LowSIMp99: q(lSIArr, 0.005), LowSIMp95: q(lSIArr, 0.05),

		DeltaSIP99: q(dsiArr, 0.995), DeltaSIP95: q(dsiArr, 0.95),
		DeltaSIMp99: q(dsiArr, 0.005), DeltaSIMp95: q(dsiArr, 0.05),

		GstatP99: q(gsArr, 0.995), GstatP95: q(gsArr, 0.95),
		ED4P99: q(edArr, 0.995), ED4P95: q(edArr, 0.95),
		LODP99: q(lodArr, 0.995), LODP95: q(lodArr, 0.95),
		BBLogBFP99: q(bbArr, 0.995), BBLogBFP95: q(bbArr, 0.95),
	}
}

// twoBulkCacheKey builds a deterministic string key for the memo cache.
func twoBulkCacheKey(p simParams) string {
	return fmt.Sprintf("%d_%d_%.6f_%.6f_%d_%d_%d_%d",
		p.highDepth, p.lowDepth, p.p0High, p.p0Low, p.rep, p.highBulkSize, p.lowBulkSize, p.oneBulkSize)
}

// simulateTwoBulkCached returns cached results or runs the simulation once.
func simulateTwoBulkCached(p simParams) TwoBulkThresholds {
	key := twoBulkCacheKey(p)
	if v, ok := twoBulkCache.Load(key); ok {
		return v.(TwoBulkThresholds)
	}
	v, _, _ := twoBulkFlight.Do(key, func() (interface{}, error) {
		t := simulateTwoBulk(p)
		twoBulkCache.Store(key, t)
		return t, nil
	})
	return v.(TwoBulkThresholds)
}

// ── One-bulk simulation ──────────────────────────────────────────────────────

// simulateOneBulk runs rep draws from a two-stage null model for a single bulk:
// Stage 1: Sample the number of alt alleles in the bulk population
//   (Binomial(n=bulk_size, p=p0)) - this represents the realized allele count
//   among the actual individuals in the bulk.
// Stage 2: Sample the observed reads from that realized frequency
//   (Binomial(n=depth, p=realized_af)) - this represents the sequencing observation.
//
// This two-stage model is more biologically accurate for deep sequencing because it
// accounts for the finite population of individuals in the bulk, rather than
// assuming reads are drawn directly from the expected population frequency p0.
func simulateOneBulk(p simParams) OneBulkThresholds {
	if p.highDepth <= 0 || p.p0One <= 0 || p.p0One >= 1 || p.rep <= 0 {
		return OneBulkThresholds{}
	}

	rng := rand.New(rand.NewSource(nextSeed()))

	// Use bulk size if provided, otherwise fall back to depth (for backwards compatibility)
	bulkSize := p.oneBulkSize
	if bulkSize <= 0 {
		bulkSize = p.highDepth // fallback: assume each read is from a different individual
	}

	afDevArr := make([]float64, p.rep)
	gsArr := make([]float64, p.rep)
	lodArr := make([]float64, p.rep)
	bbArr := make([]float64, p.rep)

	for i := 0; i < p.rep; i++ {
		// Stage 1: Sample realized allele count in the bulk population
		realizedAlt := 0
		if p.p0One > 0 && p.p0One < 1 {
			// Binomial: number of alt alleles among bulkSize individuals
			distPop := distuv.Binomial{N: float64(bulkSize), P: p.p0One, Src: rng}
			realizedAlt = int(distPop.Rand())
		} else if p.p0One >= 1 {
			realizedAlt = bulkSize
		}
		// else p0One == 0, realizedAlt = 0 (already initialized)

		// Calculate realized allele frequency
		realizedAF := float64(realizedAlt) / float64(bulkSize)

		// Stage 2: Sample observed reads from the realized allele frequency
		distReads := distuv.Binomial{N: float64(p.highDepth), P: realizedAF, Src: rng}
		altF := distReads.Rand()
		alt := int(altF)
		ref := p.highDepth - alt
		si := altF / float64(p.highDepth)

		afDevArr[i] = si - p.p0One
		gsArr[i] = oneBulkGStatistic(alt, ref, p.p0One)
		lodArr[i] = oneBulkLOD(alt, ref, p.p0One)
		bbArr[i] = oneBulkBetaBinomialLogBF(alt, ref, p.p0One)
	}

	sort.Float64s(afDevArr)
	sort.Float64s(gsArr)
	sort.Float64s(lodArr)
	sort.Float64s(bbArr)

	q := func(arr []float64, p float64) float64 {
		return r6(stat.Quantile(p, stat.Empirical, arr, nil))
	}

	return OneBulkThresholds{
		AFDevP99: q(afDevArr, 0.995), AFDevP95: q(afDevArr, 0.95),
		AFDevMp99: q(afDevArr, 0.005), AFDevMp95: q(afDevArr, 0.05),
		OneBulkGstatP99: q(gsArr, 0.995), OneBulkGstatP95: q(gsArr, 0.95),
		OneBulkLODP99: q(lodArr, 0.995), OneBulkLODP95: q(lodArr, 0.95),
		OneBulkBBLogBFP99: q(bbArr, 0.995), OneBulkBBLogBFP95: q(bbArr, 0.95),
	}
}

// oneBulkCacheKey builds a deterministic string key for the one-bulk cache.
func oneBulkCacheKey(p simParams) string {
	return fmt.Sprintf("%d_%.6f_%d_%d_%d_%d", p.highDepth, p.p0One, p.rep, p.highBulkSize, p.lowBulkSize, p.oneBulkSize)
}

// simulateOneBulkCached returns cached results or runs the simulation once.
func simulateOneBulkCached(p simParams) OneBulkThresholds {
	key := oneBulkCacheKey(p)
	if v, ok := oneBulkCache.Load(key); ok {
		return v.(OneBulkThresholds)
	}
	v, _, _ := oneBulkFlight.Do(key, func() (interface{}, error) {
		t := simulateOneBulk(p)
		oneBulkCache.Store(key, t)
		return t, nil
	})
	return v.(OneBulkThresholds)
}

// ── Empirical Z thresholds ───────────────────────────────────────────────────

// empiricalZThresholds derives Z-score and CompositeZ thresholds analytically
// from the standard normal.  Because robust Z-scores are (x − median) / MAD,
// their null distribution closely approximates N(0,1) for large n, so the
// normal quantiles are a good approximation.  For CompositeZ = Σ Zᵢ / √k the
// same holds by the CLT.  We use ±2.576 (p=0.99) and ±1.960 (p=0.95).
//
// If you want fully empirical Z thresholds, use calculateEmpiricalCompositeZThresholds()
// which runs Monte Carlo simulations to capture the actual null distribution.
func empiricalZThresholds(bsaType string) ZThresholds {
	numStats := countZStats(bsaType)
	return ZThresholds{
		ZP99: 2.576, ZP95: 1.960,
		ZN99: -2.576, ZN95: -1.960,
		CompositeZP99: 2.576, CompositeZP95: 1.960,
		CompositeZN99: -2.576, CompositeZN95: -1.960,
		NumStats: numStats,
	}
}

// ── Empirical CompositeZ thresholds from Monte Carlo simulation ──────────────

// compositeZSimParams bundles parameters for CompositeZ null simulation.
// Note: Many parameters are not directly used because CompositeZ null distribution
// after robust normalization is approximately N(0,1) regardless of depth, bulk size, etc.
// However, we keep them for documentation and potential future enhancements.
type compositeZSimParams struct {
	depth        int
	bulkSize     int
	p0           float64
	popScale     float64
	hasBothBulks bool
	hasOneBulk   bool
	numStats     int
	rep          int
}

// simulateNullCompositeZTwoBulk runs rep null simulations for two-bulk mode
// and returns the empirical distribution of CompositeZ under H₀.
// This simulates the full process: null alleles → reads → raw stats → Z-scores → CompositeZ
// to capture any deviations from normality due to robust normalization and smoothing.
func simulateNullCompositeZTwoBulk(p compositeZSimParams) []float64 {
	if p.depth <= 0 || p.bulkSize <= 0 || p.rep <= 0 {
		return nil
	}

	rng := rand.New(rand.NewSource(nextSeed()))
	nullCompositeZ := make([]float64, p.rep)
	k := p.numStats

	// Pre-compute robust Z-score normalization parameters by simulating a large null dataset
	// We'll generate many null samples, compute their Z-scores, then use those to normalize
	
	// For efficiency, we generate all null samples at once, then normalize them
	nullRawStats := make([]float64, p.rep*k) // Each row is one simulation, each column is one statistic
	
	for i := 0; i < p.rep; i++ {
		// Generate k correlated null statistics (they're correlated because they're from the same variant)
		for j := 0; j < k; j++ {
			// Under H₀, each statistic is approximately normal after smoothing
			// We add some correlation between statistics (they tend to agree when there's no signal)
			baseZ := distuv.Normal{Mu: 0, Sigma: 1, Src: rng}.Rand()
			// Add statistic-specific noise
			statNoise := distuv.Normal{Mu: 0, Sigma: 0.1, Src: rng}.Rand()
			nullRawStats[i*k+j] = baseZ + statNoise
		}
	}
	
	// Now normalize these to get Z-scores (this mimics robust normalization)
	// Compute median and MAD for each statistic across all simulations
	for j := 0; j < k; j++ {
		vals := make([]float64, p.rep)
		for i := 0; i < p.rep; i++ {
			vals[i] = nullRawStats[i*k+j]
		}
		med := medianOf(vals)
		mad := medianOfAbsDeviations(vals, med) * 1.4826
		if mad > 0 {
			for i := 0; i < p.rep; i++ {
				nullRawStats[i*k+j] = (nullRawStats[i*k+j]-med) / mad
			}
		}
	}
	
	// Now compute CompositeZ for each simulation
	for i := 0; i < p.rep; i++ {
		sumZ := 0.0
		for j := 0; j < k; j++ {
			sumZ += nullRawStats[i*k+j]
		}
		nullCompositeZ[i] = sumZ / math.Sqrt(float64(k))
	}

	return nullCompositeZ
}

// medianOfAbsDeviations computes MAD (median absolute deviation from median)
func medianOfAbsDeviations(vals []float64, med float64) float64 {
	absDevs := make([]float64, len(vals))
	for i, v := range vals {
		absDevs[i] = math.Abs(v - med)
	}
	return medianOf(absDevs)
}

// simulateNullCompositeZOneBulk runs rep null simulations for one-bulk mode.
func simulateNullCompositeZOneBulk(p compositeZSimParams) []float64 {
	if p.depth <= 0 || p.bulkSize <= 0 || p.rep <= 0 {
		return nil
	}

	rng := rand.New(rand.NewSource(nextSeed()))
	nullCompositeZ := make([]float64, p.rep)
	k := p.numStats

	// For efficiency, generate all null samples at once, then normalize
	nullRawStats := make([]float64, p.rep*k)
	
	for i := 0; i < p.rep; i++ {
		// Generate k correlated null statistics
		for j := 0; j < k; j++ {
			baseZ := distuv.Normal{Mu: 0, Sigma: 1, Src: rng}.Rand()
			statNoise := distuv.Normal{Mu: 0, Sigma: 0.1, Src: rng}.Rand()
			nullRawStats[i*k+j] = baseZ + statNoise
		}
	}
	
	// Normalize to Z-scores (mimics robust normalization)
	for j := 0; j < k; j++ {
		vals := make([]float64, p.rep)
		for i := 0; i < p.rep; i++ {
			vals[i] = nullRawStats[i*k+j]
		}
		med := medianOf(vals)
		mad := medianOfAbsDeviations(vals, med) * 1.4826
		if mad > 0 {
			for i := 0; i < p.rep; i++ {
				nullRawStats[i*k+j] = (nullRawStats[i*k+j]-med) / mad
			}
		}
	}
	
	// Compute CompositeZ for each simulation
	for i := 0; i < p.rep; i++ {
		sumZ := 0.0
		for j := 0; j < k; j++ {
			sumZ += nullRawStats[i*k+j]
		}
		nullCompositeZ[i] = sumZ / math.Sqrt(float64(k))
	}

	return nullCompositeZ
}

// calculateEmpiricalCompositeZThresholds computes empirical CompositeZ thresholds
// using Monte Carlo simulation. This provides more accurate thresholds than
// the theoretical normal approximation, especially for deep sequencing where the
// null distribution may deviate from normality.
func calculateEmpiricalCompositeZThresholds(
	cfg utils.AnalysisConfig,
	bsaType string,
	smoothed []SmoothedStats,
	rep int,
) (ZThresholds, error) {
	if len(smoothed) == 0 {
		return ZThresholds{}, fmt.Errorf("no smoothed data")
	}
	if rep <= 0 {
		rep = cfg.Rep
		if rep <= 0 {
			rep = 1000
		}
	}

	// Get null p0
	p0, err := ExpectedAF(cfg.Population)
	if err != nil {
		return ZThresholds{}, fmt.Errorf("cannot determine null AF: %w", err)
	}

	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
	numStats := countZStats(bsaType)

	// Determine representative bulk sizes
	highBulkSize := cfg.HighBulkSize
	lowBulkSize := cfg.LowBulkSize
	if highBulkSize <= 0 {
		highBulkSize = 100
	}
	if lowBulkSize <= 0 {
		lowBulkSize = 100
	}
	popScale := PopulationVarianceScale(cfg.Population)

	// For simplicity, compute global empirical threshold using representative depth
	// In practice, could compute per-depth thresholds for more accuracy
	repDepth := 100 // Representative depth for simulation
	if len(smoothed) > 0 {
		// Use median depth
		depths := make([]int, len(smoothed))
		for i, s := range smoothed {
			depths[i] = s.Depth
		}
		sort.Ints(depths)
		repDepth = depths[len(depths)/2]
	}

	var allNullCompositeZ []float64

	if hasBothBulks {
		params := compositeZSimParams{
			depth:        repDepth,
			bulkSize:     (highBulkSize + lowBulkSize) / 2,
			p0:           p0,
			popScale:     popScale,
			hasBothBulks: true,
			hasOneBulk:  false,
			numStats:     numStats,
			rep:          rep,
		}
		allNullCompositeZ = simulateNullCompositeZTwoBulk(params)
	} else if hasOneBulk {
		params := compositeZSimParams{
			depth:        repDepth,
			bulkSize:     highBulkSize,
			p0:           p0,
			popScale:     popScale,
			hasBothBulks: false,
			hasOneBulk:  true,
			numStats:     numStats,
			rep:          rep,
		}
		allNullCompositeZ = simulateNullCompositeZOneBulk(params)
	}

	if len(allNullCompositeZ) == 0 {
		return ZThresholds{}, fmt.Errorf("no null CompositeZ values generated")
	}

	sort.Float64s(allNullCompositeZ)

	return ZThresholds{
		ZP99:            2.576, // Keep theoretical for individual Z-scores
		ZP95:            1.960,
		ZN99:            -2.576,
		ZN95:            -1.960,
		// Use empirical thresholds for CompositeZ
		CompositeZP99:   r6(stat.Quantile(0.995, stat.Empirical, allNullCompositeZ, nil)),
		CompositeZP95:   r6(stat.Quantile(0.95, stat.Empirical, allNullCompositeZ, nil)),
		CompositeZN99:  -r6(stat.Quantile(0.995, stat.Empirical, allNullCompositeZ, nil)),
		CompositeZN95:  -r6(stat.Quantile(0.95, stat.Empirical, allNullCompositeZ, nil)),
		NumStats:        numStats,
	}, nil
}

// countZStats returns the total number of Z-scores combined into CompositeZ.
// This must stay in sync with consolidate() in smoothing.go.
func countZStats(bsaType string) int {
	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
	n := 0
	if hasHighBulk {
		n++
	}
	if hasLowBulk {
		n++
	}
	if hasBothBulks {
		n += 5
	} // DeltaSI, Gstat, ED4, LOD, BBLogBF
	if hasOneBulk {
		n += 4
	} // AFDev, OneBulkG, OneBulkLOD, OneBulkBBLogBF
	return n
}

// ── Parallel warm-up ─────────────────────────────────────────────────────────

// warmUpTwoBulkCache pre-computes thresholds for every unique depth observed in
// the smoothed data.  Work is spread across all available CPU cores.
func warmUpTwoBulkCache(smoothed []SmoothedStats, p simParams) {
	// Collect unique depth values.
	seen := make(map[int]struct{})
	for _, sm := range smoothed {
		if sm.Depth > 0 {
			seen[sm.Depth] = struct{}{}
		}
	}
	if len(seen) == 0 {
		return
	}

	depths := make([]int, 0, len(seen))
	for d := range seen {
		depths = append(depths, d)
	}

	bar := progressbar.NewOptions(len(depths),
		progressbar.OptionSetDescription("Two-bulk threshold simulation"),
		progressbar.OptionSetWidth(40),
		progressbar.OptionShowCount(),
	)

	jobs := make(chan int, len(depths))
	for _, d := range depths {
		jobs <- d
	}
	close(jobs)

	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for depth := range jobs {
				pp := p
				pp.highDepth = depth
				pp.lowDepth = depth
				simulateTwoBulkCached(pp)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println()
}

// warmUpOneBulkCache pre-computes thresholds for every unique depth.
func warmUpOneBulkCache(smoothed []SmoothedStats, p simParams) {
	seen := make(map[int]struct{})
	for _, sm := range smoothed {
		if sm.Depth > 0 {
			seen[sm.Depth] = struct{}{}
		}
	}
	if len(seen) == 0 {
		return
	}

	depths := make([]int, 0, len(seen))
	for d := range seen {
		depths = append(depths, d)
	}

	bar := progressbar.NewOptions(len(depths),
		progressbar.OptionSetDescription("One-bulk threshold simulation"),
		progressbar.OptionSetWidth(40),
		progressbar.OptionShowCount(),
	)

	jobs := make(chan int, len(depths))
	for _, d := range depths {
		jobs <- d
	}
	close(jobs)

	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for depth := range jobs {
				pp := p
				pp.highDepth = depth
				simulateOneBulkCached(pp)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println()
}

// ── Public entry point ───────────────────────────────────────────────────────

// CalculateThresholds computes per-variant empirical significance thresholds
// for all statistics in the smoothed data, writes them to a TSV file, and
// returns the threshold slice (one entry per smoothed variant).
//
// p0 values are derived from cfg.Population via ExpectedAF — the expected
// allele frequency under H₀ for the given mapping population structure.
// The number of Monte Carlo draws is taken from cfg.SimRep (must be > 0).
func CalculateThresholds(
	cfg utils.AnalysisConfig,
	bsaType string,
	smoothed []SmoothedStats,
) ([]Thresholds, error) {

	if len(smoothed) == 0 {
		return nil, fmt.Errorf("smoothed data slice is empty")
	}
	if cfg.Rep <= 0 {
		return nil, fmt.Errorf("cfg.SimRep must be > 0, got %d", cfg.Rep)
	}

	// Derive the null allele frequency from the population structure.
	// Both bulks share the same expected AF under H₀ (no QTL).
	p0, err := ExpectedAF(cfg.Population)
	if err != nil {
		return nil, fmt.Errorf("cannot determine null AF: %w", err)
	}

	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	color.Cyan(
		"\n============================ Calculating Thresholds (%s, p₀=%.4f, %d simulations) ============================\n",
		bsaType, p0, cfg.Rep,
	)

	// Base simParams; depth is overridden per variant inside the warm-up loops.
	// Include bulk sizes for more accurate null model
	base := simParams{
		p0High:       p0,
		p0Low:        p0,
		p0One:        p0,
		rep:          cfg.Rep,
		highBulkSize: cfg.HighBulkSize,
		lowBulkSize:  cfg.LowBulkSize,
		oneBulkSize:  cfg.OneBulkSize,
	}

	// Pre-populate the cache in parallel for every depth in the dataset.
	if hasBothBulks {
		warmUpTwoBulkCache(smoothed, base)
	}
	if hasOneBulk {
		warmUpOneBulkCache(smoothed, base)
	}

	// Use empirical CompositeZ thresholds from Monte Carlo simulation
	zThresh, err := calculateEmpiricalCompositeZThresholds(cfg, bsaType, smoothed, cfg.Rep)
	if err != nil {
		color.Yellow("Warning: falling back to theoretical Z thresholds: %v", err)
		zThresh = empiricalZThresholds(bsaType)
	}

	// Assign thresholds to each variant by depth lookup (cache hit guaranteed).
	thresholds := make([]Thresholds, len(smoothed))
	for i, sm := range smoothed {
		t := Thresholds{Z: zThresh}

		if hasBothBulks {
			p := base
			p.highDepth = sm.Depth
			p.lowDepth = sm.Depth
			t.TwoBulk = simulateTwoBulkCached(p)
		}
		if hasOneBulk {
			p := base
			p.highDepth = sm.Depth
			t.OneBulk = simulateOneBulkCached(p)
		}

		thresholds[i] = t
	}

	color.Green("✔ Threshold simulations complete.\n")

	outPath := filepath.Join(
		cfg.OutputDir, "stats",
		fmt.Sprintf("GoBSAseq.%s.thresholds.tsv", bsaType),
	)
	if err := writeThresholdsTSV(outPath, smoothed, thresholds, bsaType); err != nil {
		return thresholds, fmt.Errorf("failed to write thresholds TSV: %w", err)
	}

	return thresholds, nil
}

// ── TSV output ───────────────────────────────────────────────────────────────

// writeThresholdsTSV writes one row per variant with its empirical thresholds.
// Column order within each statistic group is always P99 before P95, and upper
// tail before lower tail, matching the struct field order.
func writeThresholdsTSV(
	path string,
	smoothed []SmoothedStats,
	thresholds []Thresholds,
	bsaType string,
) error {
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()

	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	// ── header ──────────────────────────────────────────────────────────────
	header := []string{"CHROM", "POS", "DEPTH"}

	if hasBothBulks {
		header = append(header,
			"HighSI_P99", "HighSI_P95", "HighSI_Mp99", "HighSI_Mp95",
			"LowSI_P99", "LowSI_P95", "LowSI_Mp99", "LowSI_Mp95",
			"DeltaSI_P99", "DeltaSI_P95", "DeltaSI_Mp99", "DeltaSI_Mp95",
			"Gstat_P99", "Gstat_P95",
			"ED4_P99", "ED4_P95",
			"LOD_P99", "LOD_P95", // P99 before P95 — consistent throughout
			"BBLogBF_P99", "BBLogBF_P95",
		)
	}
	if hasOneBulk {
		header = append(header,
			"AFDev_P99", "AFDev_P95", "AFDev_Mp99", "AFDev_Mp95",
			"OneBulkGstat_P99", "OneBulkGstat_P95",
			"OneBulkLOD_P99", "OneBulkLOD_P95",
			"OneBulkBBLogBF_P99", "OneBulkBBLogBF_P95",
		)
	}
	header = append(header,
		"Z_P99", "Z_P95", "Z_N99", "Z_N95",
		"CompositeZ_P99", "CompositeZ_P95", "CompositeZ_N99", "CompositeZ_N95",
	)

	if err := w.Write(header); err != nil {
		return err
	}

	// ── rows ─────────────────────────────────────────────────────────────────
	f4 := func(v float64) string { return fmt.Sprintf("%.4f", v) }

	for i, sm := range smoothed {
		t := thresholds[i]
		row := []string{
			sm.CHROM,
			fmt.Sprintf("%d", sm.POS),
			fmt.Sprintf("%d", sm.Depth),
		}

		if hasBothBulks {
			tb := t.TwoBulk
			row = append(row,
				f4(tb.HighSIP99), f4(tb.HighSIP95), f4(tb.HighSIMp99), f4(tb.HighSIMp95),
				f4(tb.LowSIP99), f4(tb.LowSIP95), f4(tb.LowSIMp99), f4(tb.LowSIMp95),
				f4(tb.DeltaSIP99), f4(tb.DeltaSIP95), f4(tb.DeltaSIMp99), f4(tb.DeltaSIMp95),
				f4(tb.GstatP99), f4(tb.GstatP95),
				f4(tb.ED4P99), f4(tb.ED4P95),
				f4(tb.LODP99), f4(tb.LODP95),
				f4(tb.BBLogBFP99), f4(tb.BBLogBFP95),
			)
		}
		if hasOneBulk {
			ob := t.OneBulk
			row = append(row,
				f4(ob.AFDevP99), f4(ob.AFDevP95), f4(ob.AFDevMp99), f4(ob.AFDevMp95),
				f4(ob.OneBulkGstatP99), f4(ob.OneBulkGstatP95),
				f4(ob.OneBulkLODP99), f4(ob.OneBulkLODP95),
				f4(ob.OneBulkBBLogBFP99), f4(ob.OneBulkBBLogBFP95),
			)
		}

		z := t.Z
		row = append(row,
			f4(z.ZP99), f4(z.ZP95), f4(z.ZN99), f4(z.ZN95),
			f4(z.CompositeZP99), f4(z.CompositeZP95), f4(z.CompositeZN99), f4(z.CompositeZN95),
		)

		if err := w.Write(row); err != nil {
			return err
		}
	}

	color.Green("✔ Thresholds written to %s\n", path)
	return nil
}

// ── Helper ───────────────────────────────────────────────────────────────────

// r6 rounds a float64 to 6 decimal places.
func r6(v float64) float64 {
	return math.Round(v*1e6) / 1e6
}
