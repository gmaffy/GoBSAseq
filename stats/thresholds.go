package stats

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
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

type Thresholds struct {
	TwoBulk TwoBulkThresholds
	OneBulk OneBulkThresholds
	Z       ZThresholds
	HMM     HMMThresholds
}

// HMMThresholds are genome-wide thresholds for the HMM LOD: the 99% and 95%
// quantiles of the maximum LOD over simulated no-QTL genomes (hmmnull.go).
// Zero when the HMM LOD is not computed.
type HMMThresholds struct {
	LODP99, LODP95 float64
}

type TwoBulkThresholds struct {
	HighSIP99, HighSIP95     float64
	HighSIMp99, HighSIMp95   float64
	LowSIP99, LowSIP95       float64
	LowSIMp99, LowSIMp95     float64
	DeltaSIP99, DeltaSIP95   float64
	DeltaSIMp99, DeltaSIMp95 float64

	GstatP99, GstatP95 float64
	ED4P99, ED4P95     float64
}

type OneBulkThresholds struct {
	AFDevP99, AFDevP95               float64
	AFDevMp99, AFDevMp95             float64 // lower tail
	OneBulkGstatP99, OneBulkGstatP95 float64
}

type ZThresholds struct {
	ZP99, ZP95                   float64
	CompositeZP99, CompositeZP95 float64
	NumStats                     int
}

// ── Simulation parameters ────────────────────────────────────────────────────

// simParams bundles the inputs to a single null simulation run.
type simParams struct {
	highDepth int
	lowDepth  int     // low-bulk depth (two-bulk); simulated separately from highDepth
	p0High    float64 // expected AF in high bulk under H₀
	p0Low     float64 // expected AF in low bulk under H₀ (usually == p0High)
	p0One     float64 // expected AF for one-bulk mode
	rep       int
	// Bulk sizes for more accurate null model (finite population sampling)
	highBulkSize int
	lowBulkSize  int
	oneBulkSize  int
	population   string
	// absDeltaSI: polarity unknown (2b, --ref-parent none), so the ΔSI
	// thresholds are for |ΔSI| (upper tail only; the lower-tail ones are NaN).
	absDeltaSI bool
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

// samplePoolAF samples the allele frequency carried by a pool of individuals.
// It models diploid genotype dosage, rather than treating each individual as a
// single allele draw. F2/F3/F4 are successive selfing generations from an F1;
// RIL is the fully fixed limit. BCnH/L are recurrent-parent backcrosses.
func samplePoolAF(rng *rand.Rand, population string, bulkSize int, p0 float64) float64 {
	if bulkSize <= 0 {
		return p0
	}
	p := strings.ToUpper(strings.TrimSpace(population))
	selfGenerations := map[string]int{"F2": 1, "F3": 2, "F4": 3}
	if g, ok := selfGenerations[p]; ok || p == "RIL" {
		het := 0.0
		if p != "RIL" {
			het = math.Pow(0.5, float64(g))
		}
		pHom := (1 - het) / 2
		copies := 0
		for range bulkSize {
			u := rng.Float64()
			switch {
			case u < pHom:
				copies += 0
			case u < pHom+het:
				copies++
			default:
				copies += 2
			}
		}
		return float64(copies) / float64(2*bulkSize)
	}
	if strings.HasPrefix(p, "BC") {
		suffix := p[2:]
		isHigh := strings.HasSuffix(suffix, "H")
		suffix = strings.TrimSuffix(strings.TrimSuffix(suffix, "H"), "L")
		isF2 := strings.HasSuffix(suffix, "F2")
		suffix = strings.TrimSuffix(suffix, "F2")
		n, err := strconv.Atoi(suffix)
		if err == nil && n >= 1 {
			// BCn: recurrent homozygote 1-q, heterozygote q, with q = 0.5^n.
			// BCnF2 (one selfing of BCn): heterozygotes segregate 1:2:1, so
			// recurrent hom = 1-q+q/4, het = q/2, donor hom = q/4.
			q := math.Pow(0.5, float64(n))
			pHet, pDonor := q, 0.0
			if isF2 {
				pHet, pDonor = q/2, q/4
			}
			recurrent := 0 // recurrent-parent allele copies in the pool
			for range bulkSize {
				u := rng.Float64()
				switch {
				case u < pDonor:
					// donor homozygote: 0 recurrent copies
				case u < pDonor+pHet:
					recurrent++
				default:
					recurrent += 2
				}
			}
			f := float64(recurrent) / float64(2*bulkSize)
			// SI is the frequency of the high-parent allele: the recurrent
			// allele for BC*H, the donor allele for BC*L.
			if !isHigh {
				f = 1 - f
			}
			return f
		}
	}
	// Unknown structures are rejected by ExpectedAF before this point. Keep a
	// conservative diploid fallback for callers that use this helper directly.
	copies := int(distuv.Binomial{N: float64(2 * bulkSize), P: p0, Src: rng}.Rand())
	return float64(copies) / float64(2*bulkSize)
}

// ── Two-bulk simulation ──────────────────────────────────────────────────────

// simulateTwoBulk runs rep draws from a two-stage null model for both bulks:
// Stage 1: Sample genotype dosage in each bulk population according to the
// selected cross. This represents the realised allele frequency among the
// actual diploid individuals in the bulk.
//
// Stage 2: Sample the observed reads from that realized frequency
//
//	(Binomial(n=depth, p=realized_af)) - this represents the sequencing observation.
//
// This two-stage model is more biologically accurate for deep sequencing because it
// accounts for the finite population of individuals in each bulk, rather than
// assuming reads are drawn directly from the expected population frequency p0.
func simulateTwoBulk(p simParams) TwoBulkThresholds {
	// Bulk sizes are validated by CalculateThresholds (cfg.BulkSizes).
	if p.highDepth <= 0 || p.lowDepth <= 0 || p.rep <= 0 || p.highBulkSize <= 0 || p.lowBulkSize <= 0 {
		return TwoBulkThresholds{}
	}

	rng := rand.New(rand.NewSource(nextSeed()))

	highBulkSize := p.highBulkSize
	lowBulkSize := p.lowBulkSize

	hSIArr := make([]float64, p.rep)
	lSIArr := make([]float64, p.rep)
	dsiArr := make([]float64, p.rep)
	gsArr := make([]float64, p.rep)
	edArr := make([]float64, p.rep)

	for i := 0; i < p.rep; i++ {
		// Stage 1: realised diploid genotype dosage in each pool.
		realizedHighAF := samplePoolAF(rng, p.population, highBulkSize, p.p0High)
		realizedLowAF := samplePoolAF(rng, p.population, lowBulkSize, p.p0Low)

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
		if p.absDeltaSI {
			dsiArr[i] = math.Abs(dsiArr[i])
		}
		gsArr[i] = GStatistic(int(hAlt), int(hRef), int(lAlt), int(lRef))
		edArr[i] = euclideanDistance4(hSI, lSI)
	}

	sort.Float64s(hSIArr)
	sort.Float64s(lSIArr)
	sort.Float64s(dsiArr)
	sort.Float64s(gsArr)
	sort.Float64s(edArr)

	q := func(arr []float64, p float64) float64 {
		return r6(stat.Quantile(p, stat.Empirical, arr, nil))
	}

	// Signed ΔSI: P99/Mp99 are the 0.5% tails (99% two-sided). |ΔSI| folds both
	// tails, so its 99th percentile gives the same 1% false-positive rate, and
	// there is no lower tail.
	dsiP99, dsiP95 := q(dsiArr, 0.995), q(dsiArr, 0.95)
	dsiMp99, dsiMp95 := q(dsiArr, 0.005), q(dsiArr, 0.05)
	if p.absDeltaSI {
		dsiP99, dsiP95 = q(dsiArr, 0.99), q(dsiArr, 0.95)
		dsiMp99, dsiMp95 = math.NaN(), math.NaN()
	}

	return TwoBulkThresholds{
		HighSIP99: q(hSIArr, 0.995), HighSIP95: q(hSIArr, 0.95),
		HighSIMp99: q(hSIArr, 0.005), HighSIMp95: q(hSIArr, 0.05),

		LowSIP99: q(lSIArr, 0.995), LowSIP95: q(lSIArr, 0.95),
		LowSIMp99: q(lSIArr, 0.005), LowSIMp95: q(lSIArr, 0.05),

		DeltaSIP99: dsiP99, DeltaSIP95: dsiP95,
		DeltaSIMp99: dsiMp99, DeltaSIMp95: dsiMp95,

		GstatP99: q(gsArr, 0.995), GstatP95: q(gsArr, 0.95),
		ED4P99: q(edArr, 0.995), ED4P95: q(edArr, 0.95),
	}
}

// twoBulkCacheKey builds a deterministic string key for the memo cache.
func twoBulkCacheKey(p simParams) string {
	return fmt.Sprintf("%d_%d_%.6f_%.6f_%d_%d_%d_%d_%s_%t",
		p.highDepth, p.lowDepth, p.p0High, p.p0Low, p.rep, p.highBulkSize, p.lowBulkSize, p.oneBulkSize,
		strings.ToUpper(strings.TrimSpace(p.population)), p.absDeltaSI)
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
// Stage 1: Sample genotype dosage in the bulk population according to the
// selected cross. This represents the realised allele frequency among the
// actual diploid individuals in the bulk.
//
// Stage 2: Sample the observed reads from that realized frequency
//
//	(Binomial(n=depth, p=realized_af)) - this represents the sequencing observation.
//
// This two-stage model is more biologically accurate for deep sequencing because it
// accounts for the finite population of individuals in the bulk, rather than
// assuming reads are drawn directly from the expected population frequency p0.
func simulateOneBulk(p simParams) OneBulkThresholds {
	// Bulk size is validated by CalculateThresholds (cfg.BulkSizes).
	if p.highDepth <= 0 || p.p0One <= 0 || p.p0One >= 1 || p.rep <= 0 || p.oneBulkSize <= 0 {
		return OneBulkThresholds{}
	}

	rng := rand.New(rand.NewSource(nextSeed()))

	bulkSize := p.oneBulkSize

	afDevArr := make([]float64, p.rep)
	gsArr := make([]float64, p.rep)

	for i := 0; i < p.rep; i++ {
		// Stage 1: realised diploid genotype dosage in the pool.
		realizedAF := samplePoolAF(rng, p.population, bulkSize, p.p0One)

		// Stage 2: Sample observed reads from the realized allele frequency
		distReads := distuv.Binomial{N: float64(p.highDepth), P: realizedAF, Src: rng}
		altF := distReads.Rand()
		alt := int(altF)
		ref := p.highDepth - alt
		si := altF / float64(p.highDepth)

		afDevArr[i] = si - p.p0One
		gsArr[i] = oneBulkGStatistic(alt, ref, p.p0One)
	}

	sort.Float64s(afDevArr)
	sort.Float64s(gsArr)

	q := func(arr []float64, p float64) float64 {
		return r6(stat.Quantile(p, stat.Empirical, arr, nil))
	}

	return OneBulkThresholds{
		AFDevP99: q(afDevArr, 0.995), AFDevP95: q(afDevArr, 0.95),
		AFDevMp99: q(afDevArr, 0.005), AFDevMp95: q(afDevArr, 0.05),
		OneBulkGstatP99: q(gsArr, 0.995), OneBulkGstatP95: q(gsArr, 0.95),
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

// ── CompositeZ ───────────────────────────────────────────────────────────────

// fallbackZThresholds is used only if the null simulation cannot run. It uses
// standard-normal quantiles (2.576 / 1.960).
func fallbackZThresholds(k int) ZThresholds {
	return ZThresholds{ZP99: 2.576, ZP95: 1.960, CompositeZP99: 2.576, CompositeZP95: 1.960, NumStats: k}
}

// countZStats returns k, the number of Z-scores consolidate combines into
// CompositeZ: two per bulk layout (|ZΔSI| + ZG, or |ZAFDev| + ZG1) plus ZHMMLOD
// when the HMM LOD is computed. TestCompositeZStatCountMatchesConsolidate
// checks it against consolidate.
func countZStats(bsaType string, hmm bool) int {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
	n := 0
	if hasBothBulks || hasOneBulk {
		n = 2
		if hmm {
			n++
		}
	}
	return n
}

// ── Parallel warm-up ─────────────────────────────────────────────────────────

// siteBulkDepths is the (high, low) depth pair a two-bulk site is simulated
// at: its own per-bulk depths, or the effective depth for both when those are
// absent.
func siteBulkDepths(sm SmoothedStats) (int, int) {
	if sm.HighDepth > 0 && sm.LowDepth > 0 {
		return sm.HighDepth, sm.LowDepth
	}
	return sm.Depth, sm.Depth
}

// warmUpTwoBulkCache pre-computes thresholds for every unique (high, low)
// depth pair observed in the smoothed data, so unequal bulk depths are
// simulated as they are. Work is spread across all available CPU cores.
func warmUpTwoBulkCache(smoothed []SmoothedStats, p simParams) {
	// Collect unique depth pairs.
	seen := make(map[[2]int]struct{})
	for _, sm := range smoothed {
		if h, l := siteBulkDepths(sm); h > 0 && l > 0 {
			seen[[2]int{h, l}] = struct{}{}
		}
	}
	if len(seen) == 0 {
		return
	}

	depths := make([][2]int, 0, len(seen))
	for d := range seen {
		depths = append(depths, d)
	}

	bar := progressbar.NewOptions(len(depths),
		progressbar.OptionSetDescription("Two-bulk threshold simulation"),
		progressbar.OptionSetWidth(40),
		progressbar.OptionShowCount(),
	)

	jobs := make(chan [2]int, len(depths))
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
				pp.highDepth, pp.lowDepth = depth[0], depth[1]
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
// CalculateThresholds returns per-site thresholds. Per-statistic thresholds
// (SI, ΔSI, G, ED⁴, AFDev, G1) come from single-SNP simulations at each site's
// depths; HMM LOD and CompositeZ thresholds come from the genome-wide null
// simulation over raw's markers (hmmnull.go).
func CalculateThresholds(
	cfg utils.AnalysisConfig,
	bsaType string,
	raw []BSAstats,
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

	highBulkSize, lowBulkSize, oneBulkSize, err := cfg.BulkSizes(bsaType)
	if err != nil {
		return nil, err
	}

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
		highBulkSize: highBulkSize,
		lowBulkSize:  lowBulkSize,
		oneBulkSize:  oneBulkSize,
		population:   cfg.Population,
		absDeltaSI:   !cfg.SignedDeltaSI(bsaType),
	}

	// Pre-populate the cache in parallel for every depth in the dataset.
	if hasBothBulks {
		warmUpTwoBulkCache(smoothed, base)
	}
	if hasOneBulk {
		warmUpOneBulkCache(smoothed, base)
	}

	// HMM LOD and CompositeZ thresholds from the genome-wide null.
	zThresh, hmmThresh, err := genomeNullThresholds(cfg, bsaType, raw)
	if err != nil {
		color.Yellow("Warning: genome-wide null failed (%v); CompositeZ falls back to normal quantiles and the HMM LOD has no threshold", err)
		zThresh = fallbackZThresholds(countZStats(bsaType, HMMAvailable(cfg, bsaType)))
		hmmThresh = HMMThresholds{LODP99: math.Inf(1), LODP95: math.Inf(1)}
	}

	// Assign thresholds to each variant by depth lookup (cache hit guaranteed).
	thresholds := make([]Thresholds, len(smoothed))
	for i, sm := range smoothed {
		t := Thresholds{Z: zThresh, HMM: hmmThresh}

		if hasBothBulks {
			p := base
			p.highDepth, p.lowDepth = siteBulkDepths(sm)
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
			"HighDepth", "LowDepth", // the depth pair the two-bulk null was simulated at
			"HighSI_P99", "HighSI_P95", "HighSI_Mp99", "HighSI_Mp95",
			"LowSI_P99", "LowSI_P95", "LowSI_Mp99", "LowSI_Mp95",
			"DeltaSI_P99", "DeltaSI_P95", "DeltaSI_Mp99", "DeltaSI_Mp95",
			"Gstat_P99", "Gstat_P95",
			"ED4_P99", "ED4_P95",
		)
	}
	if hasOneBulk {
		header = append(header,
			"AFDev_P99", "AFDev_P95", "AFDev_Mp99", "AFDev_Mp95",
			"OneBulkGstat_P99", "OneBulkGstat_P95",
		)
	}
	header = append(header,
		"HMMLOD_P99", "HMMLOD_P95", // genome-wide (max over null genomes)
		"Z_P99", "Z_P95",
		"CompositeZ_P99", "CompositeZ_P95",
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
			hd, ld := siteBulkDepths(sm)
			row = append(row,
				fmt.Sprintf("%d", hd), fmt.Sprintf("%d", ld),
				f4(tb.HighSIP99), f4(tb.HighSIP95), f4(tb.HighSIMp99), f4(tb.HighSIMp95),
				f4(tb.LowSIP99), f4(tb.LowSIP95), f4(tb.LowSIMp99), f4(tb.LowSIMp95),
				f4(tb.DeltaSIP99), f4(tb.DeltaSIP95), f4(tb.DeltaSIMp99), f4(tb.DeltaSIMp95),
				f4(tb.GstatP99), f4(tb.GstatP95),
				f4(tb.ED4P99), f4(tb.ED4P95),
			)
		}
		if hasOneBulk {
			ob := t.OneBulk
			row = append(row,
				f4(ob.AFDevP99), f4(ob.AFDevP95), f4(ob.AFDevMp99), f4(ob.AFDevMp95),
				f4(ob.OneBulkGstatP99), f4(ob.OneBulkGstatP95),
			)
		}

		z := t.Z
		row = append(row,
			f4(t.HMM.LODP99), f4(t.HMM.LODP95),
			f4(z.ZP99), f4(z.ZP95),
			f4(z.CompositeZP99), f4(z.CompositeZP95),
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
