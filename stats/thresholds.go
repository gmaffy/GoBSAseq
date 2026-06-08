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
	"sync"
	"time"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// Thresholds bundles all three families for one SmoothedStats window.
type Thresholds struct {
	TwoBulk TwoBulkThresholds
	OneBulk OneBulkThresholds
	Z       ZThresholds
}

type TwoBulkThresholds struct {
	HighSIP99, HighSIP95, HighSIMp99, HighSIMp95     float64
	LowSIP99, LowSIP95, LowSIMp99, LowSIMp95         float64
	DeltaSIP99, DeltaSIP95, DeltaSIMp99, DeltaSIMp95 float64
	GstatP99, GstatP95                               float64
	ED4P99, ED4P95                                   float64
	LODP99, LODP95                                   float64
	BBLogBFP99, BBLogBFP95                           float64
}

type OneBulkThresholds struct {
	AFDevP99, AFDevP95, AFDevMp99, AFDevMp95 float64
	OneBulkGstatP99, OneBulkGstatP95         float64
	OneBulkLODP99, OneBulkLODP95             float64
	OneBulkBBLogBFP99, OneBulkBBLogBFP95     float64
}

type ZThresholds struct {
	ZP99, ZP95, ZN99, ZN95                                     float64
	CompositeZP99, CompositeZP95, CompositeZN99, CompositeZN95 float64
	NumStats                                                   int
}

var (
	twoBulkCache sync.Map
	oneBulkCache sync.Map
)

// WriteThresholdsToTSV handles standard disk I/O serialization for reporting metrics.
// WriteThresholdsToTSV handles standard disk I/O serialization for reporting metrics.
// It maps out all simulation quantiles across both single and multi-bulk analysis frameworks.
func WriteThresholdsToTSV(outputPath string, smoothed []SmoothedStats, thresholds []Thresholds, bsaType string) error {
	if err := os.MkdirAll(filepath.Dir(outputPath), 0755); err != nil {
		return err
	}

	file, err := os.Create(outputPath)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	writer.Comma = '\t'
	defer writer.Flush()

	_, _, hasBothBulks, hasOneBulk := bulkFlags(bsaType)

	// Build out headers with strict structural parity to our underlying structs
	header := []string{"CHROM", "POS", "DEPTH"}

	if hasBothBulks {
		header = append(header,
			"HighSI_P99", "HighSI_P95", "HighSI_Mp99", "HighSI_Mp95",
			"LowSI_P99", "LowSI_P95", "LowSI_Mp99", "LowSI_Mp95",
			"DeltaSI_P99", "DeltaSI_P95", "DeltaSI_Mp99", "DeltaSI_Mp95",
			"Gstat_P99", "Gstat_P95",
			"ED4_P99", "ED4_P95",
			"LOD_P95", "LOD_P99",
			"BBLogBF_P95", "BBLogBF_P99",
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

	// Analytical normal limits included for thorough verification checks
	header = append(header, "Z_P99", "Z_P95", "CompositeZ_P99", "CompositeZ_P95")

	if err := writer.Write(header); err != nil {
		return err
	}

	// Format utility to keep floating point outputs clean and structured
	f4 := func(v float64) string { return fmt.Sprintf("%.4f", v) }

	for i, sm := range smoothed {
		t := thresholds[i]
		row := []string{
			sm.CHROM,
			fmt.Sprintf("%d", sm.POS),
			fmt.Sprintf("%d", sm.Depth),
		}

		if hasBothBulks {
			row = append(row,
				f4(t.TwoBulk.HighSIP99), f4(t.TwoBulk.HighSIP95), f4(t.TwoBulk.HighSIMp99), f4(t.TwoBulk.HighSIMp95),
				f4(t.TwoBulk.LowSIP99), f4(t.TwoBulk.LowSIP95), f4(t.TwoBulk.LowSIMp99), f4(t.TwoBulk.LowSIMp95),
				f4(t.TwoBulk.DeltaSIP99), f4(t.TwoBulk.DeltaSIP95), f4(t.TwoBulk.DeltaSIMp99), f4(t.TwoBulk.DeltaSIMp95),
				f4(t.TwoBulk.GstatP99), f4(t.TwoBulk.GstatP95),
				f4(t.TwoBulk.ED4P99), f4(t.TwoBulk.ED4P95),
				f4(t.TwoBulk.LODP99), f4(t.TwoBulk.LODP95),
				f4(t.TwoBulk.BBLogBFP99), f4(t.TwoBulk.BBLogBFP95),
			)
		}

		if hasOneBulk {
			row = append(row,
				f4(t.OneBulk.AFDevP99), f4(t.OneBulk.AFDevP95), f4(t.OneBulk.AFDevMp99), f4(t.OneBulk.AFDevMp95),
				f4(t.OneBulk.OneBulkGstatP99), f4(t.OneBulk.OneBulkGstatP95),
				f4(t.OneBulk.OneBulkLODP99), f4(t.OneBulk.OneBulkLODP95),
				f4(t.OneBulk.OneBulkBBLogBFP99), f4(t.OneBulk.OneBulkBBLogBFP95),
			)
		}

		// Push analytical reference criteria out to the final structural line elements
		row = append(row,
			f4(t.Z.ZP99), f4(t.Z.ZP95),
			f4(t.Z.CompositeZP99), f4(t.Z.CompositeZP95),
		)

		if err := writer.Write(row); err != nil {
			return err
		}
	}

	color.Green("✔ Complete empirical significance limits logs saved to: %s", outputPath)
	return nil
}

// ---------------------------------------------------------------------------
// Internal Calculations & Caching Infrastructure (Unchanged math)
// ---------------------------------------------------------------------------

func calcZThresholds(numStats int) ZThresholds {
	const (
		zSig  = 3.0
		zSugg = 2.0
	)
	if numStats < 1 {
		numStats = 1
	}
	return ZThresholds{
		ZP99: zSig, ZP95: zSugg, ZN99: -zSig, ZN95: -zSugg,
		CompositeZP99: zSig, CompositeZP95: zSugg, CompositeZN99: -zSig, CompositeZN95: -zSugg,
		NumStats: numStats,
	}
}

func numZStats(bsaType string) int {
	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := bulkFlags(bsaType)
	n := 0
	if hasHighBulk {
		n++
	}
	if hasLowBulk {
		n++
	}
	if hasBothBulks {
		n += 5
	}
	if hasOneBulk {
		n += 4
	}
	return n
}

func calcTwoBulkThresholds(highBulkDP, lowBulkDP int, highSmAF, lowSmAF float64, rep int) TwoBulkThresholds {
	if highBulkDP <= 0 || lowBulkDP <= 0 || rep <= 0 {
		return TwoBulkThresholds{}
	}
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	distHigh := distuv.Binomial{N: float64(highBulkDP), P: highSmAF, Src: rng}
	distLow := distuv.Binomial{N: float64(lowBulkDP), P: lowSmAF, Src: rng}

	highSIArr, lowSIArr, dsiArr := make([]float64, rep), make([]float64, rep), make([]float64, rep)
	gsArr, edArr, lodArr, bbArr := make([]float64, rep), make([]float64, rep), make([]float64, rep), make([]float64, rep)

	for i := 0; i < rep; i++ {
		hAlt := distHigh.Rand()
		lAlt := distLow.Rand()
		hRef, lRef := float64(highBulkDP)-hAlt, float64(lowBulkDP)-lAlt
		hSI := math.Round((hAlt/float64(highBulkDP))*1e6) / 1e6
		lSI := math.Round((lAlt/float64(lowBulkDP))*1e6) / 1e6

		highSIArr[i], lowSIArr[i] = hSI, lSI
		dsiArr[i] = math.Round((hSI-lSI)*1e6) / 1e6
		gsArr[i] = math.Round(GStatistic(int(hAlt), int(hRef), int(lAlt), int(lRef))*1e6) / 1e6
		edArr[i] = math.Round(euclideanDistance4(hSI, lSI)*1e6) / 1e6
		lodArr[i] = math.Round(lod(int(hRef), int(hAlt), int(lRef), int(lAlt))*1e6) / 1e6
		bbArr[i] = math.Round(betaBinomialLogBF(int(hAlt), int(hRef), int(lAlt), int(lRef))*1e6) / 1e6
	}

	sort.Float64s(highSIArr)
	sort.Float64s(lowSIArr)
	sort.Float64s(dsiArr)
	sort.Float64s(gsArr)
	sort.Float64s(edArr)
	sort.Float64s(lodArr)
	sort.Float64s(bbArr)

	r6 := func(v float64) float64 { return math.Round(v*1e6) / 1e6 }
	q := func(arr []float64, p float64) float64 { return r6(stat.Quantile(p, stat.Empirical, arr, nil)) }

	return TwoBulkThresholds{
		HighSIP99: q(highSIArr, 0.995), HighSIP95: q(highSIArr, 0.95), HighSIMp99: q(highSIArr, 0.005), HighSIMp95: q(highSIArr, 0.05),
		LowSIP99: q(lowSIArr, 0.995), LowSIP95: q(lowSIArr, 0.95), LowSIMp99: q(lowSIArr, 0.005), LowSIMp95: q(lowSIArr, 0.05),
		DeltaSIP99: q(dsiArr, 0.995), DeltaSIP95: q(dsiArr, 0.95), DeltaSIMp99: q(dsiArr, 0.005), DeltaSIMp95: q(dsiArr, 0.05),
		GstatP99: q(gsArr, 0.995), GstatP95: q(gsArr, 0.95),
		ED4P99: q(edArr, 0.995), ED4P95: q(edArr, 0.95),
		LODP99: q(lodArr, 0.995), LODP95: q(lodArr, 0.95),
		BBLogBFP99: q(bbArr, 0.995), BBLogBFP95: q(bbArr, 0.95),
	}
}

func calcTwoBulkCached(highBulkDP, lowBulkDP int, highSmAF, lowSmAF float64, rep int) TwoBulkThresholds {
	key := fmt.Sprintf("%d_%d_%.6f_%.6f_%d", highBulkDP, lowBulkDP, highSmAF, lowSmAF, rep)
	if v, ok := twoBulkCache.Load(key); ok {
		return v.(TwoBulkThresholds)
	}
	t := calcTwoBulkThresholds(highBulkDP, lowBulkDP, highSmAF, lowSmAF, rep)
	actual, _ := twoBulkCache.LoadOrStore(key, t)
	return actual.(TwoBulkThresholds)
}

func calcAllTwoBulkThresholds(smoothed []SmoothedStats, highSmAF, lowSmAF float64, rep int) {
	type depthKey struct{ depth int }
	seen := make(map[depthKey]bool)
	for _, sm := range smoothed {
		if sm.Depth > 0 {
			seen[depthKey{sm.Depth}] = true
		}
	}
	keys := make([]depthKey, 0, len(seen))
	for k := range seen {
		keys = append(keys, k)
	}
	if len(keys) == 0 {
		return
	}

	bar := progressbar.NewOptions(len(keys), progressbar.OptionSetDescription("Two-bulk permutation processing"), progressbar.OptionSetWidth(30), progressbar.OptionShowCount())
	keyChan := make(chan depthKey, len(keys))
	for _, k := range keys {
		keyChan <- k
	}
	close(keyChan)

	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for k := range keyChan {
				calcTwoBulkCached(k.depth, k.depth, highSmAF, lowSmAF, rep)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
}

func calcOneBulkThresholds(depth int, p0 float64, rep int) OneBulkThresholds {
	if depth <= 0 || p0 <= 0 || p0 >= 1 || rep <= 0 {
		return OneBulkThresholds{}
	}
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	dist := distuv.Binomial{N: float64(depth), P: p0, Src: rng}

	afDevArr, gsArr, lodArr, bbArr := make([]float64, rep), make([]float64, rep), make([]float64, rep), make([]float64, rep)
	for i := 0; i < rep; i++ {
		succF := dist.Rand()
		succ, fail := int(succF), depth-int(succF)
		si := math.Round((succF/float64(depth))*1e6) / 1e6
		afDevArr[i] = math.Round((si-p0)*1e6) / 1e6
		gsArr[i] = math.Round(oneBulkGStatistic(succ, fail, p0)*1e6) / 1e6
		lodArr[i] = math.Round(oneBulkLOD(succ, fail, p0)*1e6) / 1e6
		bbArr[i] = math.Round(oneBulkBetaBinomialLogBF(succ, fail, p0)*1e6) / 1e6
	}

	sort.Float64s(afDevArr)
	sort.Float64s(gsArr)
	sort.Float64s(lodArr)
	sort.Float64s(bbArr)
	r6 := func(v float64) float64 { return math.Round(v*1e6) / 1e6 }
	q := func(arr []float64, p float64) float64 { return r6(stat.Quantile(p, stat.Empirical, arr, nil)) }

	return OneBulkThresholds{
		AFDevP99: q(afDevArr, 0.995), AFDevP95: q(afDevArr, 0.95), AFDevMp99: q(afDevArr, 0.005), AFDevMp95: q(afDevArr, 0.05),
		OneBulkGstatP99: q(gsArr, 0.995), OneBulkGstatP95: q(gsArr, 0.95),
		OneBulkLODP99: q(lodArr, 0.995), OneBulkLODP95: q(lodArr, 0.95),
		OneBulkBBLogBFP99: q(bbArr, 0.995), OneBulkBBLogBFP95: q(bbArr, 0.95),
	}
}

func calcOneBulkCached(depth int, p0 float64, rep int) OneBulkThresholds {
	key := fmt.Sprintf("%d_%.6f_%d", depth, p0, rep)
	if v, ok := oneBulkCache.Load(key); ok {
		return v.(OneBulkThresholds)
	}
	t := calcOneBulkThresholds(depth, p0, rep)
	actual, _ := oneBulkCache.LoadOrStore(key, t)
	return actual.(OneBulkThresholds)
}

func calcAllOneBulkThresholds(smoothed []SmoothedStats, p0 float64, rep int) {
	seen := make(map[int]bool)
	for _, sm := range smoothed {
		if sm.Depth > 0 {
			seen[sm.Depth] = true
		}
	}
	depths := make([]int, 0, len(seen))
	for d := range seen {
		depths = append(depths, d)
	}
	if len(depths) == 0 {
		return
	}

	bar := progressbar.NewOptions(len(depths), progressbar.OptionSetDescription("One-bulk permutation processing"), progressbar.OptionSetWidth(30), progressbar.OptionShowCount())
	depthChan := make(chan int, len(depths))
	for _, d := range depths {
		depthChan <- d
	}
	close(depthChan)

	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for d := range depthChan {
				calcOneBulkCached(d, p0, rep)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
}

func CalculateThresholds(cfg utils.AnalysisConfig, bsaType string, smoothed []SmoothedStats, highSmAF, lowSmAF float64, p0 float64, rep int) ([]Thresholds, error) {
	if len(smoothed) == 0 {
		return nil, fmt.Errorf("smoothed data slice is empty")
	}
	if rep <= 0 {
		return nil, fmt.Errorf("invalid permutation replication count: %d", rep)
	}

	_, _, hasBothBulks, hasOneBulk := bulkFlags(bsaType)

	color.Cyan("\n============================ Calculating Thresholds (%d simulations) ============================\n", rep)

	if hasBothBulks {
		calcAllTwoBulkThresholds(smoothed, highSmAF, lowSmAF, rep)
	}
	if hasOneBulk && p0 > 0 && p0 < 1 {
		calcAllOneBulkThresholds(smoothed, p0, rep)
	}

	numStats := numZStats(bsaType)
	thresholdResults := make([]Thresholds, len(smoothed))

	for i := range smoothed {
		sm := smoothed[i]
		t := Thresholds{}

		if hasBothBulks {
			t.TwoBulk = calcTwoBulkCached(sm.Depth, sm.Depth, highSmAF, lowSmAF, rep)
		}
		if hasOneBulk && p0 > 0 && p0 < 1 {
			t.OneBulk = calcOneBulkCached(sm.Depth, p0, rep)
		}
		t.Z = calcZThresholds(numStats)

		thresholdResults[i] = t
	}

	color.Green("✔ Threshold simulations complete.")

	// Export generated arrays directly out to a TSV file
	tsvPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.thresholds.tsv"), bsaType)
	if err := WriteThresholdsToTSV(tsvPath, smoothed, thresholdResults, bsaType); err != nil {
		return thresholdResults, fmt.Errorf("failed to save TSV logs: %w", err)
	}

	return thresholdResults, nil
}
