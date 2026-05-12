// Package twobulk contains plotting and QTL-detection for GoBSAseq two-bulk analysis.
//
// GenerateHtmlPlotsAndQTL now produces THREE separate HTML files:
//
//	GoBSAseq_IndividualPlots.html   – one raw-value line chart per statistic per chromosome
//	GoBSAseq_NormalizedOverlay.html – per-chromosome threshold-relative overlay
//	GoBSAseq_RobustZScore.html      – genome-wide robust Z-score view, one chart per chromosome
//
// Normalization change (v2):
//
//	The old approach divided each value by its per-chromosome average p99 threshold.
//	That confounds biological signal with sequencing depth: chromosomes with lower
//	coverage have a higher p99 threshold, so identical peaks shrink on low-coverage
//	chromosomes.
//
//	The new approach uses a genome-wide robust Z-score:
//	    z = (x − median_genome) / (MAD_genome × 1.4826)
//	where MAD is the median absolute deviation and 1.4826 makes it consistent with σ
//	under normality.  The top 1 % of values per statistic are excluded when estimating
//	the background median/MAD so that genuine QTL peaks do not inflate the spread.
//	Reference lines are drawn at z = ±2 (suggestive) and z = ±3 (significant).
package twobulk

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// ---------------------------------------------------------------------------
// Domain types
// ---------------------------------------------------------------------------

// BSAstats holds the raw statistics for a single SNP position.
type BSAstats struct {
	CHROM      string
	POS        int64
	REF        string
	ALT        string
	HighParGT  []int
	LowParGT   []int
	HighBulkGT []int
	HighBulkAD string
	LowBulkGT  []int
	LowBulkAD  string

	HighBulkL int
	HighBulkH int
	LowBulkL  int
	LowBulkH  int
	HighSI    float64
	LowSI     float64

	DeltaSI float64
	Gstat   float64
	ED      float64
	LOD     float64
	BBLogBF float64

	Depth int
}

// SmoothedStats holds the window-averaged statistics for one genomic window.
type SmoothedStats struct {
	CHROM          string
	POS            int64
	DeltaSI        float64
	Gstat          float64
	ED             float64
	LOD            float64
	BBLogBF        float64
	HighSI         float64
	LowSI          float64
	NumSNPs        int
	MeanHighBulkDP int
	MeanLowBulkDP  int
}

// Thresholds holds the permutation-derived significance levels for each statistic.
type Thresholds struct {
	DsiP99  float64
	DsiP95  float64
	DsiMp99 float64
	DsiMp95 float64

	GsP99 float64
	GsP95 float64

	EdP99 float64
	EdP95 float64

	LodP99 float64
	LodP95 float64

	BbP99 float64
	BbP95 float64

	HighP99  float64
	HighP95  float64
	HighMp99 float64
	HighMp95 float64

	LowP99  float64
	LowP95  float64
	LowMp99 float64
	LowMp95 float64
}

// QTLRecord holds the detected QTL interval and its peak value.
type QTLRecord struct {
	Chrom string
	Start int64
	Stop  int64
	Peak  float64
	Stat  string
	CI    string
}

// BRMBlock holds one BRM-style allele-frequency-difference block interval.
type BRMBlock struct {
	Chrom     string
	Start     int64
	Stop      int64
	PeakPos   int64
	Peak      float64
	Threshold float64
}

// ---------------------------------------------------------------------------
// Caching
// ---------------------------------------------------------------------------

var thresholdCache sync.Map

// ---------------------------------------------------------------------------
// Analysis pipeline
// ---------------------------------------------------------------------------

// RunTwoBulkTwoParents is the main entry point for the two-bulk two-parent analysis.
func RunTwoBulkTwoParents(cfg utils.AnalysisConfig) {
	highParIdx := cfg.HighParentIdx
	highParDP := cfg.HighParentDepth
	lowParIdx := cfg.LowParentIdx
	lowParDP := cfg.LowParentDepth
	highBulkIdx := cfg.HighBulkIdx
	highBulkDP := cfg.HighBulkDepth
	lowBulkIdx := cfg.LowBulkIdx
	lowBulkDP := cfg.LowBulkDepth
	vcfRdr := cfg.Rdr
	outDir := cfg.OutputDir

	windowSize := int64(cfg.WindowSize)
	stepSize := int64(cfg.StepSize)
	rep := cfg.Rep
	pop := cfg.Population

	highSmAF := utils.SimulateAF(pop, float64(cfg.HighBulkSize), cfg.Rep)
	lowSmAF := utils.SimulateAF(pop, float64(cfg.LowBulkSize), cfg.Rep)

	overallStart := time.Now()
	variantChan := make(chan *vcfgo.Variant, 10000)
	statsChan := make(chan BSAstats, 10000)
	rawWriteChan := make(chan BSAstats, 10000)

	numWorkers := runtime.NumCPU()
	var workerWG sync.WaitGroup

	var rawWG sync.WaitGroup
	rawWG.Add(1)
	go func() {
		defer rawWG.Done()
		if err := writeRawTSV(filepath.Join(outDir, "GoBSAseq.raw.tsv"), rawWriteChan); err != nil {
			color.Red("Error writing raw TSV: %v", err)
		}
	}()

	color.Cyan("============================ Calculating Statistics =============================\n\n")
	bar := progressbar.Default(-1, "Processing variants")

	for i := 0; i < numWorkers; i++ {
		workerWG.Add(1)
		go func() {
			defer workerWG.Done()
			for variant := range variantChan {
				if !GoodVariants(variant, highParIdx, highParDP, lowParIdx, lowParDP, highBulkIdx, highBulkDP, lowBulkIdx, lowBulkDP) {
					_ = bar.Add(1)
					continue
				}

				lpGT := variant.Samples[lowParIdx].GT
				hbGT := variant.Samples[highBulkIdx].GT
				lbGT := variant.Samples[lowBulkIdx].GT
				hbDP, lbDP := variant.Samples[highBulkIdx].DP, variant.Samples[lowBulkIdx].DP

				if hbDP == 0 || lbDP == 0 {
					_ = bar.Add(1)
					continue
				}

				hbRefDep, _ := variant.Samples[highBulkIdx].RefDepth()
				hbAltDeps, _ := variant.Samples[highBulkIdx].AltDepths()
				lbRefDep, _ := variant.Samples[lowBulkIdx].RefDepth()
				lbAltDeps, _ := variant.Samples[lowBulkIdx].AltDepths()

				if len(hbAltDeps) == 0 || len(lbAltDeps) == 0 {
					_ = bar.Add(1)
					continue
				}

				var hbL, hbH, lbL, lbH int
				if hbGT[0] == lpGT[0] {
					hbL, hbH = hbRefDep, hbAltDeps[0]
				} else {
					hbL, hbH = hbAltDeps[0], hbRefDep
				}
				if lbGT[0] == lpGT[0] {
					lbL, lbH = lbRefDep, lbAltDeps[0]
				} else {
					lbL, lbH = lbAltDeps[0], lbRefDep
				}

				hSI := float64(hbH) / float64(hbDP)
				lSI := float64(lbH) / float64(lbDP)
				minDepth := hbDP
				if lbDP < minDepth {
					minDepth = lbDP
				}

				s := BSAstats{
					CHROM:      variant.Chromosome,
					POS:        int64(variant.Pos),
					REF:        variant.Reference,
					ALT:        variant.Alt()[0],
					HighParGT:  variant.Samples[highParIdx].GT,
					LowParGT:   variant.Samples[lowParIdx].GT,
					HighBulkGT: variant.Samples[highBulkIdx].GT,
					HighBulkAD: fmt.Sprintf("%d,%d", hbRefDep, hbAltDeps[0]),
					LowBulkGT:  variant.Samples[lowBulkIdx].GT,
					LowBulkAD:  fmt.Sprintf("%d,%d", lbRefDep, lbAltDeps[0]),
					HighBulkL:  hbL,
					HighBulkH:  hbH,
					LowBulkL:   lbL,
					LowBulkH:   lbH,
					HighSI:     math.Round(hSI*1e6) / 1e6,
					LowSI:      math.Round(lSI*1e6) / 1e6,
					DeltaSI:    math.Round((hSI-lSI)*1e6) / 1e6,
					Gstat:      math.Round(GStatistic(hbH, hbL, lbH, lbL)*1e6) / 1e6,
					ED:         math.Round(math.Abs(hSI-lSI)*1e6) / 1e6,
					LOD:        math.Round(lod(hbL, hbH, lbL, lbH)*1e6) / 1e6,
					BBLogBF:    math.Round(betaBinomialLogBF(hbH, hbL, lbH, lbL)*1e6) / 1e6,
					Depth:      minDepth,
				}
				statsChan <- s
				rawWriteChan <- s
				_ = bar.Add(1)
			}
		}()
	}

	go func() {
		for v := vcfRdr.Read(); v != nil; v = vcfRdr.Read() {
			variantChan <- v
		}
		close(variantChan)
	}()

	go func() {
		workerWG.Wait()
		close(statsChan)
		close(rawWriteChan)
	}()

	chromStats := make(map[string][]BSAstats)
	for s := range statsChan {
		chromStats[s.CHROM] = append(chromStats[s.CHROM], s)
	}
	_ = bar.Finish()
	rawWG.Wait()

	color.Cyan("\n============================ Smoothing Statistics =============================\n\n")

	var allSmoothed []SmoothedStats
	for chrom, stats := range chromStats {
		color.Yellow("Smoothing %s: %d SNPs", chrom, len(stats))
		smoothed := smoothChromosome(stats, windowSize, stepSize)
		allSmoothed = append(allSmoothed, smoothed...)
	}

	color.Cyan("\n============================ Calculating Thresholds (%d simulations per depth pair) ==============================\n\n", rep)
	calcAllThresholds(allSmoothed, highSmAF, lowSmAF, rep)
	color.Green("\nThreshold calculations complete.")

	color.Cyan("\n=========================================== Writing Smoothed TSV =================================================\n\n")
	smoothTSV := filepath.Join(outDir, "GoBSAseq.smooth.tsv")
	if err := writeSmoothedTSV(smoothTSV, allSmoothed, highSmAF, lowSmAF, rep); err != nil {
		color.Red("Error writing smoothed TSV: %v", err)
	} else {
		color.Green("Wrote %d smoothed windows to %s", len(allSmoothed), smoothTSV)
	}
	color.Green("Raw stats written to %s", filepath.Join(outDir, "GoBSAseq.raw.tsv"))
	color.Green("\nTotal time: %s\n", time.Since(overallStart).Round(time.Second))

	color.Cyan("\n============================ Generating HTML Plots & QTLs ========================================\n\n")
	htmlFile := filepath.Join(outDir, "GoBSAseq_InteractivePlots.html")
	qtlFile := filepath.Join(outDir, "GoBSAseq_QTL.tsv")

	if err := GenerateHtmlPlotsAndQTL(
		allSmoothed,
		highSmAF, lowSmAF,
		cfg.HighBulkSize, cfg.LowBulkSize,
		pop, cfg.Alphas, rep,
		htmlFile, qtlFile,
	); err != nil {
		color.Red("Error generating Plots and QTLs: %v", err)
	} else {
		color.Green("Interactive HTML plots written under %s", outDir)
		color.Green("QTL tabular results written to %s", qtlFile)
	}
}

// ---------------------------------------------------------------------------
// Threshold simulation and caching
// ---------------------------------------------------------------------------

func calcThresholds(highBulkDP, lowBulkDP int, highSmAF, lowSmAF float64, rep int) Thresholds {
	if highBulkDP <= 0 || lowBulkDP <= 0 || rep <= 0 {
		return Thresholds{}
	}

	src := rand.NewSource(time.Now().UnixNano())
	rng := rand.New(src)
	distHigh := distuv.Binomial{N: float64(highBulkDP), P: highSmAF, Src: rng}
	distLow := distuv.Binomial{N: float64(lowBulkDP), P: lowSmAF, Src: rng}

	highSIArr := make([]float64, rep)
	lowSIArr := make([]float64, rep)
	dsiArr := make([]float64, rep)
	gsArr := make([]float64, rep)
	edArr := make([]float64, rep)
	lodArr := make([]float64, rep)
	bbArr := make([]float64, rep)

	for i := 0; i < rep; i++ {
		hAlt := distHigh.Rand()
		lAlt := distLow.Rand()
		hRef := float64(highBulkDP) - hAlt
		lRef := float64(lowBulkDP) - lAlt

		hSI := math.Round((hAlt/float64(highBulkDP))*1e6) / 1e6
		lSI := math.Round((lAlt/float64(lowBulkDP))*1e6) / 1e6

		highSIArr[i] = hSI
		lowSIArr[i] = lSI
		dsiArr[i] = math.Round((hSI-lSI)*1e6) / 1e6
		gsArr[i] = math.Round(GStatistic(int(hAlt), int(hRef), int(lAlt), int(lRef))*1e6) / 1e6
		edArr[i] = math.Round(math.Abs(hSI-lSI)*1e6) / 1e6
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

	return Thresholds{
		HighP99:  math.Round(stat.Quantile(0.995, stat.Empirical, highSIArr, nil)*1e6) / 1e6,
		HighP95:  math.Round(stat.Quantile(0.95, stat.Empirical, highSIArr, nil)*1e6) / 1e6,
		HighMp99: math.Round(stat.Quantile(0.005, stat.Empirical, highSIArr, nil)*1e6) / 1e6,
		HighMp95: math.Round(stat.Quantile(0.05, stat.Empirical, highSIArr, nil)*1e6) / 1e6,

		LowP99:  math.Round(stat.Quantile(0.995, stat.Empirical, lowSIArr, nil)*1e6) / 1e6,
		LowP95:  math.Round(stat.Quantile(0.95, stat.Empirical, lowSIArr, nil)*1e6) / 1e6,
		LowMp99: math.Round(stat.Quantile(0.005, stat.Empirical, lowSIArr, nil)*1e6) / 1e6,
		LowMp95: math.Round(stat.Quantile(0.05, stat.Empirical, lowSIArr, nil)*1e6) / 1e6,

		DsiP99:  math.Round(stat.Quantile(0.995, stat.Empirical, dsiArr, nil)*1e6) / 1e6,
		DsiP95:  math.Round(stat.Quantile(0.95, stat.Empirical, dsiArr, nil)*1e6) / 1e6,
		DsiMp99: math.Round(stat.Quantile(0.005, stat.Empirical, dsiArr, nil)*1e6) / 1e6,
		DsiMp95: math.Round(stat.Quantile(0.05, stat.Empirical, dsiArr, nil)*1e6) / 1e6,

		GsP99: math.Round(stat.Quantile(0.995, stat.Empirical, gsArr, nil)*1e6) / 1e6,
		GsP95: math.Round(stat.Quantile(0.95, stat.Empirical, gsArr, nil)*1e6) / 1e6,

		EdP99: math.Round(stat.Quantile(0.995, stat.Empirical, edArr, nil)*1e6) / 1e6,
		EdP95: math.Round(stat.Quantile(0.95, stat.Empirical, edArr, nil)*1e6) / 1e6,

		LodP99: math.Round(stat.Quantile(0.995, stat.Empirical, lodArr, nil)*1e6) / 1e6,
		LodP95: math.Round(stat.Quantile(0.95, stat.Empirical, lodArr, nil)*1e6) / 1e6,

		BbP99: math.Round(stat.Quantile(0.995, stat.Empirical, bbArr, nil)*1e6) / 1e6,
		BbP95: math.Round(stat.Quantile(0.95, stat.Empirical, bbArr, nil)*1e6) / 1e6,
	}
}

func calcThresholdsCached(highBulkDP, lowBulkDP int, highSmAF, lowSmAF float64, rep int) Thresholds {
	key := fmt.Sprintf("%d_%d_%.6f_%.6f_%d", highBulkDP, lowBulkDP, highSmAF, lowSmAF, rep)
	if v, ok := thresholdCache.Load(key); ok {
		return v.(Thresholds)
	}
	t := calcThresholds(highBulkDP, lowBulkDP, highSmAF, lowSmAF, rep)
	actual, _ := thresholdCache.LoadOrStore(key, t)
	return actual.(Thresholds)
}

func calcAllThresholds(allSmoothed []SmoothedStats, highSmAF, lowSmAF float64, rep int) {
	type depthPair struct{ h, l int }
	seen := make(map[depthPair]bool)
	for _, sm := range allSmoothed {
		if sm.MeanHighBulkDP > 0 && sm.MeanLowBulkDP > 0 {
			seen[depthPair{sm.MeanHighBulkDP, sm.MeanLowBulkDP}] = true
		}
	}

	pairs := make([]depthPair, 0, len(seen))
	for dp := range seen {
		pairs = append(pairs, dp)
	}
	if len(pairs) == 0 {
		return
	}

	bar := progressbar.NewOptions(len(pairs),
		progressbar.OptionSetDescription("Computing thresholds"),
		progressbar.OptionSetWidth(40),
		progressbar.OptionShowCount(),
		progressbar.OptionSetTheme(progressbar.Theme{
			Saucer:        "=",
			SaucerHead:    ">",
			SaucerPadding: " ",
			BarStart:      "[",
			BarEnd:        "]",
		}),
	)

	pairChan := make(chan depthPair, len(pairs))
	for _, dp := range pairs {
		pairChan <- dp
	}
	close(pairChan)

	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for dp := range pairChan {
				calcThresholdsCached(dp.h, dp.l, highSmAF, lowSmAF, rep)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println()
}

// ---------------------------------------------------------------------------
// Statistic helpers
// ---------------------------------------------------------------------------

func isHomozygous(gt []int) bool {
	if len(gt) == 0 {
		return false
	}
	first := gt[0]
	if first < 0 {
		return false
	}
	for _, a := range gt[1:] {
		if a != first {
			return false
		}
	}
	return true
}

func GStatistic(highBulkAlt, highBulkRef, lowBulkAlt, lowBulkRef int) float64 {
	hba := float64(highBulkAlt) + 0.5
	hbr := float64(highBulkRef) + 0.5
	lba := float64(lowBulkAlt) + 0.5
	lbr := float64(lowBulkRef) + 0.5

	highBulkTotal := hba + hbr
	lowBulkTotal := lba + lbr
	total := highBulkTotal + lowBulkTotal
	if total == 0 {
		return 0
	}

	expHighAlt := highBulkTotal * (hba + lba) / total
	expHighRef := highBulkTotal * (hbr + lbr) / total
	expLowAlt := lowBulkTotal * (hba + lba) / total
	expLowRef := lowBulkTotal * (hbr + lbr) / total

	g := 0.0
	if hba > 0 {
		g += hba * math.Log(hba/expHighAlt)
	}
	if hbr > 0 {
		g += hbr * math.Log(hbr/expHighRef)
	}
	if lba > 0 {
		g += lba * math.Log(lba/expLowAlt)
	}
	if lbr > 0 {
		g += lbr * math.Log(lbr/expLowRef)
	}
	return 2 * g
}

func logBeta(a, b float64) float64 {
	la, _ := math.Lgamma(a)
	lb, _ := math.Lgamma(b)
	lab, _ := math.Lgamma(a + b)
	return la + lb - lab
}

func lod(ref1, alt1, ref2, alt2 int) float64 {
	n1 := float64(ref1 + alt1)
	n2 := float64(ref2 + alt2)
	total := n1 + n2
	if n1 == 0 || n2 == 0 || total == 0 {
		return 0
	}

	const eps = 1e-10
	clamp := func(p float64) float64 {
		if p <= 0 {
			return eps
		}
		if p >= 1 {
			return 1 - eps
		}
		return p
	}

	p1 := clamp(float64(alt1) / n1)
	p2 := clamp(float64(alt2) / n2)
	p0 := clamp(float64(alt1+alt2) / total)

	logLik := func(k, n, p float64) float64 {
		if n == 0 {
			return 0
		}
		return k*math.Log(p) + (n-k)*math.Log(1-p)
	}

	logL1 := logLik(float64(alt1), n1, p1) + logLik(float64(alt2), n2, p2)
	logL0 := logLik(float64(alt1), n1, p0) + logLik(float64(alt2), n2, p0)
	return (logL1 - logL0) / math.Log(10)
}

func betaBinomialLogBF(highSucc, highFail, lowSucc, lowFail int) float64 {
	alphaH, betaH := 0.5, 0.5
	alphaL, betaL := 0.5, 0.5

	logAlt := logBeta(alphaH+float64(highSucc), betaH+float64(highFail)) - logBeta(alphaH, betaH)
	logAlt += logBeta(alphaL+float64(lowSucc), betaL+float64(lowFail)) - logBeta(alphaL, betaL)

	alpha0, beta0 := 0.5, 0.5
	logNull := logBeta(alpha0+float64(highSucc+lowSucc), beta0+float64(highFail+lowFail)) - logBeta(alpha0, beta0)
	return (logAlt - logNull) / math.Log(10)
}

// GoodVariants filters to biallelic, fully called, homozygous divergent parents
// with sufficient parent and bulk depth.
func GoodVariants(v *vcfgo.Variant, highPar, highParDP, lowPar, lowParDP, highBulk, highBulkDP, lowBulk, lowBulkDP int) bool {
	indices := []int{highPar, lowPar, highBulk, lowBulk}
	if len(v.Alt()) != 1 {
		return false
	}

	for _, idx := range indices {
		if idx < 0 || idx >= len(v.Samples) {
			return false
		}
		s := v.Samples[idx]
		if len(s.GT) == 0 {
			return false
		}
		for _, allele := range s.GT {
			if allele < 0 {
				return false
			}
		}
	}

	hpGT, lpGT := v.Samples[highPar].GT, v.Samples[lowPar].GT
	hpDP, lpDP := v.Samples[highPar].DP, v.Samples[lowPar].DP
	hbDP, lbDP := v.Samples[highBulk].DP, v.Samples[lowBulk].DP

	if !isHomozygous(hpGT) || !isHomozygous(lpGT) || hpGT[0] == lpGT[0] {
		return false
	}
	return hpDP >= highParDP && lpDP >= lowParDP && hbDP >= highBulkDP && lbDP >= lowBulkDP
}

// ---------------------------------------------------------------------------
// Smoothing and TSV output
// ---------------------------------------------------------------------------

func tricubeWeight(d, D float64) float64 {
	if D <= 0 || d >= D {
		return 0
	}
	x := 1 - math.Pow(d/D, 3)
	return x * x * x
}

func smoothChromosome(stats []BSAstats, windowSize int64, step int64) []SmoothedStats {
	if len(stats) == 0 || windowSize <= 0 || step <= 0 {
		return nil
	}

	sort.Slice(stats, func(i, j int) bool { return stats[i].POS < stats[j].POS })

	var smoothed []SmoothedStats
	chrom := stats[0].CHROM
	minPos := stats[0].POS
	maxPos := stats[len(stats)-1].POS

	for center := minPos; center <= maxPos; center += step {
		windowStart := center - windowSize/2
		windowEnd := center + windowSize/2

		var (
			sumDeltaSI, sumED                float64
			sumWeightDeltaSI, sumWeightED    float64
			sumGstat, sumLOD, sumBBLogBF     float64
			sumHighSI, sumLowSI, sumWeightSI float64
			sumHighDP, sumLowDP              float64
			nSNPs                            int
		)

		for _, s := range stats {
			if s.POS < windowStart || s.POS > windowEnd {
				continue
			}
			nSNPs++

			d := math.Abs(float64(s.POS - center))
			w := tricubeWeight(d, float64(windowSize)/2)
			depthWeight := math.Sqrt(float64(s.Depth))
			wStat := w * depthWeight

			sumDeltaSI += s.DeltaSI * wStat
			sumWeightDeltaSI += wStat
			sumED += s.ED * wStat
			sumWeightED += wStat
			sumGstat += s.Gstat * wStat
			sumLOD += s.LOD * wStat
			sumBBLogBF += s.BBLogBF * wStat
			sumHighSI += s.HighSI * wStat
			sumLowSI += s.LowSI * wStat
			sumWeightSI += wStat
			sumHighDP += float64(s.HighBulkL + s.HighBulkH)
			sumLowDP += float64(s.LowBulkL + s.LowBulkH)
		}

		if nSNPs == 0 {
			continue
		}

		sm := SmoothedStats{
			CHROM:          chrom,
			POS:            center,
			NumSNPs:        nSNPs,
			MeanHighBulkDP: int(sumHighDP / float64(nSNPs)),
			MeanLowBulkDP:  int(sumLowDP / float64(nSNPs)),
		}
		if sumWeightDeltaSI > 0 {
			sm.DeltaSI = sumDeltaSI / sumWeightDeltaSI
			sm.Gstat = sumGstat / sumWeightDeltaSI
			sm.LOD = sumLOD / sumWeightDeltaSI
			sm.BBLogBF = sumBBLogBF / sumWeightDeltaSI
		}
		if sumWeightED > 0 {
			sm.ED = sumED / sumWeightED
		}
		if sumWeightSI > 0 {
			sm.HighSI = sumHighSI / sumWeightSI
			sm.LowSI = sumLowSI / sumWeightSI
		}

		smoothed = append(smoothed, sm)
	}

	return smoothed
}

func writeRawTSV(filename string, statsChan <-chan BSAstats) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	fmt.Fprintln(w, "CHROM\tPOS\tREF\tALT\tHighParGT\tLowParGT\tHighBulkGT\tHighBulkAD\tLowBulkGT\tLowBulkAD\tHighBulkL\tHighBulkH\tLowBulkL\tLowBulkH\tHighSI\tLowSI\tDeltaSI\tGstat\tED\tLOD\tBBLogBF\tDepth")
	for s := range statsChan {
		fmt.Fprintf(w, "%s\t%d\t%s\t%s\t%v\t%v\t%v\t%s\t%v\t%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
			s.CHROM, s.POS, s.REF, s.ALT,
			s.HighParGT, s.LowParGT, s.HighBulkGT, s.HighBulkAD, s.LowBulkGT, s.LowBulkAD,
			s.HighBulkL, s.HighBulkH, s.LowBulkL, s.LowBulkH,
			s.HighSI, s.LowSI, s.DeltaSI, s.Gstat, s.ED, s.LOD, s.BBLogBF, s.Depth)
	}
	return nil
}

func writeSmoothedTSV(filename string, data []SmoothedStats, highSmAF, lowSmAF float64, rep int) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	header := "CHROM\tPOS\tHighSI\tLowSI\tDeltaSI\tGstat\tED\tLOD\tBBLogBF\tNumSNPs\tMeanHighDP\tMeanLowDP" +
		"\tHighSI_p99\tHighSI_p95\tHighSI_m_p99\tHighSI_m_p95" +
		"\tLowSI_p99\tLowSI_p95\tLowSI_m_p99\tLowSI_m_p95" +
		"\tDeltaSI_p99\tDeltaSI_p95\tDeltaSI_m_p99\tDeltaSI_m_p95" +
		"\tGstat_p99\tGstat_p95" +
		"\tED_p99\tED_p95" +
		"\tLOD_p99\tLOD_p95" +
		"\tBBLogBF_p99\tBBLogBF_p95"
	fmt.Fprintln(w, header)

	for _, d := range data {
		t := calcThresholdsCached(d.MeanHighBulkDP, d.MeanLowBulkDP, highSmAF, lowSmAF, rep)
		row := fmt.Sprintf(
			"%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d"+
				"\t%.6f\t%.6f\t%.6f\t%.6f"+
				"\t%.6f\t%.6f\t%.6f\t%.6f"+
				"\t%.6f\t%.6f\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f",
			d.CHROM, d.POS,
			d.HighSI, d.LowSI, d.DeltaSI, d.Gstat, d.ED, d.LOD, d.BBLogBF,
			d.NumSNPs, d.MeanHighBulkDP, d.MeanLowBulkDP,
			t.HighP99, t.HighP95, t.HighMp99, t.HighMp95,
			t.LowP99, t.LowP95, t.LowMp99, t.LowMp95,
			t.DsiP99, t.DsiP95, t.DsiMp99, t.DsiMp95,
			t.GsP99, t.GsP95,
			t.EdP99, t.EdP95,
			t.LodP99, t.LodP95,
			t.BbP99, t.BbP95,
		)
		fmt.Fprintln(w, row)
	}
	return nil
}

// ---------------------------------------------------------------------------
// Robust Z-score helpers
// ---------------------------------------------------------------------------

// robustBackground computes the median and MAD of vals, excluding the top
// trimFrac proportion of values so that genuine QTL peaks do not inflate the
// spread estimate.  trimFrac = 0.01 excludes the top 1 %.
func robustBackground(vals []float64, trimFrac float64) (median, mad float64) {
	if len(vals) == 0 {
		return 0, 0
	}
	sorted := make([]float64, len(vals))
	copy(sorted, vals)
	sort.Float64s(sorted)

	// trim the top trimFrac of values
	cutIdx := int(math.Round(float64(len(sorted)) * (1.0 - trimFrac)))
	if cutIdx < 1 {
		cutIdx = 1
	}
	trimmed := sorted[:cutIdx]

	median = quantile(trimmed, 0.5)

	// MAD = median(|x − median|)
	devs := make([]float64, len(trimmed))
	for i, v := range trimmed {
		devs[i] = math.Abs(v - median)
	}
	sort.Float64s(devs)
	mad = quantile(devs, 0.5)
	return median, mad
}

// quantile returns the p-th quantile of a pre-sorted slice.
func quantile(sorted []float64, p float64) float64 {
	n := len(sorted)
	if n == 0 {
		return 0
	}
	if p <= 0 {
		return sorted[0]
	}
	if p >= 1 {
		return sorted[n-1]
	}
	pos := p * float64(n-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	if lo == hi {
		return sorted[lo]
	}
	frac := pos - float64(lo)
	return sorted[lo]*(1-frac) + sorted[hi]*frac
}

// robustZScore normalises vals to genome-wide robust Z-scores.
// The 1.4826 factor makes the MAD-based scale consistent with σ under normality.
func robustZScore(vals []float64, median, mad float64) []float64 {
	out := make([]float64, len(vals))
	scale := mad * 1.4826
	if scale == 0 {
		// Flat signal — return zeros rather than NaN.
		return out
	}
	for i, v := range vals {
		out[i] = (v - median) / scale
	}
	return out
}

// ---------------------------------------------------------------------------
// Genome-wide normalisation state
// ---------------------------------------------------------------------------

// genomeWideNorms collects per-statistic genome-wide robust Z parameters.
type genomeWideNorms struct {
	hiMed, hiMAD   float64
	liMed, liMAD   float64
	dsiMed, dsiMAD float64
	gsMed, gsMAD   float64
	edMed, edMAD   float64
	lodMed, lodMAD float64
	bblMed, bblMAD float64
}

// computeGenomeWideNorms performs a first pass over all smoothed data to
// derive robust background statistics for every statistic.
func computeGenomeWideNorms(allSmoothed []SmoothedStats) genomeWideNorms {
	hi := make([]float64, 0, len(allSmoothed))
	li := make([]float64, 0, len(allSmoothed))
	dsi := make([]float64, 0, len(allSmoothed))
	gs := make([]float64, 0, len(allSmoothed))
	ed := make([]float64, 0, len(allSmoothed))
	lod := make([]float64, 0, len(allSmoothed))
	bbl := make([]float64, 0, len(allSmoothed))

	for _, s := range allSmoothed {
		hi = append(hi, s.HighSI)
		li = append(li, s.LowSI)
		dsi = append(dsi, s.DeltaSI)
		gs = append(gs, s.Gstat)
		ed = append(ed, s.ED)
		lod = append(lod, s.LOD)
		bbl = append(bbl, s.BBLogBF)
	}

	const trim = 0.01 // exclude top 1 % when estimating background
	n := genomeWideNorms{}
	n.hiMed, n.hiMAD = robustBackground(hi, trim)
	n.liMed, n.liMAD = robustBackground(li, trim)
	n.dsiMed, n.dsiMAD = robustBackground(dsi, trim)
	n.gsMed, n.gsMAD = robustBackground(gs, trim)
	n.edMed, n.edMAD = robustBackground(ed, trim)
	n.lodMed, n.lodMAD = robustBackground(lod, trim)
	n.bblMed, n.bblMAD = robustBackground(bbl, trim)
	return n
}

// ---------------------------------------------------------------------------
// Shared chart style constants
// ---------------------------------------------------------------------------

const (
	chartTheme  = types.ThemeWesteros
	chartWidth  = "900px"
	chartHeight = "500px"

	// Significance Z thresholds for the overlay plots.
	zSig  = 3.0 // ~p99 equivalent
	zSugg = 2.0 // ~p95 equivalent

	defaultBRMAlpha = 0.05
)

// posFormatter is the shared Mb/kb x-axis label formatter (JS).
const posFormatter = `function(value) {
	if (value >= 1000000) { return (value / 1000000).toFixed(2) + ' Mb'; }
	if (value >= 1000)    { return (value / 1000).toFixed(1)    + ' kb'; }
	return value;
}`

// commonGlobalOpts returns the ECharts global option setters shared by every chart.
func commonGlobalOpts(title, subtitle, yLabel, width, height string, bidirectional bool) []charts.GlobalOpts {
	yMin := opts.Float(0.0)
	if bidirectional {
		yMin = nil // let ECharts auto-range for DeltaSI
	}
	return []charts.GlobalOpts{
		charts.WithInitializationOpts(opts.Initialization{
			Theme:  chartTheme,
			Width:  width,
			Height: height,
		}),
		charts.WithTitleOpts(opts.Title{
			Title:    title,
			Subtitle: subtitle,
			Left:     "center",
			Top:      "1%",
		}),
		charts.WithXAxisOpts(opts.XAxis{
			Name:         "Genomic Position",
			NameLocation: "middle",
			NameGap:      35,
			AxisLabel: &opts.AxisLabel{
				Rotate:    30,
				Formatter: opts.FuncOpts(posFormatter),
			},
		}),
		charts.WithYAxisOpts(opts.YAxis{
			Name:         yLabel,
			NameLocation: "middle",
			NameGap:      55,
			Min:          yMin,
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
		}),
		charts.WithDataZoomOpts(
			opts.DataZoom{Type: "slider", XAxisIndex: []int{0}, Start: 0, End: 100},
			opts.DataZoom{Type: "inside", XAxisIndex: []int{0}},
		),
		charts.WithLegendOpts(opts.Legend{
			Show:   opts.Bool(true),
			Top:    "9%",
			Left:   "center",
			Type:   "scroll",
			Orient: "horizontal",
		}),
		charts.WithToolboxOpts(opts.Toolbox{
			Show:  opts.Bool(true),
			Right: "2%",
			Feature: &opts.ToolBoxFeature{
				SaveAsImage: &opts.ToolBoxFeatureSaveAsImage{Show: opts.Bool(true), Title: "Save PNG"},
				DataZoom:    &opts.ToolBoxFeatureDataZoom{Show: opts.Bool(true), Title: map[string]string{"zoom": "Zoom", "back": "Reset"}},
				Restore:     &opts.ToolBoxFeatureRestore{Show: opts.Bool(true), Title: "Reset"},
			},
		}),
		charts.WithGridOpts(opts.Grid{
			Left:         "8%",
			Right:        "4%",
			Top:          "20%",
			Bottom:       "14%",
			ContainLabel: opts.Bool(true),
		}),
	}
}

// ---------------------------------------------------------------------------
// Individual (raw-value) line chart
// ---------------------------------------------------------------------------

// createInteractiveLineChart renders one statistic for one chromosome with
// permutation-threshold reference lines.  Pass hasNegativeThresh=true for
// DeltaSI to draw the valley thresholds as well.
func createInteractiveLineChart(
	title string,
	x []int64,
	y []float64,
	t99, t95 float64,
	tm99, tm95 float64,
	hasNegativeThresh bool,
	brmBlocks []BRMBlock,
) *charts.Line {

	subtitle := fmt.Sprintf("p99 threshold: %.4f  |  p95 threshold: %.4f  |  shaded: BRM blocks", t99, t95)

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Value", chartWidth, chartHeight, hasNegativeThresh)...)

	// Tooltip: flag values that cross the p99/p95 threshold.
	line.SetGlobalOptions(charts.WithTooltipOpts(opts.Tooltip{
		Show:        opts.Bool(true),
		Trigger:     "axis",
		AxisPointer: &opts.AxisPointer{Type: "cross"},
		Formatter: opts.FuncOpts(`function(params) {
			let pos = params[0].axisValue;
			let posFmt = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
			let result = '<strong>Position: ' + posFmt + '</strong><br/>';
			let t99val = null, t95val = null;
			params.forEach(function(p) {
				if (p.seriesName === 'p99') t99val = parseFloat(p.value);
				if (p.seriesName === 'p95') t95val = parseFloat(p.value);
			});
			params.forEach(function(item) {
				let val = parseFloat(item.value);
				if (isNaN(val)) return;
				let sig = '';
				if (item.seriesName === 'Statistic') {
					if (t99val !== null && val > t99val)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ p99</span>';
					else if (t95val !== null && val > t95val) sig = ' <span style="color:#f39c12">● p95</span>';
				}
				result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(5) + sig + '<br/>';
			});
			return result;
		}`),
	}))

	// Build data arrays.
	n := len(y)
	yData := make([]opts.LineData, n)
	y99 := make([]opts.LineData, n)
	y95 := make([]opts.LineData, n)
	var ym99, ym95 []opts.LineData
	if hasNegativeThresh {
		ym99 = make([]opts.LineData, n)
		ym95 = make([]opts.LineData, n)
	}
	for i, v := range y {
		yData[i] = opts.LineData{Value: v}
		y99[i] = opts.LineData{Value: t99}
		y95[i] = opts.LineData{Value: t95}
		if hasNegativeThresh {
			ym99[i] = opts.LineData{Value: tm99}
			ym95[i] = opts.LineData{Value: tm95}
		}
	}

	statOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#1f77b4"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	statOpts = append(statOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	line.SetXAxis(positionLabels(x)).
		AddSeries("Statistic", yData, statOpts...).
		AddSeries("p99", y99,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("p95", y95,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	if hasNegativeThresh {
		line.
			AddSeries("p99 valley", ym99,
				charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
				charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
			).
			AddSeries("p95 valley", ym95,
				charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
				charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
			)
	}

	return line
}

// ---------------------------------------------------------------------------
// Robust Z-score overlay chart (per chromosome, all stats)
// ---------------------------------------------------------------------------

// statColor returns the colorblind-friendly hex for each stat series name.
func statColor(name string) string {
	switch name {
	case "HighSI":
		return "#1f77b4" // blue
	case "LowSI":
		return "#ff7f0e" // orange
	case "DeltaSI":
		return "#2ca02c" // green  (thicker; primary stat)
	case "Gstat":
		return "#17becf" // cyan
	case "ED":
		return "#d62728" // red
	case "LOD":
		return "#9467bd" // purple
	case "BBLogBF":
		return "#8c564b" // brown
	}
	return "#7f7f7f"
}

// createRobustZOverlayChart builds an overlay of all seven robust Z-score
// series for one chromosome, with z=±2 and z=±3 reference lines.
func createRobustZOverlayChart(
	chrom string,
	x []int64,
	hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ []float64,
	brmBlocks []BRMBlock,
) *charts.Line {

	title := chrom + " — Robust Z-score Overlay"
	subtitle := "Genome-wide robust Z-score (background median+MAD, top 1% trimmed). " +
		"z = ±2 suggestive · z = ±3 significant. Shaded bands: BRM blocks."

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Robust Z-score", chartWidth, chartHeight, true)...)

	// Override y-axis for Z-score semantics and fix-point label at ±2, ±3.
	line.SetGlobalOptions(
		charts.WithYAxisOpts(opts.YAxis{
			Name:         "Robust Z-score",
			NameLocation: "middle",
			NameGap:      55,
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
			AxisLabel: &opts.AxisLabel{
				Formatter: opts.FuncOpts(`function(v) {
					let m = {3:'z=3 ★', 2:'z=2 ●', 0:'0', '-2':'z=-2 ●', '-3':'z=-3 ★'};
					let k = parseFloat(v.toFixed(1));
					return m[k] !== undefined ? m[k] : v.toFixed(1);
				}`),
			},
		}),
		charts.WithTooltipOpts(opts.Tooltip{
			Show:        opts.Bool(true),
			Trigger:     "axis",
			AxisPointer: &opts.AxisPointer{Type: "cross"},
			Formatter: opts.FuncOpts(`function(params) {
				let pos = params[0].axisValue;
				let posStr = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
				let result = '<strong>' + posStr + '</strong><br/>';
				let statSeries = ['HighSI','LowSI','DeltaSI','Gstat','ED','LOD','BBLogBF'];
				params.forEach(function(item) {
					if (statSeries.indexOf(item.seriesName) === -1) return;
					let val = parseFloat(item.value);
					if (isNaN(val)) return;
					let sig = '';
					if (Math.abs(val) >= 3.0)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ significant</span>';
					else if (Math.abs(val) >= 2.0)  sig = ' <span style="color:#f39c12">● suggestive</span>';
					result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
				});
				return result;
			}`),
		}),
	)

	n := len(x)

	// Threshold reference lines (flat horizontal series).
	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}
	zero := mkRef(0)
	z2p := mkRef(zSugg)
	z3p := mkRef(zSig)
	z2n := mkRef(-zSugg)
	z3n := mkRef(-zSig)

	zeroOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#bdc3c7"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	zeroOpts = append(zeroOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	// Reference lines first so they render behind data.
	line.SetXAxis(positionLabels(x)).
		AddSeries("z=0", zero, zeroOpts...).
		AddSeries("z=+2 (sugg.)", z2p,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=+3 (sig.)", z3p,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=-2 (sugg.)", z2n,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=-3 (sig.)", z3n,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	// Data series.
	type seriesDef struct {
		name  string
		data  []float64
		width float32
	}
	series := []seriesDef{
		{"HighSI", hiZ, 2.0},
		{"LowSI", liZ, 2.0},
		{"DeltaSI", dsiZ, 3.0}, // thicker — primary BSA stat
		{"Gstat", gsZ, 2.0},
		{"ED", edZ, 2.0},
		{"LOD", lodZ, 2.0},
		{"BBLogBF", bblZ, 2.0},
	}

	for _, s := range series {
		col := statColor(s.name)
		ld := floatSliceToLineData(s.data)
		line.AddSeries(s.name, ld,
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: s.width, Color: col}),
			charts.WithItemStyleOpts(opts.ItemStyle{Color: col, Opacity: opts.Float(0)}),
		)
	}

	return line
}

// brmBlockMarkAreaOpts converts BRM intervals into low-opacity vertical bands.
// The bands are attached to an existing series so no extra legend item is added.
func brmBlockMarkAreaOpts(blocks []BRMBlock, x []int64) []charts.SeriesOpts {
	if len(blocks) == 0 || len(x) == 0 {
		return nil
	}

	xLabels := positionLabels(x)
	areas := make([][]opts.MarkAreaData, 0, len(blocks))
	for _, b := range blocks {
		startIdx := sort.Search(len(x), func(i int) bool { return x[i] >= b.Start })
		stopIdx := sort.Search(len(x), func(i int) bool { return x[i] > b.Stop }) - 1
		if startIdx >= len(x) || stopIdx < 0 {
			continue
		}
		if stopIdx < startIdx {
			stopIdx = startIdx
		}

		areas = append(areas, []opts.MarkAreaData{
			{
				Name:  fmt.Sprintf("BRM block %.4f", b.Peak),
				XAxis: xLabels[startIdx],
			},
			{XAxis: xLabels[stopIdx]},
		})
	}
	if len(areas) == 0 {
		return nil
	}

	return []charts.SeriesOpts{
		charts.WithMarkAreaData(areas...),
		charts.WithMarkAreaStyleOpts(opts.MarkAreaStyle{
			Label:     &opts.Label{Show: opts.Bool(false)},
			ItemStyle: &opts.ItemStyle{Color: "rgba(243, 156, 18, 0.22)"},
		}),
	}
}

// ---------------------------------------------------------------------------
// Normalized overlay chart (per chromosome, threshold-relative — kept for
// backward compatibility but now clearly labelled)
// ---------------------------------------------------------------------------

// createNormalizedOverlayChart is retained for the legacy overlay page.
// Each statistic is divided by its per-chromosome average p99 threshold so
// that p99 = 1.0 and p95 ≈ 0.95.  See createRobustZOverlayChart for the
// depth-independent genome-wide version.
func createNormalizedOverlayChart(
	chrom string,
	x []int64,
	hi, li, dsi, gs, ed, lod, bbl []float64,
	avgHp99, avgHp95, avgLp99, avgLp95 float64,
	avgDp99, avgDp95, avgDMp99, avgDMp95 float64,
	avgGs99, avgGs95, avgEp99, avgEp95, avgLodp99, avgLodp95 float64,
	avgBbp99, avgBbp95 float64,
	brmBlocks []BRMBlock,
) *charts.Line {

	title := chrom + " — Threshold-Relative Overlay"
	subtitle := "Values divided by per-chromosome avg p99 threshold (p99=1.0; p95 varies by statistic). " +
		"Note: depth-dependent — see Robust Z-score page. Shaded bands: BRM blocks."

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Threshold-relative value", chartWidth, chartHeight, true)...)

	tooltipFormatter := fmt.Sprintf(`function(params) {
		let pos = params[0].axisValue;
		let posStr = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
		let result = '<strong>' + posStr + '</strong><br/>';
		let p95Pos = {
			'HighSI': %.6f,
			'LowSI': %.6f,
			'DeltaSI': %.6f,
			'Gstat': %.6f,
			'ED': %.6f,
			'LOD': %.6f,
			'BBLogBF': %.6f
		};
		let p95Neg = {'DeltaSI': %.6f};
		let statSeries = ['HighSI','LowSI','DeltaSI','Gstat','ED','LOD','BBLogBF'];
		params.forEach(function(item) {
			if (statSeries.indexOf(item.seriesName) === -1) return;
			let val = parseFloat(item.value);
			if (isNaN(val)) return;
			let p95 = val < 0 ? (p95Neg[item.seriesName] || p95Pos[item.seriesName] || 0) : (p95Pos[item.seriesName] || 0);
			let sig = '';
			if (Math.abs(val) >= 1.0) sig = ' <span style="color:#e74c3c;font-weight:bold">★ p99</span>';
			else if (p95 > 0 && Math.abs(val) >= p95) sig = ' <span style="color:#f39c12">● p95</span>';
			result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
		});
		return result;
	}`,
		thresholdRatio(avgHp95, avgHp99),
		thresholdRatio(avgLp95, avgLp99),
		thresholdRatio(avgDp95, avgDp99),
		thresholdRatio(avgGs95, avgGs99),
		thresholdRatio(avgEp95, avgEp99),
		thresholdRatio(avgLodp95, avgLodp99),
		thresholdRatio(avgBbp95, avgBbp99),
		thresholdRatio(math.Abs(avgDMp95), math.Abs(avgDMp99)),
	)

	// Override Y-axis labels to show threshold semantics.
	line.SetGlobalOptions(
		charts.WithYAxisOpts(opts.YAxis{
			Name:         "Threshold-relative value",
			NameLocation: "middle",
			NameGap:      55,
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
			AxisLabel: &opts.AxisLabel{
				Formatter: opts.FuncOpts(`function(v) {
					let fv = parseFloat(v.toFixed(3));
					if (fv === 1.0)   return 'p99 (+)';
					if (fv === -1.0)  return 'p99 (-)';
					if (fv === 0.0)   return '0';
					return v.toFixed(2);
				}`),
			},
		}),
		charts.WithTooltipOpts(opts.Tooltip{
			Show:        opts.Bool(true),
			Trigger:     "axis",
			AxisPointer: &opts.AxisPointer{Type: "cross"},
			Formatter:   opts.FuncOpts(tooltipFormatter),
		}),
	)

	n := len(x)
	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}

	baselineOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#bdc3c7"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	baselineOpts = append(baselineOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	line.SetXAxis(positionLabels(x)).
		AddSeries("z=0 baseline", mkRef(0), baselineOpts...).
		AddSeries("p99 (+)", mkRef(1.0),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("p99 (-)", mkRef(-1.0),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	type sd struct {
		name string
		data []float64
		w    float32
	}
	series := []sd{
		{"HighSI", normalizeToThreshold(hi, avgHp99, false), 2.0},
		{"LowSI", normalizeToThreshold(li, avgLp99, false), 2.0},
		{"DeltaSI", normalizeDeltaSI(dsi, avgDp99, avgDMp99), 3.0},
		{"Gstat", normalizeToThreshold(gs, avgGs99, false), 2.0},
		{"ED", normalizeToThreshold(ed, avgEp99, false), 2.0},
		{"LOD", normalizeToThreshold(lod, avgLodp99, false), 2.0},
		{"BBLogBF", normalizeToThreshold(bbl, avgBbp99, false), 2.0},
	}
	for _, s := range series {
		col := statColor(s.name)
		line.AddSeries(s.name, floatSliceToLineData(s.data),
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: s.w, Color: col}),
			charts.WithItemStyleOpts(opts.ItemStyle{Color: col, Opacity: opts.Float(0)}),
		)
	}

	return line
}

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

// floatSliceToLineData converts a []float64 to the ECharts line-data format.
func floatSliceToLineData(vals []float64) []opts.LineData {
	ld := make([]opts.LineData, len(vals))
	for i, v := range vals {
		ld[i] = opts.LineData{Value: v}
	}
	return ld
}

// positionLabels keeps the chart on a category x-axis while making markArea
// boundaries match actual category names instead of being interpreted as indexes.
func positionLabels(x []int64) []string {
	labels := make([]string, len(x))
	for i, v := range x {
		labels[i] = fmt.Sprintf("%d", v)
	}
	return labels
}

// normalizeToThreshold divides each value by the reference threshold.
// If invert is true the sign is flipped (for valley detection).
func normalizeToThreshold(vals []float64, ref float64, invert bool) []float64 {
	out := make([]float64, len(vals))
	if ref == 0 {
		return out
	}
	sign := 1.0
	if invert {
		sign = -1.0
	}
	for i, v := range vals {
		out[i] = sign * v / ref
	}
	return out
}

func thresholdRatio(numerator, denominator float64) float64 {
	if denominator == 0 {
		return 0
	}
	return math.Abs(numerator / denominator)
}

// normalizeDeltaSI maps positive DeltaSI values to [0, +1] using the positive
// p99 threshold and negative values to [0, -1] using the negative threshold.
func normalizeDeltaSI(dsi []float64, p99, mp99 float64) []float64 {
	out := make([]float64, len(dsi))
	for i, v := range dsi {
		if v >= 0 {
			if p99 != 0 {
				out[i] = v / p99
			}
		} else {
			if mp99 != 0 {
				out[i] = v / math.Abs(mp99) // preserves negative sign
			}
		}
	}
	return out
}

// detectQTLs identifies consecutive positions crossing a threshold and reports
// their interval bounds and peak value.
func detectQTLs(chrom string, x []int64, y []float64, threshold float64, statName, ci string, isValley bool) []QTLRecord {
	if len(x) == 0 || len(x) != len(y) || threshold == 0 {
		return nil
	}

	var qtls []QTLRecord
	inQTL := false
	var start, stop int64
	var peak float64

	for i, val := range y {
		condition := val > threshold
		if isValley {
			condition = val < threshold
		}

		if condition {
			if !inQTL {
				inQTL = true
				start = x[i]
				peak = val
				continue
			}
			if (isValley && val < peak) || (!isValley && val > peak) {
				peak = val
			}
			continue
		}

		if inQTL {
			stop = x[i-1]
			qtls = append(qtls, QTLRecord{
				Chrom: chrom,
				Start: start,
				Stop:  stop,
				Peak:  peak,
				Stat:  statName,
				CI:    ci,
			})
			inQTL = false
		}
	}

	if inQTL {
		stop = x[len(x)-1]
		qtls = append(qtls, QTLRecord{
			Chrom: chrom,
			Start: start,
			Stop:  stop,
			Peak:  peak,
			Stat:  statName,
			CI:    ci,
		})
	}
	return qtls
}

// calculateBRMBlocks applies BRM's block-threshold idea to the existing
// smoothed windows: AFD is DeltaSI, AFP is the mean of the two bulk SIs, and
// consecutive significant windows are merged into a single shaded interval.
func calculateBRMBlocks(chrom string, stats []SmoothedStats, highBulkSize, lowBulkSize, popLevel int, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || highBulkSize <= 0 || lowBulkSize <= 0 {
		return nil
	}

	n1 := float64(highBulkSize)
	n2 := float64(lowBulkSize)
	popScale := math.Pow(2, float64(popLevel))
	varianceScale := (n1 + n2) / (popScale * n1 * n2)

	var blocks []BRMBlock
	inBlock := false
	startIdx := 0
	peakIdx := 0
	peak := 0.0
	peakThreshold := 0.0

	emitBlock := func(startIdx, stopIdx, peakIdx int, peak, threshold float64) {
		start := stats[startIdx].POS
		if startIdx > 0 {
			start = (stats[startIdx-1].POS + stats[startIdx].POS) / 2
		}
		stop := stats[stopIdx].POS
		if stopIdx < len(stats)-1 {
			stop = (stats[stopIdx].POS + stats[stopIdx+1].POS) / 2
		}
		if stop < start {
			stop = start
		}
		blocks = append(blocks, BRMBlock{
			Chrom:     chrom,
			Start:     start,
			Stop:      stop,
			PeakPos:   stats[peakIdx].POS,
			Peak:      peak,
			Threshold: threshold,
		})
	}

	for i, s := range stats {
		afp := (s.HighSI + s.LowSI) / 2
		if afp < 0 {
			afp = 0
		}
		if afp > 1 {
			afp = 1
		}
		threshold := uAlpha * math.Sqrt(varianceScale*afp*(1-afp))
		significant := threshold > 0 && math.Abs(s.DeltaSI) >= threshold

		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = s.DeltaSI
				peakThreshold = threshold
				continue
			}
			if math.Abs(s.DeltaSI) > math.Abs(peak) {
				peakIdx = i
				peak = s.DeltaSI
				peakThreshold = threshold
			}
			continue
		}

		if inBlock {
			stopIdx := i - 1
			emitBlock(startIdx, stopIdx, peakIdx, peak, peakThreshold)
			inBlock = false
		}
	}

	if inBlock {
		stopIdx := len(stats) - 1
		emitBlock(startIdx, stopIdx, peakIdx, peak, peakThreshold)
	}
	return blocks
}

// writeHTMLPage creates (or truncates) a file and renders a page into it.
// It closes the file explicitly and returns any close error.
func writeHTMLPage(page *components.Page, path string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create %s: %w", path, err)
	}
	if err := page.Render(f); err != nil {
		_ = f.Close()
		return fmt.Errorf("render %s: %w", path, err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("close %s: %w", path, err)
	}
	return nil
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

// GenerateHtmlPlotsAndQTL processes smoothed window statistics, detects QTLs,
// and writes three HTML files plus a QTL TSV into the same directory as
// htmlOutFile.  The htmlOutFile name itself is not used — only its directory.
//
// Output files:
//
//	GoBSAseq_IndividualPlots.html   – raw-value charts with permutation thresholds
//	GoBSAseq_NormalizedOverlay.html – threshold-relative overlay (legacy method)
//	GoBSAseq_RobustZScore.html      – genome-wide robust Z-score overlay (recommended)
//	<qtlOutFile>                    – TSV of detected QTL intervals
//	GoBSAseq_BRMBlocks.tsv          – TSV of BRM-style block intervals used for plot shading
func GenerateHtmlPlotsAndQTL(
	allSmoothed []SmoothedStats,
	highSmAF, lowSmAF float64,
	highBulkSize, lowBulkSize int,
	population string,
	alphas []float64,
	rep int,
	htmlOutFile, qtlOutFile string,
) error {

	outDir := filepath.Dir(htmlOutFile)

	// ------------------------------------------------------------------
	// Pass 1 — genome-wide robust Z parameters (across all chromosomes)
	// ------------------------------------------------------------------
	norms := computeGenomeWideNorms(allSmoothed)

	// ------------------------------------------------------------------
	// Pass 2 — group by chromosome, build charts, detect QTLs
	// ------------------------------------------------------------------
	byChr := make(map[string][]SmoothedStats)
	for _, s := range allSmoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	var allQTLs []QTLRecord
	var allBRMBlocks []BRMBlock

	popLevel := 0
	if population == "F2" {
		popLevel = 1
	}
	brmAlpha := defaultBRMAlpha
	for _, alpha := range alphas {
		if alpha > brmAlpha && alpha < 1 {
			brmAlpha = alpha
		}
	}
	brmUAlpha := distuv.UnitNormal.Quantile(1 - brmAlpha/2)

	// Three separate pages.
	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics"

	normalizedPage := components.NewPage()
	normalizedPage.SetLayout(components.PageFlexLayout)
	normalizedPage.PageTitle = "GoBSAseq — Threshold-Relative Overlay"

	robustZPage := components.NewPage()
	robustZPage.SetLayout(components.PageFlexLayout)
	robustZPage.PageTitle = "GoBSAseq — Robust Z-score Overlay"

	for _, chrom := range chroms {
		stats := byChr[chrom]
		n := float64(len(stats))
		if n == 0 {
			continue
		}

		// Average permutation thresholds across windows of this chromosome.
		var (
			sumHp99, sumHp95     float64
			sumLp99, sumLp95     float64
			sumDp99, sumDp95     float64
			sumDMp99, sumDMp95   float64
			sumGs99, sumGs95     float64
			sumEp99, sumEp95     float64
			sumLodp99, sumLodp95 float64
			sumBbp99, sumBbp95   float64
		)
		for _, s := range stats {
			t := calcThresholdsCached(s.MeanHighBulkDP, s.MeanLowBulkDP, highSmAF, lowSmAF, rep)
			sumHp99 += t.HighP99
			sumHp95 += t.HighP95
			sumLp99 += t.LowP99
			sumLp95 += t.LowP95
			sumDp99 += t.DsiP99
			sumDp95 += t.DsiP95
			sumDMp99 += t.DsiMp99
			sumDMp95 += t.DsiMp95
			sumGs99 += t.GsP99
			sumGs95 += t.GsP95
			sumEp99 += t.EdP99
			sumEp95 += t.EdP95
			sumLodp99 += t.LodP99
			sumLodp95 += t.LodP95
			sumBbp99 += t.BbP99
			sumBbp95 += t.BbP95
		}
		avgHp99, avgHp95 := sumHp99/n, sumHp95/n
		avgLp99, avgLp95 := sumLp99/n, sumLp95/n
		avgDp99, avgDp95 := sumDp99/n, sumDp95/n
		avgDMp99, avgDMp95 := sumDMp99/n, sumDMp95/n
		avgGs99, avgGs95 := sumGs99/n, sumGs95/n
		avgEp99, avgEp95 := sumEp99/n, sumEp95/n
		avgLodp99, avgLodp95 := sumLodp99/n, sumLodp95/n
		avgBbp99, avgBbp95 := sumBbp99/n, sumBbp95/n

		// Extract per-chromosome data arrays.
		x := make([]int64, 0, len(stats))
		hi := make([]float64, 0, len(stats))
		li := make([]float64, 0, len(stats))
		dsi := make([]float64, 0, len(stats))
		gs := make([]float64, 0, len(stats))
		ed := make([]float64, 0, len(stats))
		lod := make([]float64, 0, len(stats))
		bbl := make([]float64, 0, len(stats))

		for _, s := range stats {
			x = append(x, s.POS)
			hi = append(hi, s.HighSI)
			li = append(li, s.LowSI)
			dsi = append(dsi, s.DeltaSI)
			gs = append(gs, s.Gstat)
			ed = append(ed, s.ED)
			lod = append(lod, s.LOD)
			bbl = append(bbl, s.BBLogBF)
		}

		// The permutation-based QTL TSV remains for continuity. Plot tinting now
		// uses BRM blocks, which are based on DeltaSI and pool-size variance.
		var chromQTLs []QTLRecord
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, hi, avgHp99, "HighSI", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, hi, avgHp95, "HighSI", "95", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, li, avgLp99, "LowSI", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, li, avgLp95, "LowSI", "95", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, dsi, avgDp99, "DeltaSI", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, dsi, avgDp95, "DeltaSI", "95", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, dsi, avgDMp99, "DeltaSI", "99", true)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, dsi, avgDMp95, "DeltaSI", "95", true)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, gs, avgGs99, "Gstat", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, gs, avgGs95, "Gstat", "95", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, ed, avgEp99, "ED", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, ed, avgEp95, "ED", "95", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, lod, avgLodp99, "LOD", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, lod, avgLodp95, "LOD", "95", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, bbl, avgBbp99, "BBLogBF", "99", false)...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, bbl, avgBbp95, "BBLogBF", "95", false)...)
		allQTLs = append(allQTLs, chromQTLs...)

		chromBRMBlocks := calculateBRMBlocks(chrom, stats, highBulkSize, lowBulkSize, popLevel, brmUAlpha)
		allBRMBlocks = append(allBRMBlocks, chromBRMBlocks...)

		// Page 1 — individual raw-value charts.
		individualPage.AddCharts(
			createInteractiveLineChart(chrom+" HighSI", x, hi, avgHp99, avgHp95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" LowSI", x, li, avgLp99, avgLp95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" DeltaSI", x, dsi, avgDp99, avgDp95, avgDMp99, avgDMp95, true, chromBRMBlocks),
			createInteractiveLineChart(chrom+" Gstat", x, gs, avgGs99, avgGs95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" ED", x, ed, avgEp99, avgEp95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" LOD", x, lod, avgLodp99, avgLodp95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" BBLogBF", x, bbl, avgBbp99, avgBbp95, 0, 0, false, chromBRMBlocks),
		)

		// Page 2 — threshold-relative overlay (legacy).
		normalizedPage.AddCharts(createNormalizedOverlayChart(
			chrom, x, hi, li, dsi, gs, ed, lod, bbl,
			avgHp99, avgHp95, avgLp99, avgLp95,
			avgDp99, avgDp95, avgDMp99, avgDMp95,
			avgGs99, avgGs95, avgEp99, avgEp95, avgLodp99, avgLodp95,
			avgBbp99, avgBbp95,
			chromBRMBlocks,
		))

		// Page 3 — genome-wide robust Z-score overlay.
		hiZ := robustZScore(hi, norms.hiMed, norms.hiMAD)
		liZ := robustZScore(li, norms.liMed, norms.liMAD)
		dsiZ := robustZScore(dsi, norms.dsiMed, norms.dsiMAD)
		gsZ := robustZScore(gs, norms.gsMed, norms.gsMAD)
		edZ := robustZScore(ed, norms.edMed, norms.edMAD)
		lodZ := robustZScore(lod, norms.lodMed, norms.lodMAD)
		bblZ := robustZScore(bbl, norms.bblMed, norms.bblMAD)
		robustZPage.AddCharts(createRobustZOverlayChart(chrom, x, hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))
	}

	// ------------------------------------------------------------------
	// Write HTML files
	// ------------------------------------------------------------------
	if err := writeHTMLPage(individualPage,
		filepath.Join(outDir, "GoBSAseq_IndividualPlots.html")); err != nil {
		return err
	}
	if err := writeHTMLPage(normalizedPage,
		filepath.Join(outDir, "GoBSAseq_NormalizedOverlay.html")); err != nil {
		return err
	}
	if err := writeHTMLPage(robustZPage,
		filepath.Join(outDir, "GoBSAseq_RobustZScore.html")); err != nil {
		return err
	}

	// ------------------------------------------------------------------
	// Write QTL TSV
	// ------------------------------------------------------------------
	fTsv, err := os.Create(qtlOutFile)
	if err != nil {
		return fmt.Errorf("create qtl file: %w", err)
	}
	fmt.Fprintf(fTsv, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\n")
	for _, q := range allQTLs {
		fmt.Fprintf(fTsv, "%s\t%d\t%d\t%.6f\t%s\t%s\n",
			q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI)
	}
	if err := fTsv.Close(); err != nil {
		return fmt.Errorf("close qtl file: %w", err)
	}

	fBRM, err := os.Create(filepath.Join(outDir, "GoBSAseq_BRMBlocks.tsv"))
	if err != nil {
		return fmt.Errorf("create brm blocks file: %w", err)
	}
	fmt.Fprintf(fBRM, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK_DELTA_SI\tBRM_THRESHOLD\n")
	for _, b := range allBRMBlocks {
		fmt.Fprintf(fBRM, "%s\t%d\t%d\t%d\t%.6f\t%.6f\n",
			b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.Threshold)
	}
	if err := fBRM.Close(); err != nil {
		return fmt.Errorf("close brm blocks file: %w", err)
	}

	return nil
}
