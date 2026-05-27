package oneBulk

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

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/go-echarts/go-echarts/v2/types"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const (
	minSNPsPerWindow  = 5
	maxGapWindows     = 3
	consensusMinStats = 3
	afpFloor          = 0.05
	chartTheme        = types.ThemeWesteros
	chartWidth        = "900px"
	chartHeight       = "500px"
	zSig              = 3.0 // ~p99 equivalent
	zSugg             = 2.0 // ~p95 equivalent
	defaultBRMAlpha   = 0.05
)

type OneBulkStats struct {
	CHROM            string
	POS              int64
	REF              string
	ALT              string
	HighParGT        []int
	LowParGT         []int
	BulkGT           []int
	BulkAD           string
	BulkSusAlleleCnt int
	BulkResAlleleCnt int
	SI               float64 // ALT / (ALT+REF)
	AbsSI            float64 // |SI - 0.5|
	Gstat            float64 // one-bulk G vs. uniform
	LOD              float64 // one-bulk LOD vs. p=0.5
	BBLogBF          float64 // BF vs. p=0.5
	ED               float64 // |SI-0.5|^4
	Depth            int
}

type Thresholds struct {
	SIP99  float64
	SIP95  float64
	SIMp99 float64
	SIMp95 float64

	AbsSIP99 float64
	AbsSIP95 float64

	GsP99 float64
	GsP95 float64

	EdP99 float64
	EdP95 float64

	LodP99 float64
	LodP95 float64

	BbP99 float64
	BbP95 float64
}

type SmoothedStats struct {
	CHROM      string
	POS        int64
	DeltaSI    float64
	Gstat      float64
	ED         float64
	LOD        float64
	BBLogBF    float64
	SI         float64
	AbsSI      float64
	NumSNPs    int
	MeanBulkDP int

	// per-window threshold lookup (set during smoothing, used in detectQTLs)
	thresholds Thresholds
}

// BRMBlock holds one one-bulk BRM-style segregation-deviation interval.
type BRMBlock struct {
	Chrom      string
	Start      int64
	Stop       int64
	PeakPos    int64
	Peak       float64
	ExpectedSI float64
	Threshold  float64
}

func smoothChromosome(stats []OneBulkStats, windowSize int64, step int64, bulkSmAF float64, rep int) []SmoothedStats {
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
			// DeltaSI, Gstat, LOD, BBLogBF share spatial+depth weight.
			sumGstat, sumWeightGs    float64
			sumLOD, sumWeightLod     float64
			sumBBLogBF, sumWeightBB  float64
			sumED, sumWeightED       float64
			sumSI, sumWeightSI       float64
			sumAbsSI, sumWeightAbsSI float64
			sumBulkDP                float64
			nSNPs                    int
		)

		for _, s := range stats {
			if s.POS < windowStart || s.POS > windowEnd {
				continue
			}
			nSNPs++

			d := math.Abs(float64(s.POS - center))
			w := utils.TricubeWeight(d, float64(windowSize)/2)
			depthWeight := math.Sqrt(float64(s.Depth))
			wStat := w * depthWeight

			sumGstat += s.Gstat * wStat
			sumWeightGs += wStat
			sumLOD += s.LOD * wStat
			sumWeightLod += wStat
			sumBBLogBF += s.BBLogBF * wStat
			sumWeightBB += wStat
			sumED += s.ED * wStat
			sumWeightED += wStat
			sumSI += s.SI * wStat
			sumWeightSI += wStat

			sumAbsSI += s.AbsSI * wStat
			sumWeightAbsSI += wStat
			sumBulkDP += float64(s.BulkResAlleleCnt + s.BulkSusAlleleCnt)

		}

		// Skip sparse windows — they produce unreliable signal.
		if nSNPs < minSNPsPerWindow {
			continue
		}

		sm := SmoothedStats{
			CHROM:      chrom,
			POS:        center,
			NumSNPs:    nSNPs,
			MeanBulkDP: int(sumBulkDP / float64(nSNPs)),
		}

		if sumWeightGs > 0 {
			sm.Gstat = sumGstat / sumWeightGs
		}
		if sumWeightLod > 0 {
			sm.LOD = sumLOD / sumWeightLod
		}
		if sumWeightBB > 0 {
			sm.BBLogBF = sumBBLogBF / sumWeightBB
		}
		if sumWeightED > 0 {
			sm.ED = sumED / sumWeightED
		}
		if sumWeightSI > 0 {
			sm.SI = sumSI / sumWeightSI

		}
		if sumWeightAbsSI > 0 {
			sm.AbsSI = sumAbsSI / sumWeightAbsSI
		}

		sm.thresholds = calcThresholdsCached(sm.MeanBulkDP, bulkSmAF, rep)
		smoothed = append(smoothed, sm)
	}

	return smoothed
}

func GStatisticOneBulk(gt []int) float64 {
	alt, ref := gt[0], gt[1]
	a := float64(alt) + 0.5
	r := float64(ref) + 0.5
	total := a + r
	exp := total / 2.0
	g := 0.0
	if a > 0 {
		g += a * math.Log(a/exp)
	}
	if r > 0 {
		g += r * math.Log(r/exp)
	}
	return 2 * g
}

func LodOneBulk(gt []int) float64 {
	alt, ref := gt[0], gt[0]
	total := alt + ref
	if total == 0 {
		return 0
	}
	const eps = 1e-10
	p := math.Max(eps, math.Min(1-eps, float64(alt)/float64(total)))
	p0 := 0.5
	return float64(total) * (p*math.Log(p/p0) + (1-p)*math.Log((1-p)/(1-p0)))
}

func logBeta(a, b float64) float64 {
	la, _ := math.Lgamma(a)
	lb, _ := math.Lgamma(b)
	lab, _ := math.Lgamma(a + b)
	return la + lb - lab
}

func BetaBinomialOneBulk(gt []int) float64 {
	alt, ref := gt[0], gt[1]
	a := float64(alt)
	r := float64(ref)
	total := a + r
	if total == 0 {
		return 0
	}
	// BF = logP(data | p estimated) - logP(data | p=0.5)
	logNull := total * math.Log(0.5)
	logAlt := logBeta(a+1, r+1) - logBeta(1, 1) // Beta(1,1) = uniform prior
	return logAlt - logNull
}

func calcThresholds(bulkDP int, bulkSmAF float64, rep int) Thresholds {
	if bulkDP <= 0 || rep <= 0 {
		return Thresholds{}
	}

	src := rand.NewSource(time.Now().UnixNano())
	rng := rand.New(src)
	dist := distuv.Binomial{N: float64(bulkDP), P: bulkSmAF, Src: rng}

	bulkSIArr := make([]float64, rep)

	absSiArr := make([]float64, rep)
	gsArr := make([]float64, rep)
	edArr := make([]float64, rep)
	lodArr := make([]float64, rep)
	bbArr := make([]float64, rep)

	for i := 0; i < rep; i++ {
		alt := dist.Rand()
		ref := float64(bulkDP) - alt

		si := math.Round((alt/float64(bulkDP))*1e6) / 1e6
		bulkSIArr[i] = si

		absSiArr[i] = math.Abs(si - 0.5) //math.Round((hSI-lSI)*1e6) / 1e6
		gsArr[i] = math.Round(GStatisticOneBulk([]int{int(alt), int(ref)})*1e6) / 1e6
		edArr[i] = math.Pow(math.Abs(si-0.5), 4)
		lodArr[i] = math.Round(LodOneBulk([]int{int(alt), int(ref)})*1e6) / 1e6
		bbArr[i] = math.Round(BetaBinomialOneBulk([]int{int(alt), int(ref)})*1e6) / 1e6
	}

	sort.Float64s(bulkSIArr)
	sort.Float64s(absSiArr)
	sort.Float64s(gsArr)
	sort.Float64s(edArr)
	sort.Float64s(lodArr)
	sort.Float64s(bbArr)

	return Thresholds{
		SIP99:  math.Round(stat.Quantile(0.995, stat.Empirical, bulkSIArr, nil)*1e6) / 1e6,
		SIP95:  math.Round(stat.Quantile(0.95, stat.Empirical, bulkSIArr, nil)*1e6) / 1e6,
		SIMp99: math.Round(stat.Quantile(0.005, stat.Empirical, bulkSIArr, nil)*1e6) / 1e6,
		SIMp95: math.Round(stat.Quantile(0.05, stat.Empirical, bulkSIArr, nil)*1e6) / 1e6,

		AbsSIP99: math.Round(stat.Quantile(0.995, stat.Empirical, absSiArr, nil)*1e6) / 1e6,
		AbsSIP95: math.Round(stat.Quantile(0.95, stat.Empirical, absSiArr, nil)*1e6) / 1e6,

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

var thresholdCache sync.Map

func calcThresholdsCached(bulkDP int, bulkSmAF float64, rep int) Thresholds {
	key := fmt.Sprintf("%d_%.6f_%d", bulkDP, bulkSmAF, rep)
	if v, ok := thresholdCache.Load(key); ok {
		return v.(Thresholds)
	}
	t := calcThresholds(bulkDP, bulkSmAF, rep)
	actual, _ := thresholdCache.LoadOrStore(key, t)
	return actual.(Thresholds)
}

func calcAllThresholds(allSmoothed []SmoothedStats, bulkSmAF float64, rep int) {
	type bulkDepth struct{ h int }
	seen := make(map[bulkDepth]bool)
	for _, sm := range allSmoothed {
		if sm.MeanBulkDP > 0 {
			seen[bulkDepth{sm.MeanBulkDP}] = true
		}
	}

	pairs := make([]bulkDepth, 0, len(seen))
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

	pairChan := make(chan bulkDepth, len(pairs))
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
				calcThresholdsCached(dp.h, bulkSmAF, rep)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println()
}

// calculateBRMBlocksOneBulk applies a BRM-style threshold to one-bulk windows.
// With no second bulk, the tested signal is the signed deviation of SI from the
// expected segregation frequency rather than a two-bulk allele-frequency
// difference. The variance model is therefore the one-sample analogue:
//
//	Var(SI - p0) = p0 * (1-p0) / (2^popLevel * bulkSize)
func calculateBRMBlocksOneBulk(chrom string, stats []SmoothedStats, bulkSize, popLevel int, expectedSI, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || bulkSize <= 0 || uAlpha <= 0 {
		return nil
	}

	p0 := expectedSI
	if math.IsNaN(p0) || math.IsInf(p0, 0) || p0 <= 0 || p0 >= 1 {
		p0 = 0.5
	}
	if p0 < afpFloor {
		p0 = afpFloor
	}
	if p0 > 1-afpFloor {
		p0 = 1 - afpFloor
	}

	n := float64(bulkSize)
	popScale := math.Pow(2, float64(popLevel))
	threshold := uAlpha * math.Sqrt((p0*(1-p0))/(popScale*n))
	if threshold <= 0 || math.IsNaN(threshold) || math.IsInf(threshold, 0) {
		return nil
	}

	var blocks []BRMBlock
	inBlock := false
	startIdx := 0
	peakIdx := 0
	peak := 0.0

	emitBlock := func(startIdx, stopIdx, peakIdx int, peak float64) {
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
			Chrom:      chrom,
			Start:      start,
			Stop:       stop,
			PeakPos:    stats[peakIdx].POS,
			Peak:       peak,
			ExpectedSI: p0,
			Threshold:  threshold,
		})
	}

	for i, s := range stats {
		deviation := s.SI - p0
		significant := math.Abs(deviation) >= threshold

		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = deviation
				continue
			}
			if math.Abs(deviation) > math.Abs(peak) {
				peakIdx = i
				peak = deviation
			}
			continue
		}

		if inBlock {
			emitBlock(startIdx, i-1, peakIdx, peak)
			inBlock = false
		}
	}

	if inBlock {
		emitBlock(startIdx, len(stats)-1, peakIdx, peak)
	}
	return blocks
}

func RunTwoParentsLowBulk(cfg utils.AnalysisConfig, hfcfg utils.HardFilterConfig) error {
	rdr := cfg.Rdr
	highParIdx := cfg.HighParentIdx
	lowParIdx := cfg.LowParentIdx
	lowBulkIdx := cfg.LowBulkIdx
	outDir := cfg.OutputDir

	windowSize := int64(cfg.WindowSize)
	stepSize := int64(cfg.StepSize)
	rep := cfg.Rep
	pop := cfg.Population

	overallStart := time.Now()

	//-------------------------------------- Remove problematic fields ---------------------------------------------- //
	for _, id := range []string{"PGT", "PID"} {
		delete(rdr.Header.SampleFormats, id)
	}

	//-------------------------------------- Header for writing ----------------------------------------------------- //
	fmt.Println(rdr.Header.SampleNames[highParIdx], rdr.Header.SampleNames[lowParIdx], rdr.Header.SampleNames[lowBulkIdx], highParIdx, lowParIdx, lowBulkIdx)
	sampleNames := []string{rdr.Header.SampleNames[highParIdx], rdr.Header.SampleNames[lowParIdx], rdr.Header.SampleNames[lowBulkIdx]}

	writerHeader := *rdr.Header
	writerHeader.SampleNames = sampleNames

	// ================================================= Run ======================================================== //

	// -------------------------------------------- Open output files ----------------------------------------------- //

	// --------------------------------------------- Raw BSAseq tsv ------------------------------------------------- //
	err := os.MkdirAll(filepath.Join(outDir, "stats"), 0755)
	if err != nil {
		return err
	}

	fmt.Println("Writing output to: ", filepath.Join(outDir, "stats"))
	rawFile := filepath.Join(outDir, "stats", "GoBSAseq.raw.tsv")
	rawHandle, err := os.Create(rawFile)
	if err != nil {
		return err
	}
	defer rawHandle.Close()

	rawWriter := bufio.NewWriter(rawHandle)
	defer rawWriter.Flush()

	_, err = fmt.Fprintln(rawWriter, "CHROM\tPOS\tREF\tALT\tHighParGT\tLowParGT\tBulkGT\tBulkAD\tSI\tAbsSI\tGstat\tED4\tLOD\tBBLogBF\tDepth")
	if err != nil {
		return err
	}

	// ----------------------------------------- Smoothed BSAseq tsv ------------------------------------------------ //

	smoothFile := filepath.Join(outDir, "stats", "GoBSAseq.smooth.tsv")
	smoothHandle, err := os.Create(smoothFile)
	if err != nil {
		return err
	}
	defer smoothHandle.Close()

	smoothWriter := bufio.NewWriter(smoothHandle)
	defer smoothWriter.Flush()

	header := "CHROM\tPOS\tSI\tAbsSI\tGstat\tED4\tLOD\tBBLogBF\tNumSNPs\tMeanBulkDP" +
		"\tSI_p99\tSI_p95\tSI_m_p99\tSI_m_p95" +
		"\tAbsSI_p99\tAbsSI_p95" +
		"\tGstat_p99\tGstat_p95" +
		"\tED4_p99\tED4_p95" +
		"\tLOD_p99\tLOD_p95" +
		"\tBBLogBF_p99\tBBLogBF_p95"
	fmt.Fprintln(smoothWriter, header)

	badVariant := 0
	chromStats := make(map[string][]OneBulkStats)
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		// ------------------------------------- Biallelic & non missing -------------------------------------------- //
		alts := v.Alt()
		if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
			badVariant++
			continue
		}

		passed := utils.PassesHardFilter(v, hfcfg) && v.Samples[cfg.HighParentIdx].DP >= cfg.HighParentDepth && v.Samples[cfg.LowParentIdx].DP >= cfg.LowParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
		if passed {
			highPar := v.Samples[cfg.HighParentIdx]
			lowPar := v.Samples[cfg.LowParentIdx]
			lowBulk := v.Samples[cfg.LowBulkIdx]

			bulkRefDep, _ := lowBulk.RefDepth()
			bulkAltDeps, _ := lowBulk.AltDepths()
			//fmt.Println(highPar.GT, lowPar.GT, lowBulk.GT, highPar.DP, lowPar.DP, lowBulk.DP)
			var bulkSusAlleleCount int
			var bulkResAlleleCount int
			if lowBulk.GT[0] == lowPar.GT[0] {
				bulkSusAlleleCount = bulkRefDep
				bulkResAlleleCount = bulkAltDeps[0]
			} else {
				bulkSusAlleleCount = bulkAltDeps[0]
				bulkResAlleleCount = bulkRefDep
			}

			SI := float64(bulkSusAlleleCount) / float64(lowBulk.DP)
			fmt.Println(SI)
			s := OneBulkStats{
				CHROM:            v.Chromosome,
				POS:              int64(v.Pos),
				REF:              v.Reference,
				ALT:              v.Alt()[0],
				HighParGT:        highPar.GT,
				LowParGT:         lowPar.GT,
				BulkGT:           lowBulk.GT,
				BulkSusAlleleCnt: bulkSusAlleleCount,
				BulkResAlleleCnt: bulkResAlleleCount,
				BulkAD:           fmt.Sprintf("%d,%d", bulkRefDep, bulkAltDeps[0]),
				SI:               SI,
				AbsSI:            math.Abs(SI - 0.5),
				ED:               math.Pow(math.Abs(SI-0.5), 4),
				Gstat:            math.Round(GStatisticOneBulk(lowBulk.GT)*1e6) / 1e6,
				LOD:              math.Round(LodOneBulk(lowBulk.GT)*1e6) / 1e6,
				BBLogBF:          math.Round(BetaBinomialOneBulk(lowBulk.GT)*1e6) / 1e6,
				Depth:            cfg.LowBulkDepth,
			}
			chromStats[s.CHROM] = append(chromStats[s.CHROM], s)
			fmt.Fprintf(rawWriter, "%s\t%d\t%s\t%s\t%v\t%v\t%v\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
				s.CHROM, s.POS, s.REF, s.ALT, s.HighParGT, s.LowParGT, s.BulkGT, s.BulkAD,
				s.SI, s.AbsSI, s.Gstat, s.ED, s.LOD, s.BBLogBF, s.Depth)
		}
	}

	// ---------------------------------------------------------------------------------------------------------
	// Stage 2 — smoothing
	// ---------------------------------------------------------------------------------------------------------
	color.Cyan("\n============================ Smoothing Statistics =============================\n\n")
	bulkSmAF := utils.SimulateAF(pop, float64(cfg.LowBulkSize), rep)
	popLevel := 0
	if pop == "F2" {
		popLevel = 1
	}
	brmAlpha := defaultBRMAlpha
	for _, alpha := range cfg.Alphas {
		if alpha > brmAlpha && alpha < 1 {
			brmAlpha = alpha
		}
	}
	brmUAlpha := distuv.UnitNormal.Quantile(1 - brmAlpha/2)

	var allSmoothed []SmoothedStats
	var allBRMBlocks []BRMBlock
	color.Cyan("\n============================ Smoothing & Calculating Thresholds (%d simulations per depth pair) ==============================\n\n", rep)
	for chrom, stats := range chromStats {
		color.Yellow("Smoothing %s: %d SNPs", chrom, len(stats))
		smoothed := smoothChromosome(stats, windowSize, stepSize, bulkSmAF, rep)
		for _, d := range smoothed {
			t := d.thresholds
			fmt.Fprintf(smoothWriter,
				"%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d"+
					"\t%.6f\t%.6f\t%.6f\t%.6f"+
					"\t%.6f\t%.6f"+
					"\t%.6f\t%.6f"+
					"\t%.6f\t%.6f"+
					"\t%.6f\t%.6f"+
					"\t%.6f\t%.6f\n",
				d.CHROM, d.POS,
				d.SI, d.AbsSI, d.Gstat, d.ED, d.LOD, d.BBLogBF,
				d.NumSNPs, d.MeanBulkDP,
				t.SIP99, t.SIP95, t.SIMp99, t.SIMp95,
				t.AbsSIP99, t.AbsSIP95,
				t.GsP99, t.GsP95,
				t.EdP99, t.EdP95,
				t.LodP99, t.LodP95,
				t.BbP99, t.BbP95,
			)
		}
		allSmoothed = append(allSmoothed, smoothed...)
		allBRMBlocks = append(allBRMBlocks, calculateBRMBlocksOneBulk(chrom, smoothed, cfg.LowBulkSize, popLevel, bulkSmAF, brmUAlpha)...)
	}
	color.Green("\nSmoothing, threshold calculations, and smoothed TSV complete.")
	color.Green("Raw stats written to %s", filepath.Join(outDir, "GoBSAseq.raw.tsv"))
	color.Green("\nTotal time: %s\n", time.Since(overallStart).Round(time.Second))

	// -----------------------------------------------------------------------
	// Stage 3 — BRM analysis
	// -----------------------------------------------------------------------
	sort.Slice(allBRMBlocks, func(i, j int) bool {
		if allBRMBlocks[i].Chrom == allBRMBlocks[j].Chrom {
			return allBRMBlocks[i].Start < allBRMBlocks[j].Start
		}
		return allBRMBlocks[i].Chrom < allBRMBlocks[j].Chrom
	})

	brmFile := filepath.Join(outDir, "stats", "GoBSAseq.onebulk.BRMBlocks.tsv")
	brmHandle, err := os.Create(brmFile)
	if err != nil {
		return fmt.Errorf("create one-bulk brm blocks file: %w", err)
	}
	brmWriter := bufio.NewWriter(brmHandle)
	fmt.Fprintln(brmWriter, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK_SI_DEVIATION\tEXPECTED_SI\tBRM_THRESHOLD\tVALIDATION")
	for _, b := range allBRMBlocks {
		validation := "PASS"
		if math.Abs(b.Peak) < b.Threshold {
			validation = "FAIL"
		}
		fmt.Fprintf(brmWriter, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%s\n",
			b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.ExpectedSI, b.Threshold, validation)
	}
	if err := brmWriter.Flush(); err != nil {
		_ = brmHandle.Close()
		return fmt.Errorf("flush one-bulk brm blocks file: %w", err)
	}
	if err := brmHandle.Close(); err != nil {
		return fmt.Errorf("close one-bulk brm blocks file: %w", err)
	}
	color.Green("One-bulk BRM-style blocks written to %s", brmFile)

	return nil

}
