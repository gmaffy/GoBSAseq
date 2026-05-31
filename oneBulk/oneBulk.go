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
	"strings"
	"sync"
	"time"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/genespace"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
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
	SI               float64 // resistant-parent allele count / observed bulk allele depth
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

	// Genome-wide robust Z-scores (set after all chromosomes are smoothed)
	SIZ      float64
	AbsSIZ   float64
	GstatZ   float64
	EDZ      float64
	LODZ     float64
	BBLogBFZ float64

	// max-|Z| composite signal
	CompositeZ float64

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

// QTLRecord holds the detected QTL interval and its peak value.
type QTLRecord struct {
	Chrom  string
	Start  int64
	Stop   int64
	Peak   float64
	Stat   string
	CI     string
	Source string
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

func GStatisticOneBulk(alt, ref float64) float64 {
	//alt, ref := gt[1], gt[0]
	a := alt + 0.5
	r := ref + 0.5
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

func LodOneBulk(alt, ref float64) float64 {
	//alt, ref := gt[1], gt[0]
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

func BetaBinomialOneBulk(alt, ref float64) float64 {
	//alt, ref := gt[1], gt[0]

	total := alt + ref
	if total == 0 {
		return 0
	}
	// BF = logP(data | p estimated) - logP(data | p=0.5)
	logNull := total * math.Log(0.5)
	logAlt := logBeta(alt+1, ref+1) - logBeta(1, 1) // Beta(1,1) = uniform prior
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
		gsArr[i] = math.Round(GStatisticOneBulk(alt, ref)*1e6) / 1e6
		edArr[i] = math.Pow(math.Abs(si-0.5), 4)
		lodArr[i] = math.Round(LodOneBulk(alt, ref)*1e6) / 1e6
		bbArr[i] = math.Round(BetaBinomialOneBulk(alt, ref)*1e6) / 1e6
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

// ---------------------------------------------------------------------------
// Robust Z-score helpers (adapted from twoBulk)
// ---------------------------------------------------------------------------

// zQuantile returns the p-th quantile of a pre-sorted slice.
func zQuantile(sorted []float64, p float64) float64 {
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

// robustBackground computes the median and MAD of vals, excluding the top
// trimFrac proportion of values so that genuine QTL peaks do not inflate
// the spread estimate.
func robustBackground(vals []float64, trimFrac float64) (median, mad float64) {
	if len(vals) == 0 {
		return 0, 0
	}
	sorted := make([]float64, len(vals))
	copy(sorted, vals)
	sort.Float64s(sorted)

	cutIdx := int(math.Round(float64(len(sorted)) * (1.0 - trimFrac)))
	if cutIdx < 1 {
		cutIdx = 1
	}
	trimmed := sorted[:cutIdx]

	median = zQuantile(trimmed, 0.5)

	devs := make([]float64, len(trimmed))
	for i, v := range trimmed {
		devs[i] = math.Abs(v - median)
	}
	sort.Float64s(devs)
	mad = zQuantile(devs, 0.5)
	return median, mad
}

// robustZScore normalises vals to genome-wide robust Z-scores.
// The 1.4826 factor makes the MAD-based scale consistent with σ under normality.
func robustZScoreSlice(vals []float64, median, mad float64) []float64 {
	out := make([]float64, len(vals))
	scale := mad * 1.4826
	if scale == 0 {
		return out
	}
	for i, v := range vals {
		out[i] = (v - median) / scale
	}
	return out
}

// computeZScores calculates genome-wide robust Z-scores for each statistic
// across all smoothed windows and assigns them (plus a max-|Z| composite)
// back onto each SmoothedStats entry in-place.
func computeZScores(allSmoothed []SmoothedStats) {
	if len(allSmoothed) == 0 {
		return
	}
	const trim = 0.01

	collect := func(fn func(SmoothedStats) float64) []float64 {
		v := make([]float64, len(allSmoothed))
		for i, s := range allSmoothed {
			v[i] = fn(s)
		}
		return v
	}

	siMed, siMAD := robustBackground(collect(func(s SmoothedStats) float64 { return s.SI }), trim)
	absMed, absMAD := robustBackground(collect(func(s SmoothedStats) float64 { return s.AbsSI }), trim)
	gsMed, gsMAD := robustBackground(collect(func(s SmoothedStats) float64 { return s.Gstat }), trim)
	edMed, edMAD := robustBackground(collect(func(s SmoothedStats) float64 { return s.ED }), trim)
	lodMed, lodMAD := robustBackground(collect(func(s SmoothedStats) float64 { return s.LOD }), trim)
	bblMed, bblMAD := robustBackground(collect(func(s SmoothedStats) float64 { return s.BBLogBF }), trim)

	siZ := robustZScoreSlice(collect(func(s SmoothedStats) float64 { return s.SI }), siMed, siMAD)
	absZ := robustZScoreSlice(collect(func(s SmoothedStats) float64 { return s.AbsSI }), absMed, absMAD)
	gsZ := robustZScoreSlice(collect(func(s SmoothedStats) float64 { return s.Gstat }), gsMed, gsMAD)
	edZ := robustZScoreSlice(collect(func(s SmoothedStats) float64 { return s.ED }), edMed, edMAD)
	lodZ := robustZScoreSlice(collect(func(s SmoothedStats) float64 { return s.LOD }), lodMed, lodMAD)
	bblZ := robustZScoreSlice(collect(func(s SmoothedStats) float64 { return s.BBLogBF }), bblMed, bblMAD)

	for i := range allSmoothed {
		allSmoothed[i].SIZ = siZ[i]
		allSmoothed[i].AbsSIZ = absZ[i]
		allSmoothed[i].GstatZ = gsZ[i]
		allSmoothed[i].EDZ = edZ[i]
		allSmoothed[i].LODZ = lodZ[i]
		allSmoothed[i].BBLogBFZ = bblZ[i]
		allSmoothed[i].CompositeZ = math.Max(math.Abs(siZ[i]),
			math.Max(math.Abs(absZ[i]),
				math.Max(math.Abs(gsZ[i]),
					math.Max(math.Abs(edZ[i]),
						math.Max(math.Abs(lodZ[i]), math.Abs(bblZ[i]))))))
	}
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

	// -----------------------------------------------------------------------
	// Stage 0 — hard filtering
	// -----------------------------------------------------------------------
	color.Cyan("============================ GATK Hard Filtering (One bulk two parents) ============================\n\n")

	filteredVcfPath := filepath.Join(outDir, "GoBSAseq.hard_filtered.vcf.gz")
	badVcfPath := filepath.Join(outDir, "GoBSAseq.bad_variants.vcf.gz")

	passedVariants, original, hardFiltered, err := utils.HardFilterVcf(rdr, filteredVcfPath, badVcfPath, cfg, hfcfg, 2)
	if err != nil {
		return fmt.Errorf("hard filter error: %w", err)
	}

	// Update indices to match subsetted samples: [HighParent, LowParent, LowBulk].
	cfg.HighParentIdx = 0
	cfg.LowParentIdx = 1
	cfg.LowBulkIdx = 2
	highParIdx = 0
	lowParIdx = 1
	lowBulkIdx = 2

	color.Green("Original variants: %v\nHard filtered variants: %v", original, hardFiltered)

	//-------------------------------------- Remove problematic fields ---------------------------------------------- //
	for _, id := range []string{"PGT", "PID"} {
		delete(rdr.Header.SampleFormats, id)
	}

	//-------------------------------------- Header for writing ----------------------------------------------------- //
	fmt.Println(cfg.HighParentName, cfg.LowParentName, cfg.LowBulkName, highParIdx, lowParIdx, lowBulkIdx)
	sampleNames := []string{cfg.HighParentName, cfg.LowParentName, cfg.LowBulkName}

	writerHeader := *rdr.Header
	writerHeader.SampleNames = sampleNames

	// ================================================= Run ======================================================== //

	// -------------------------------------------- Open output files ----------------------------------------------- //

	// --------------------------------------------- Raw BSAseq tsv ------------------------------------------------- //
	err = os.MkdirAll(filepath.Join(outDir, "stats"), 0755)
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
		"\tBBLogBF_p99\tBBLogBF_p95" +
		"\tSI_Z\tAbsSI_Z\tGstat_Z\tED4_Z\tLOD_Z\tBBLogBF_Z\tComposite_Z"
	fmt.Fprintln(smoothWriter, header)

	badVariant := 0
	chromStats := make(map[string][]OneBulkStats)
	for _, v := range passedVariants {
		// ------------------------------------- Biallelic & non missing -------------------------------------------- //
		alts := v.Alt()
		if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
			badVariant++
			continue
		}

		realAltIdx := -1
		for i, alt := range alts {
			if !(alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<')) {
				realAltIdx = i
				break
			}
		}
		if realAltIdx < 0 {
			continue
		}

		highPar := v.Samples[cfg.HighParentIdx]
		lowPar := v.Samples[cfg.LowParentIdx]
		lowBulk := v.Samples[cfg.LowBulkIdx]

		bulkRefDep, _ := lowBulk.RefDepth()
		bulkAltDeps, _ := lowBulk.AltDepths()
		if len(bulkAltDeps) <= realAltIdx {
			continue
		}
		//fmt.Println(highPar.GT, lowPar.GT, lowBulk.GT, highPar.DP, lowPar.DP, lowBulk.DP)
		var bulkSusAlleleCount int
		var bulkResAlleleCount int
		if lowPar.GT[0] == 0 {
			bulkSusAlleleCount = bulkRefDep
			bulkResAlleleCount = bulkAltDeps[realAltIdx]
		} else {
			bulkSusAlleleCount = bulkAltDeps[realAltIdx]
			bulkResAlleleCount = bulkRefDep
		}

		observedDepth := bulkSusAlleleCount + bulkResAlleleCount
		if observedDepth == 0 {
			continue
		}

		SI := float64(bulkResAlleleCount) / float64(observedDepth)
		//fmt.Println(SI)
		s := OneBulkStats{
			CHROM:            v.Chromosome,
			POS:              int64(v.Pos),
			REF:              v.Reference,
			ALT:              v.Alt()[realAltIdx],
			HighParGT:        highPar.GT,
			LowParGT:         lowPar.GT,
			BulkGT:           lowBulk.GT,
			BulkSusAlleleCnt: bulkSusAlleleCount,
			BulkResAlleleCnt: bulkResAlleleCount,
			BulkAD:           fmt.Sprintf("%d,%d", bulkRefDep, bulkAltDeps[realAltIdx]),
			SI:               SI,
			AbsSI:            math.Abs(SI - 0.5),
			ED:               math.Pow(math.Abs(SI-0.5), 4),
			Gstat:            math.Round(GStatisticOneBulk(float64(bulkResAlleleCount), float64(bulkSusAlleleCount))*1e6) / 1e6,
			LOD:              math.Round(LodOneBulk(float64(bulkResAlleleCount), float64(bulkSusAlleleCount))*1e6) / 1e6,
			BBLogBF:          math.Round(BetaBinomialOneBulk(float64(bulkResAlleleCount), float64(bulkSusAlleleCount))*1e6) / 1e6,
			Depth:            observedDepth,
		}
		chromStats[s.CHROM] = append(chromStats[s.CHROM], s)
		fmt.Fprintf(rawWriter, "%s\t%d\t%s\t%s\t%v\t%v\t%v\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
			s.CHROM, s.POS, s.REF, s.ALT, s.HighParGT, s.LowParGT, s.BulkGT, s.BulkAD,
			s.SI, s.AbsSI, s.Gstat, s.ED, s.LOD, s.BBLogBF, s.Depth)
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
	brmAlpha := cfg.BrmAlpha
	brmUAlpha := distuv.UnitNormal.Quantile(1 - brmAlpha/2)

	var allSmoothed []SmoothedStats
	var allBRMBlocks []BRMBlock
	color.Cyan("\n============================ Smoothing & Calculating Thresholds (%d simulations per depth pair) ==============================\n\n", rep)
	for chrom, stats := range chromStats {
		color.Yellow("Smoothing %s: %d SNPs", chrom, len(stats))
		smoothed := smoothChromosome(stats, windowSize, stepSize, bulkSmAF, rep)
		allSmoothed = append(allSmoothed, smoothed...)
		allBRMBlocks = append(allBRMBlocks, calculateBRMBlocksOneBulk(chrom, smoothed, cfg.LowBulkSize, popLevel, bulkSmAF, brmUAlpha)...)
	}
	color.Green("\nSmoothing, threshold calculations, and smoothed TSV complete.")

	// Compute genome-wide robust Z-scores and composite Z across all windows.
	color.Cyan("\n============================ Computing Robust Z-scores =============================\n\n")
	computeZScores(allSmoothed)
	color.Green("Z-score normalisation complete.")

	// Write smoothed TSV (z-scores are now populated).
	color.Cyan("\n============================ Writing Smoothed TSV =============================\n\n")
	for _, d := range allSmoothed {
		t := d.thresholds
		fmt.Fprintf(smoothWriter,
			"%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d"+
				"\t%.6f\t%.6f\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f"+
				"\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
			d.CHROM, d.POS,
			d.SI, d.AbsSI, d.Gstat, d.ED, d.LOD, d.BBLogBF,
			d.NumSNPs, d.MeanBulkDP,
			t.SIP99, t.SIP95, t.SIMp99, t.SIMp95,
			t.AbsSIP99, t.AbsSIP95,
			t.GsP99, t.GsP95,
			t.EdP99, t.EdP95,
			t.LodP99, t.LodP95,
			t.BbP99, t.BbP95,
			d.SIZ, d.AbsSIZ, d.GstatZ, d.EDZ, d.LODZ, d.BBLogBFZ, d.CompositeZ,
		)
	}
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

	// ---------------------------------------------------------------------------------------------------------
	// Stage 4: Plotting and QTL detection
	// ---------------------------------------------------------------------------------------------------------
	color.Cyan("\n============================ Generating HTML Plots & QTLs ========================================\n\n")

	finalQTLs, err := GenerateHtmlPlotsAndQTLOneBulk(
		allSmoothed, bulkSmAF, cfg.LowBulkSize, pop, cfg.Alphas, rep, outDir,
	)
	if err != nil {
		color.Red("Error generating Plots and QTLs: %v", err)
		return fmt.Errorf("plots and QTLs error: %w", err)
	}
	color.Green("HTML plots written to %s", filepath.Join(outDir, "plots"))
	color.Green("QTL tabular results written to %s", filepath.Join(outDir, "GoBSAseq_QTL.tsv"))

	// ---------------------------------------------------------------------------------------------------------
	// Stage 5: Gene space analysis
	// ---------------------------------------------------------------------------------------------------------
	color.Cyan("\n============================ Performing Gene Space Analysis ========================================\n\n")

	geneSpaceEnabled := cfg.SnpEffDB != "" && cfg.GeneDesc != "" && cfg.Prg != "" && cfg.Gff != ""
	if !geneSpaceEnabled {
		color.Yellow("Skipping gene space analysis: required parameters were not provided.")
		return nil
	}

	color.Green("Step 1: Annotating genes with SnpEff (creating Super VCF table)")
	_, hasEFF := cfg.Rdr.Header.Infos["EFF"]
	if hasEFF {
		color.Yellow("EFF column is already present; no Super VCF table was generated, so gene space analysis will be skipped.")
		return nil
	}
	color.Yellow("EFF column is not present; annotating filtered VCF.")

	err, annotatedTsvFiles := annotation.CreateSuperVcf([]string{filteredVcfPath}, cfg.SnpEffDB, true, cfg.GeneDesc, cfg.Prg)
	if err != nil {
		color.Red("Failed variant annotation with SNPEFF: %s", err)
		return fmt.Errorf("snpeff annotation: %w", err)
	}
	if len(annotatedTsvFiles) == 0 {
		color.Yellow("Skipping gene space analysis: SnpEff annotation did not produce an annotated TSV.")
		return nil
	}
	color.Green("SNPEFF annotation complete.")
	color.Green("Annotated TSV files: %v", annotatedTsvFiles)

	color.Cyan("Performing Gene space analysis for %d QTL intervals ...", len(finalQTLs))
	for _, qtl := range finalQTLs {
		color.Blue("GeneSpace interval: %s:%d-%d", qtl.Chrom, qtl.Start, qtl.Stop)
		_, err = genespace.GeneSpace(
			cfg.Gff, annotatedTsvFiles[0], qtl.Chrom,
			int(qtl.Start), int(qtl.Stop), []string{cfg.HighParentName},
			[]string{cfg.LowParentName}, cfg.GeneDesc, cfg.Prg, outDir)
		if err != nil {
			color.Red("Failed gene space analysis: %s", err)
			return fmt.Errorf("gene space: %w", err)
		}
	}
	color.Green("Gene space analysis complete.")

	return nil
}

// ---------------------------------------------------------------------------
// QTL Detection Helper Functions
// ---------------------------------------------------------------------------

func detectQTLs(chrom string, x []int64, y []float64, threshold float64, statName, ci string, isValley bool, source string) []QTLRecord {
	if len(x) == 0 || len(x) != len(y) || threshold == 0 {
		return nil
	}
	type run struct {
		start, stop int64
		peak        float64
	}
	var runs []run
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
			} else if (!isValley && val > peak) || (isValley && val < peak) {
				peak = val
			}
			stop = x[i]
			continue
		}
		if inQTL {
			runs = append(runs, run{start, stop, peak})
			inQTL = false
		}
	}
	if inQTL {
		runs = append(runs, run{start, stop, peak})
	}
	if len(runs) == 0 {
		return nil
	}

	merged := []run{runs[0]}
	for i := 1; i < len(runs); i++ {
		prev := &merged[len(merged)-1]
		prevStopIdx := sort.Search(len(x), func(j int) bool { return x[j] >= prev.stop })
		nextStartIdx := sort.Search(len(x), func(j int) bool { return x[j] >= runs[i].start })
		if nextStartIdx-prevStopIdx-1 <= maxGapWindows {
			prev.stop = runs[i].stop
			if (!isValley && runs[i].peak > prev.peak) || (isValley && runs[i].peak < prev.peak) {
				prev.peak = runs[i].peak
			}
		} else {
			merged = append(merged, runs[i])
		}
	}

	qtls := make([]QTLRecord, 0, len(merged))
	for _, r := range merged {
		qtls = append(qtls, QTLRecord{Chrom: chrom, Start: r.start, Stop: r.stop, Peak: r.peak, Stat: statName, CI: ci, Source: source})
	}
	return qtls
}

func detectQTLsAdaptive(chrom string, x []int64, y, thresholds []float64, statName, ci string, isValley bool, source string) []QTLRecord {
	if len(x) == 0 || len(x) != len(y) || len(x) != len(thresholds) {
		return nil
	}
	type run struct {
		start, stop int64
		peak        float64
	}
	var runs []run
	inQTL := false
	var start, stop int64
	var peak float64

	for i, val := range y {
		t := thresholds[i]
		if t == 0 {
			if inQTL {
				runs = append(runs, run{start, stop, peak})
				inQTL = false
			}
			continue
		}
		condition := val > t
		if isValley {
			condition = val < t
		}
		if condition {
			if !inQTL {
				inQTL = true
				start = x[i]
				peak = val
			} else if (!isValley && val > peak) || (isValley && val < peak) {
				peak = val
			}
			stop = x[i]
			continue
		}
		if inQTL {
			runs = append(runs, run{start, stop, peak})
			inQTL = false
		}
	}
	if inQTL {
		runs = append(runs, run{start, stop, peak})
	}
	if len(runs) == 0 {
		return nil
	}

	merged := []run{runs[0]}
	for i := 1; i < len(runs); i++ {
		prev := &merged[len(merged)-1]
		prevStopIdx := sort.Search(len(x), func(j int) bool { return x[j] >= prev.stop })
		nextStartIdx := sort.Search(len(x), func(j int) bool { return x[j] >= runs[i].start })
		if nextStartIdx-prevStopIdx-1 <= maxGapWindows {
			prev.stop = runs[i].stop
			if (!isValley && runs[i].peak > prev.peak) || (isValley && runs[i].peak < prev.peak) {
				prev.peak = runs[i].peak
			}
		} else {
			merged = append(merged, runs[i])
		}
	}

	qtls := make([]QTLRecord, 0, len(merged))
	for _, r := range merged {
		qtls = append(qtls, QTLRecord{Chrom: chrom, Start: r.start, Stop: r.stop, Peak: r.peak, Stat: statName, CI: ci, Source: source})
	}
	return qtls
}

func detectConsensusQTLsOneBulk(chrom string, stats []SmoothedStats) []QTLRecord {
	type hit struct {
		pos   int64
		count int
		fired []string
	}
	var hits []hit
	for _, s := range stats {
		t := s.thresholds
		var fired []string
		if s.SI > t.SIP99 || s.SI < t.SIMp99 {
			fired = append(fired, "SI")
		}
		if s.AbsSI > t.AbsSIP99 {
			fired = append(fired, "AbsSI")
		}
		if s.Gstat > t.GsP99 {
			fired = append(fired, "Gstat")
		}
		if s.ED > t.EdP99 {
			fired = append(fired, "ED4")
		}
		if s.LOD > t.LodP99 {
			fired = append(fired, "LOD")
		}
		if s.BBLogBF > t.BbP99 {
			fired = append(fired, "BBLogBF")
		}
		if len(fired) >= consensusMinStats {
			hits = append(hits, hit{s.POS, len(fired), fired})
		}
	}
	if len(hits) == 0 {
		return nil
	}

	var qtls []QTLRecord
	start := hits[0].pos
	stop := hits[0].pos
	maxCount := hits[0].count
	uniqueStats := make(map[string]bool)
	for _, st := range hits[0].fired {
		uniqueStats[st] = true
	}

	emit := func() {
		allFired := make([]string, 0, len(uniqueStats))
		for st := range uniqueStats {
			allFired = append(allFired, st)
		}
		sort.Strings(allFired)
		qtls = append(qtls, QTLRecord{
			Chrom: chrom, Start: start, Stop: stop, Peak: float64(maxCount),
			Stat: "Consensus", CI: "99", Source: strings.Join(allFired, ","),
		})
	}

	for i := 1; i < len(hits); i++ {
		prevIdx := sort.Search(len(stats), func(j int) bool { return stats[j].POS >= stop })
		nextIdx := sort.Search(len(stats), func(j int) bool { return stats[j].POS >= hits[i].pos })
		if nextIdx-prevIdx-1 <= maxGapWindows {
			stop = hits[i].pos
			if hits[i].count > maxCount {
				maxCount = hits[i].count
			}
			for _, st := range hits[i].fired {
				uniqueStats[st] = true
			}
		} else {
			emit()
			start = hits[i].pos
			stop = hits[i].pos
			maxCount = hits[i].count
			uniqueStats = make(map[string]bool)
			for _, st := range hits[i].fired {
				uniqueStats[st] = true
			}
		}
	}
	emit()
	return qtls
}

func intersectQTLsWithBRM(qtls []QTLRecord, brm []BRMBlock, targetSource string) []QTLRecord {
	var hc []QTLRecord
	for _, q := range qtls {
		for _, b := range brm {
			if q.Chrom == b.Chrom && q.Start <= b.Stop && q.Stop >= b.Start {
				hc = append(hc, QTLRecord{
					Chrom: q.Chrom, Start: q.Start, Stop: q.Stop, Peak: q.Peak,
					Stat: q.Stat, CI: q.CI, Source: targetSource,
				})
				break
			}
		}
	}
	return hc
}

// ---------------------------------------------------------------------------
// Charting / Plotting Helper Functions
// ---------------------------------------------------------------------------

func commonGlobalOpts(title, subtitle, yLabel string, hasNegativeThresh bool) []charts.GlobalOpts {
	yMin := opts.Float(0)
	if hasNegativeThresh {
		yMin = nil
	}
	posFormatter := `function(v) {
		if (v >= 1e9) return (v/1e9).toFixed(2)+' Gb';
		if (v >= 1e6) return (v/1e6).toFixed(2)+' Mb';
		if (v >= 1e3) return (v/1e3).toFixed(1)+' kb';
		return v;
	}`
	return []charts.GlobalOpts{
		charts.WithInitializationOpts(opts.Initialization{
			Theme:  chartTheme,
			Width:  chartWidth,
			Height: chartHeight,
		}),
		charts.WithTitleOpts(opts.Title{
			Title:    title,
			Subtitle: subtitle,
			Left:     "center",
			Top:      "2%",
		}),
		charts.WithXAxisOpts(opts.XAxis{
			Name:         "Genomic Position",
			NameLocation: "middle",
			NameGap:      35,
			AxisLabel:    &opts.AxisLabel{Rotate: 30, Formatter: opts.FuncOpts(posFormatter)},
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
		charts.WithLegendOpts(opts.Legend{Show: opts.Bool(true), Top: "9%", Left: "center", Type: "scroll"}),
		charts.WithToolboxOpts(opts.Toolbox{
			Show:  opts.Bool(true),
			Right: "2%",
			Feature: &opts.ToolBoxFeature{
				SaveAsImage: &opts.ToolBoxFeatureSaveAsImage{Show: opts.Bool(true)},
				DataZoom:    &opts.ToolBoxFeatureDataZoom{Show: opts.Bool(true)},
				Restore:     &opts.ToolBoxFeatureRestore{Show: opts.Bool(true)},
			},
		}),
		charts.WithGridOpts(opts.Grid{Left: "8%", Right: "4%", Top: "20%", Bottom: "14%", ContainLabel: opts.Bool(true)}),
	}
}

func createInteractiveLineChart(title string, x []int64, y []float64, t99, t95, tm99, tm95 float64, hasNegativeThresh bool, brmBlocks []BRMBlock) *charts.Line {
	subtitle := fmt.Sprintf("p99 threshold: %.4f  |  p95 threshold: %.4f  |  shaded: BRM blocks", t99, t95)
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Value", hasNegativeThresh)...)

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
		AddSeries("p99", y99, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"})).
		AddSeries("p95", y95, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}))

	if hasNegativeThresh {
		line.AddSeries("p99 valley", ym99, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"})).
			AddSeries("p95 valley", ym95, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}))
	}
	return line
}

func createRobustZOverlayChartOneBulk(chrom string, x []int64, siZ, absZ, gsZ, edZ, lodZ, bblZ []float64, brmBlocks []BRMBlock) *charts.Line {
	title := chrom + " — Robust Z-score Overlay"
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, "Single-bulk robust Z-scores. z=±2 suggestive, z=±3 significant.", "Robust Z-score", true)...)

	n := len(x)
	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}
	zeroOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#bdc3c7"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	zeroOpts = append(zeroOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	line.SetXAxis(positionLabels(x)).
		AddSeries("z=0", mkRef(0), zeroOpts...).
		AddSeries("z=+3", mkRef(zSig), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#e74c3c"})).
		AddSeries("z=-3", mkRef(-zSig), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#e74c3c"}))

	type seriesDef struct {
		name string
		data []float64
		col  string
	}
	for _, s := range []seriesDef{
		{"SI", siZ, "#1f77b4"},
		{"AbsSI", absZ, "#ff7f0e"},
		{"Gstat", gsZ, "#17becf"},
		{"ED4", edZ, "#d62728"},
		{"LOD", lodZ, "#9467bd"},
		{"BBLogBF", bblZ, "#8c564b"},
	} {
		line.AddSeries(s.name, floatSliceToLineData(s.data),
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: 2, Color: s.col}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)
	}
	return line
}

func createCompositeSignalChartOneBulk(chrom string, x []int64, siZ, absZ, gsZ, edZ, lodZ, bblZ []float64, brmBlocks []BRMBlock) *charts.Line {
	n := len(x)
	composite := make([]float64, n)
	for i := range composite {
		composite[i] = math.Max(math.Abs(siZ[i]),
			math.Max(math.Abs(absZ[i]),
				math.Max(math.Abs(gsZ[i]),
					math.Max(math.Abs(edZ[i]),
						math.Max(math.Abs(lodZ[i]), math.Abs(bblZ[i]))))))
	}
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(chrom+" — goplot", "Max |Z| across all statistics.", "max |Z-score|", false)...)

	compositeData := floatSliceToLineData(composite)
	compositeOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#2ca02c"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	compositeOpts = append(compositeOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}
	line.SetXAxis(positionLabels(x)).
		AddSeries("Composite", compositeData, compositeOpts...).
		AddSeries("z=2", mkRef(zSugg), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#f39c12"})).
		AddSeries("z=3", mkRef(zSig), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#e74c3c"}))
	return line
}

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
		areas = append(areas, []opts.MarkAreaData{{XAxis: xLabels[startIdx]}, {XAxis: xLabels[stopIdx]}})
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

func floatSliceToLineData(vals []float64) []opts.LineData {
	ld := make([]opts.LineData, len(vals))
	for i, v := range vals {
		ld[i] = opts.LineData{Value: v}
	}
	return ld
}

func positionLabels(x []int64) []string {
	labels := make([]string, len(x))
	for i, v := range x {
		labels[i] = fmt.Sprintf("%d", v)
	}
	return labels
}

func writeHTMLPage(page *components.Page, path string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create %s: %w", path, err)
	}
	if err := page.Render(f); err != nil {
		_ = f.Close()
		return fmt.Errorf("render %s: %w", path, err)
	}
	return f.Close()
}

// ---------------------------------------------------------------------------
// Main Plot + QTL Entry Point (Single Bulk)
// ---------------------------------------------------------------------------

func GenerateHtmlPlotsAndQTLOneBulk(
	allSmoothed []SmoothedStats,
	bulkSmAF float64,
	bulkSize int,
	population string,
	alphas []float64,
	rep int,
	outDir string,
) ([]QTLRecord, error) {

	plotsDir := filepath.Join(outDir, "plots")
	if err := os.MkdirAll(plotsDir, 0755); err != nil {
		return nil, fmt.Errorf("create plots directory: %w", err)
	}

	const trim = 0.01
	collectStat := func(fn func(SmoothedStats) float64) []float64 {
		v := make([]float64, len(allSmoothed))
		for i, s := range allSmoothed {
			v[i] = fn(s)
		}
		return v
	}

	siMed, siMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.SI }), trim)
	absMed, absMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.AbsSI }), trim)
	gsMed, gsMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.Gstat }), trim)
	edMed, edMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.ED }), trim)
	lodMed, lodMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.LOD }), trim)
	bblMed, bblMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.BBLogBF }), trim)

	byChr := make(map[string][]SmoothedStats)
	for _, s := range allSmoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	var allQTLs, allConsensusQTLs, allMaxZQTLs []QTLRecord
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

	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics (Single Bulk)"

	robustZPage := components.NewPage()
	robustZPage.SetLayout(components.PageFlexLayout)
	robustZPage.PageTitle = "GoBSAseq — Robust Z-score Overlay (Single Bulk)"

	compositePage := components.NewPage()
	compositePage.SetLayout(components.PageFlexLayout)
	compositePage.PageTitle = "goplot"

	for _, chrom := range chroms {
		stats := byChr[chrom]
		if len(stats) == 0 {
			continue
		}

		nf := float64(len(stats))
		var (
			sumSp99, sumSp95, sumSMp99, sumSMp95 float64
			sumAp99, sumAp95                     float64
			sumGs99, sumGs95                     float64
			sumEp99, sumEp95                     float64
			sumLod99, sumLod95                   float64
			sumBb99, sumBb95                     float64
		)
		for _, s := range stats {
			t := s.thresholds
			sumSp99 += t.SIP99
			sumSp95 += t.SIP95
			sumSMp99 += t.SIMp99
			sumSMp95 += t.SIMp95
			sumAp99 += t.AbsSIP99
			sumAp95 += t.AbsSIP95
			sumGs99 += t.GsP99
			sumGs95 += t.GsP95
			sumEp99 += t.EdP99
			sumEp95 += t.EdP95
			sumLod99 += t.LodP99
			sumLod95 += t.LodP95
			sumBb99 += t.BbP99
			sumBb95 += t.BbP95
		}

		n := len(stats)
		x := make([]int64, n)
		si := make([]float64, n)
		abs := make([]float64, n)
		gs := make([]float64, n)
		ed := make([]float64, n)
		lod := make([]float64, n)
		bbl := make([]float64, n)
		siT99, siTM99 := make([]float64, n), make([]float64, n)
		siT95, siTM95 := make([]float64, n), make([]float64, n)
		absT99, absT95 := make([]float64, n), make([]float64, n)
		gsT99, gsT95 := make([]float64, n), make([]float64, n)
		edT99, edT95 := make([]float64, n), make([]float64, n)
		lodT99, lodT95 := make([]float64, n), make([]float64, n)
		bblT99, bblT95 := make([]float64, n), make([]float64, n)

		for i, s := range stats {
			x[i] = s.POS
			si[i], abs[i], gs[i], ed[i], lod[i], bbl[i] = s.SI, s.AbsSI, s.Gstat, s.ED, s.LOD, s.BBLogBF
			t := s.thresholds
			siT99[i], siTM99[i] = t.SIP99, t.SIMp99
			siT95[i], siTM95[i] = t.SIP95, t.SIMp95
			absT99[i], absT95[i] = t.AbsSIP99, t.AbsSIP95
			gsT99[i], gsT95[i] = t.GsP99, t.GsP95
			edT99[i], edT95[i] = t.EdP99, t.EdP95
			lodT99[i], lodT95[i] = t.LodP99, t.LodP95
			bblT99[i], bblT95[i] = t.BbP99, t.BbP95
		}

		var chromQTLs []QTLRecord
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siT99, "SI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siTM99, "SI", "99", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siT95, "SI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siTM95, "SI", "95", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, abs, absT99, "AbsSI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, abs, absT95, "AbsSI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, gs, gsT99, "Gstat", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, gs, gsT95, "Gstat", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, ed, edT99, "ED4", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, ed, edT95, "ED4", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, lod, lodT99, "LOD", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, lod, lodT95, "LOD", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, bbl, bblT99, "BBLogBF", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, bbl, bblT95, "BBLogBF", "95", false, "Permutation")...)

		siZ := robustZScoreSlice(si, siMed, siMAD)
		absZ := robustZScoreSlice(abs, absMed, absMAD)
		gsZ := robustZScoreSlice(gs, gsMed, gsMAD)
		edZ := robustZScoreSlice(ed, edMed, edMAD)
		lodZ := robustZScoreSlice(lod, lodMed, lodMAD)
		bblZ := robustZScoreSlice(bbl, bblMed, bblMAD)

		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, siZ, zSig, "SI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, siZ, -zSig, "SI_Z", "z3", true, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, absZ, zSig, "AbsSI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, gsZ, zSig, "Gstat_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, edZ, zSig, "ED4_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, lodZ, zSig, "LOD_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, bblZ, zSig, "BBLogBF_Z", "z3", false, "ZScore")...)

		cQTLs := detectConsensusQTLsOneBulk(chrom, stats)
		allConsensusQTLs = append(allConsensusQTLs, cQTLs...)

		chromBRMBlocks := calculateBRMBlocksOneBulk(chrom, stats, bulkSize, popLevel, bulkSmAF, brmUAlpha)
		allBRMBlocks = append(allBRMBlocks, chromBRMBlocks...)

		composite := make([]float64, n)
		for i := range composite {
			composite[i] = math.Max(math.Abs(siZ[i]),
				math.Max(math.Abs(absZ[i]),
					math.Max(math.Abs(gsZ[i]),
						math.Max(math.Abs(edZ[i]),
							math.Max(math.Abs(lodZ[i]), math.Abs(bblZ[i]))))))
		}
		maxZQTLs := detectQTLs(chrom, x, composite, zSig, "Composite_Z", "z3", false, "MaxZ")
		allMaxZQTLs = append(allMaxZQTLs, maxZQTLs...)
		allMaxZQTLs = append(allMaxZQTLs, intersectQTLsWithBRM(maxZQTLs, chromBRMBlocks, "CompositeHighConfidence")...)

		chromQTLs = append(chromQTLs, intersectQTLsWithBRM(chromQTLs, chromBRMBlocks, "HighConfidence")...)
		allQTLs = append(allQTLs, chromQTLs...)

		robustZPage.AddCharts(createRobustZOverlayChartOneBulk(chrom, x, siZ, absZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))
		individualPage.AddCharts(
			createInteractiveLineChart(chrom+" SI", x, si, sumSp99/nf, sumSp95/nf, sumSMp99/nf, sumSMp95/nf, true, chromBRMBlocks),
			createInteractiveLineChart(chrom+" AbsSI", x, abs, sumAp99/nf, sumAp95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" Gstat", x, gs, sumGs99/nf, sumGs95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" ED4", x, ed, sumEp99/nf, sumEp95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" LOD", x, lod, sumLod99/nf, sumLod95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" BBLogBF", x, bbl, sumBb99/nf, sumBb95/nf, 0, 0, false, chromBRMBlocks),
		)
		compositePage.AddCharts(createCompositeSignalChartOneBulk(chrom, x, siZ, absZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))
	}

	if err := writeHTMLPage(individualPage, filepath.Join(plotsDir, "GoBSAseq_IndividualPlots.html")); err != nil {
		return nil, err
	}
	if err := writeHTMLPage(robustZPage, filepath.Join(plotsDir, "GoBSAseq_RobustZScore.html")); err != nil {
		return nil, err
	}
	if err := writeHTMLPage(compositePage, filepath.Join(plotsDir, "GoBSAseq_final.html")); err != nil {
		return nil, err
	}

	qtlOutFile := filepath.Join(outDir, "GoBSAseq_QTL.tsv")
	fTsv, err := os.Create(qtlOutFile)
	if err != nil {
		return nil, fmt.Errorf("create qtl file: %w", err)
	}
	fmt.Fprintf(fTsv, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\tSOURCE\n")
	for _, q := range allQTLs {
		fmt.Fprintf(fTsv, "%s\t%d\t%d\t%.6f\t%s\t%s\t%s\n", q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI, q.Source)
	}
	_ = fTsv.Close()

	fCons, err := os.Create(filepath.Join(outDir, "GoBSAseq_QTL_CONSENSUS.tsv"))
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(fCons, "CHROM\tSTART\tSTOP\t#STATS\tSTAT\tCI\tSTATS\n")
	for _, q := range allConsensusQTLs {
		fmt.Fprintf(fCons, "%s\t%d\t%d\t%d\t%s\t%s\t%s\n", q.Chrom, q.Start, q.Stop, int(q.Peak), q.Stat, q.CI, q.Source)
	}
	_ = fCons.Close()

	fMaxZ, err := os.Create(filepath.Join(outDir, "GoBSAseq_QTL_MAX_Z.tsv"))
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(fMaxZ, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\tSOURCE\n")
	for _, q := range allMaxZQTLs {
		fmt.Fprintf(fMaxZ, "%s\t%d\t%d\t%.6f\t%s\t%s\t%s\n", q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI, q.Source)
	}
	_ = fMaxZ.Close()

	fBRM, err := os.Create(filepath.Join(outDir, "GoBSAseq_BRMBlocks.tsv"))
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(fBRM, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK_ABS_SI\tBRM_THRESHOLD\n")
	for _, b := range allBRMBlocks {
		fmt.Fprintf(fBRM, "%s\t%d\t%d\t%d\t%.6f\t%.6f\n", b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.Threshold)
	}
	_ = fBRM.Close()

	return allMaxZQTLs, nil
}
