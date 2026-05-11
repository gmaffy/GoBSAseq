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
// Constants
// ---------------------------------------------------------------------------

const (
	// minSNPsPerWindow is the minimum number of SNPs required for a smoothed
	// window to be emitted.  Windows below this threshold (e.g. pericentromeric
	// deserts) produce unstable averages and can generate spurious QTL peaks.
	minSNPsPerWindow = 5

	// maxGapWindows controls gap-bridging in detectQTLs: two QTL runs separated
	// by this many or fewer sub-threshold windows are merged into a single record.
	maxGapWindows = 3

	// consensusMinStats is the minimum number of statistics that must
	// simultaneously exceed their thresholds for a window to be called a
	// consensus QTL.
	consensusMinStats = 3

	// afpFloor is the minimum allele-frequency product used in calculateBRMBlocks
	// to prevent the variance threshold from approaching zero near fixation.
	afpFloor = 0.05

	chartTheme  = types.ThemeWesteros
	chartWidth  = "900px"
	chartHeight = "500px"

	// Significance Z thresholds for the overlay plots.
	zSig  = 3.0 // ~p99 equivalent
	zSugg = 2.0 // ~p95 equivalent

	defaultBRMAlpha = 0.05
)

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
	ED      float64 // ED^4 power (Magwene et al.) – distinct from |DeltaSI|
	LOD     float64
	BBLogBF float64

	Depth int
}

// HardFilterStats summarises the outcome of a hard-filter pass.
type HardFilterStats struct {
	Total        int
	Passed       int
	SNPs         int
	INDELs       int
	SNPPass      int
	INDELPass    int
	FilterCounts map[string]int // per-filter-name fail counts
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

	// per-window threshold lookup (set during smoothing, used in detectQTLs)
	thresholds Thresholds
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
	// Source identifies the detection method: "Permutation", "ZScore", "Consensus", "HighConfidence"
	Source string
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

// ================================== BSA-seq statistics ======================================================== //

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

// euclideanDistance4 computes the ED^4 formulation from Magwene et al. (2011),
// which provides a statistic that is genuinely distinct from |DeltaSI|:
//
//	ED = sqrt((hSI - lSI)^2)  →  ED^4 = (hSI - lSI)^4
//
// Raising to the 4th power amplifies large differences while strongly
// suppressing noise near zero, increasing separation between signal and background.
func euclideanDistance4(hSI, lSI float64) float64 {
	d := hSI - lSI
	return d * d * d * d
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
// Smoothing and TSV output
// ---------------------------------------------------------------------------

func tricubeWeight(d, D float64) float64 {
	if D <= 0 || d >= D {
		return 0
	}
	x := 1 - math.Pow(d/D, 3)
	return x * x * x
}

// smoothChromosome computes weighted sliding-window averages over sorted BSA
// statistics.  Each statistic uses its own normalised weight accumulator so
// that Gstat, LOD, and BBLogBF are not biased by the DeltaSI weight
// distribution.  Windows with fewer than minSNPsPerWindow SNPs are skipped to
// avoid unstable estimates in low-density genomic regions.
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
			// DeltaSI, Gstat, LOD, BBLogBF share spatial+depth weight.
			sumDeltaSI, sumWeightDSI float64
			sumGstat, sumWeightGs    float64
			sumLOD, sumWeightLod     float64
			sumBBLogBF, sumWeightBB  float64
			// ED has its own accumulator since it is now ED^4.
			sumED, sumWeightED float64
			// HighSI / LowSI use their own accumulator.
			sumHighSI, sumLowSI, sumWeightSI float64

			sumHighDP, sumLowDP float64
			nSNPs               int
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
			sumWeightDSI += wStat
			sumGstat += s.Gstat * wStat
			sumWeightGs += wStat
			sumLOD += s.LOD * wStat
			sumWeightLod += wStat
			sumBBLogBF += s.BBLogBF * wStat
			sumWeightBB += wStat
			sumED += s.ED * wStat
			sumWeightED += wStat
			sumHighSI += s.HighSI * wStat
			sumLowSI += s.LowSI * wStat
			sumWeightSI += wStat
			sumHighDP += float64(s.HighBulkL + s.HighBulkH)
			sumLowDP += float64(s.LowBulkL + s.LowBulkH)
		}

		// Skip sparse windows — they produce unreliable signal.
		if nSNPs < minSNPsPerWindow {
			continue
		}

		sm := SmoothedStats{
			CHROM:          chrom,
			POS:            center,
			NumSNPs:        nSNPs,
			MeanHighBulkDP: int(sumHighDP / float64(nSNPs)),
			MeanLowBulkDP:  int(sumLowDP / float64(nSNPs)),
		}
		if sumWeightDSI > 0 {
			sm.DeltaSI = sumDeltaSI / sumWeightDSI
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

	fmt.Fprintln(w, "CHROM\tPOS\tREF\tALT\tHighParGT\tLowParGT\tHighBulkGT\tHighBulkAD\tLowBulkGT\tLowBulkAD\tHighBulkL\tHighBulkH\tLowBulkL\tLowBulkH\tHighSI\tLowSI\tDeltaSI\tGstat\tED4\tLOD\tBBLogBF\tDepth")
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

	header := "CHROM\tPOS\tHighSI\tLowSI\tDeltaSI\tGstat\tED4\tLOD\tBBLogBF\tNumSNPs\tMeanHighDP\tMeanLowDP" +
		"\tHighSI_p99\tHighSI_p95\tHighSI_m_p99\tHighSI_m_p95" +
		"\tLowSI_p99\tLowSI_p95\tLowSI_m_p99\tLowSI_m_p95" +
		"\tDeltaSI_p99\tDeltaSI_p95\tDeltaSI_m_p99\tDeltaSI_m_p95" +
		"\tGstat_p99\tGstat_p95" +
		"\tED4_p99\tED4_p95" +
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

	cutIdx := int(math.Round(float64(len(sorted)) * (1.0 - trimFrac)))
	if cutIdx < 1 {
		cutIdx = 1
	}
	trimmed := sorted[:cutIdx]

	median = quantile(trimmed, 0.5)

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
		return out
	}
	for i, v := range vals {
		out[i] = (v - median) / scale
	}
	return out
}



// ---------------------------------------------------------------------------
// QTL detection
// ---------------------------------------------------------------------------

// detectQTLs identifies consecutive windows exceeding a threshold and merges
// runs separated by ≤ maxGapWindows sub-threshold windows (gap-bridging).
// When isValley is true the sense is inverted (for DeltaSI valley calls).
// source labels the QTL records for traceability in the output TSV.
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
			} else {
				if (!isValley && val > peak) || (isValley && val < peak) {
					peak = val
				}
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

	// Gap-bridging: merge runs separated by ≤ maxGapWindows windows.
	merged := []run{runs[0]}
	for i := 1; i < len(runs); i++ {
		prev := &merged[len(merged)-1]
		// Find the index of prev.stop and runs[i].start in x to count the gap.
		prevStopIdx := sort.Search(len(x), func(j int) bool { return x[j] >= prev.stop })
		nextStartIdx := sort.Search(len(x), func(j int) bool { return x[j] >= runs[i].start })
		gap := nextStartIdx - prevStopIdx - 1
		if gap <= maxGapWindows {
			// Merge: extend the previous run and keep the more extreme peak.
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
		qtls = append(qtls, QTLRecord{
			Chrom:  chrom,
			Start:  r.start,
			Stop:   r.stop,
			Peak:   r.peak,
			Stat:   statName,
			CI:     ci,
			Source: source,
		})
	}
	return qtls
}

// detectConsensusQTLs finds windows where ≥ consensusMinStats statistics
// simultaneously exceed their per-window thresholds and reports the merged
// genomic intervals.
func detectConsensusQTLs(chrom string, stats []SmoothedStats) []QTLRecord {
	type hit struct {
		pos   int64
		count int
	}

	hits := make([]hit, 0, len(stats))
	for _, s := range stats {
		t := s.thresholds
		count := 0
		if s.HighSI > t.HighP99 || s.HighSI < t.HighMp99 {
			count++
		}
		if s.LowSI > t.LowP99 || s.LowSI < t.LowMp99 {
			count++
		}
		if s.DeltaSI > t.DsiP99 || s.DeltaSI < t.DsiMp99 {
			count++
		}
		if s.Gstat > t.GsP99 {
			count++
		}
		if s.ED > t.EdP99 {
			count++
		}
		if s.LOD > t.LodP99 {
			count++
		}
		if s.BBLogBF > t.BbP99 {
			count++
		}
		if count >= consensusMinStats {
			hits = append(hits, hit{s.POS, count})
		}
	}

	if len(hits) == 0 {
		return nil
	}

	// Merge consecutive hit positions (gap-bridge by maxGapWindows implicitly
	// because hits with count < threshold are naturally absent).
	var qtls []QTLRecord
	start := hits[0].pos
	stop := hits[0].pos
	maxCount := hits[0].count

	for i := 1; i < len(hits); i++ {
		// Use index distance in the full stats slice as proxy for gap size.
		prevIdx := sort.Search(len(stats), func(j int) bool { return stats[j].POS >= stop })
		nextIdx := sort.Search(len(stats), func(j int) bool { return stats[j].POS >= hits[i].pos })
		gap := nextIdx - prevIdx - 1
		if gap <= maxGapWindows {
			stop = hits[i].pos
			if hits[i].count > maxCount {
				maxCount = hits[i].count
			}
		} else {
			qtls = append(qtls, QTLRecord{
				Chrom:  chrom,
				Start:  start,
				Stop:   stop,
				Peak:   float64(maxCount),
				Stat:   "Consensus",
				CI:     "99",
				Source: "Consensus",
			})
			start = hits[i].pos
			stop = hits[i].pos
			maxCount = hits[i].count
		}
	}
	qtls = append(qtls, QTLRecord{
		Chrom:  chrom,
		Start:  start,
		Stop:   stop,
		Peak:   float64(maxCount),
		Stat:   "Consensus",
		CI:     "99",
		Source: "Consensus",
	})
	return qtls
}

// intersectQTLsWithBRM marks permutation-called QTLs that overlap at least one
// BRM block as high-confidence and returns a deduplicated list.
func intersectQTLsWithBRM(qtls []QTLRecord, brm []BRMBlock) []QTLRecord {
	var hc []QTLRecord
	for _, q := range qtls {
		if q.Source == "Consensus" || q.Source == "HighConfidence" {
			continue
		}
		for _, b := range brm {
			if b.Chrom == q.Chrom && b.Start <= q.Stop && b.Stop >= q.Start {
				hc = append(hc, QTLRecord{
					Chrom:  q.Chrom,
					Start:  q.Start,
					Stop:   q.Stop,
					Peak:   q.Peak,
					Stat:   q.Stat,
					CI:     q.CI,
					Source: "HighConfidence",
				})
				break
			}
		}
	}
	return hc
}

// ---------------------------------------------------------------------------
// BRM block detection
// ---------------------------------------------------------------------------

// calculateBRMBlocks applies BRM's block-threshold idea to smoothed windows.
// An AFP floor of afpFloor prevents the variance threshold from approaching
// zero near chromosomal fixation points, which would otherwise cause
// hypersensitive spurious blocks.
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
		// Clamp to [afpFloor, 1-afpFloor] to prevent near-zero variance at fixation.
		if afp < afpFloor {
			afp = afpFloor
		}
		if afp > 1-afpFloor {
			afp = 1 - afpFloor
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

// ---------------------------------------------------------------------------
// Shared chart style constants and helpers
// ---------------------------------------------------------------------------

const posFormatter = `function(value) {
	if (value >= 1000000) { return (value / 1000000).toFixed(2) + ' Mb'; }
	if (value >= 1000)    { return (value / 1000).toFixed(1)    + ' kb'; }
	return value;
}`

func commonGlobalOpts(title, subtitle, yLabel, width, height string, bidirectional bool) []charts.GlobalOpts {
	yMin := opts.Float(0.0)
	if bidirectional {
		yMin = nil
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
// Individual (raw-value) line chart — one series per statistic, per chromosome
// ---------------------------------------------------------------------------

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
// Composite best-signal chart — max |Z| across all 7 statistics per window.
// A single line that peaks wherever any statistic fires strongly, with z=2/z=3
// reference lines and BRM block shading identical to the Z-score overlay.
// ---------------------------------------------------------------------------

func createCompositeSignalChart(chrom string, x []int64, hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ []float64, brmBlocks []BRMBlock) *charts.Line {
	n := len(x)
	composite := make([]float64, n)
	for i := range composite {
		composite[i] = math.Max(math.Abs(hiZ[i]),
			math.Max(math.Abs(liZ[i]),
				math.Max(math.Abs(dsiZ[i]),
					math.Max(math.Abs(gsZ[i]),
						math.Max(math.Abs(edZ[i]),
							math.Max(math.Abs(lodZ[i]), math.Abs(bblZ[i])))))))
	}

	title := chrom + " — Composite Signal (max |Z|)"
	subtitle := "Max absolute robust Z-score across all 7 statistics per window. " +
		"z = 2 suggestive · z = 3 significant. Shaded bands: BRM blocks."

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "max |Z-score|", chartWidth, chartHeight, false)...)

	line.SetGlobalOptions(
		charts.WithTooltipOpts(opts.Tooltip{
			Show:        opts.Bool(true),
			Trigger:     "axis",
			AxisPointer: &opts.AxisPointer{Type: "cross"},
			Formatter: opts.FuncOpts(`function(params) {
				let pos = params[0].axisValue;
				let posStr = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
				let result = '<strong>' + posStr + '</strong><br/>';
				params.forEach(function(item) {
					let val = parseFloat(item.value);
					if (isNaN(val)) return;
					let sig = '';
					if (item.seriesName === 'Composite') {
						if (val >= 3.0)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ z≥3 significant</span>';
						else if (val >= 2.0) sig = ' <span style="color:#f39c12">● z≥2 suggestive</span>';
					}
					result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
				});
				return result;
			}`),
		}),
	)

	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}

	compositeData := make([]opts.LineData, n)
	for i, v := range composite {
		compositeData[i] = opts.LineData{Value: v}
	}

	compositeOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#2ca02c"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	compositeOpts = append(compositeOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	line.SetXAxis(positionLabels(x)).
		AddSeries("Composite", compositeData, compositeOpts...).
		AddSeries("z=2 (sugg.)", mkRef(zSugg),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=3 (sig.)", mkRef(zSig),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	return line
}

// ---------------------------------------------------------------------------
// Robust Z-score overlay chart
// ---------------------------------------------------------------------------

func statColor(name string) string {
	switch name {
	case "HighSI":
		return "#1f77b4"
	case "LowSI":
		return "#ff7f0e"
	case "DeltaSI":
		return "#2ca02c"
	case "Gstat":
		return "#17becf"
	case "ED":
		return "#d62728"
	case "LOD":
		return "#9467bd"
	case "BBLogBF":
		return "#8c564b"
	}
	return "#7f7f7f"
}

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

	type seriesDef struct {
		name  string
		data  []float64
		width float32
	}
	series := []seriesDef{
		{"HighSI", hiZ, 2.0},
		{"LowSI", liZ, 2.0},
		{"DeltaSI", dsiZ, 3.0},
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

// ---------------------------------------------------------------------------
// BRM mark-area helpers
// ---------------------------------------------------------------------------

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
// Utility helpers
// ---------------------------------------------------------------------------

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
	if err := f.Close(); err != nil {
		return fmt.Errorf("close %s: %w", path, err)
	}
	return nil
}

// ---------------------------------------------------------------------------
// Main plot + QTL entry point
// ---------------------------------------------------------------------------

// GenerateHtmlPlotsAndQTL processes smoothed window statistics, detects QTLs
// via four complementary methods, and writes three HTML files plus QTL TSVs.
//
// Output files:
//
//	GoBSAseq_IndividualPlots.html  – raw-value charts with permutation thresholds (7 per chromosome)
//	GoBSAseq_RobustZScore.html     – genome-wide robust Z-score overlay (all stats overlaid)
//	GoBSAseq_CompositeSignal.html  – single max-|Z| composite signal per chromosome
//	<qtlOutFile>                   – TSV of all QTL intervals (all detection methods)
//	GoBSAseq_BRMBlocks.tsv        – TSV of BRM-style block intervals used for plot shading
func GenerateHtmlPlotsAndQTL(allSmoothed []SmoothedStats, highSmAF, lowSmAF float64, highBulkSize, lowBulkSize int,
	population string, alphas []float64, rep int, htmlOutFile, qtlOutFile string) error {

	outDir := filepath.Dir(htmlOutFile)

	// -----------------------------------------------------------------------
	// Pass 1 — compute genome-wide robust Z normalisation parameters.
	// Each statistic's median and MAD are derived from all windows across all
	// chromosomes, with the top 1% trimmed so true QTL peaks don't inflate
	// the spread estimate.
	// -----------------------------------------------------------------------
	const trim = 0.01
	collectStat := func(fn func(SmoothedStats) float64) []float64 {
		v := make([]float64, len(allSmoothed))
		for i, s := range allSmoothed {
			v[i] = fn(s)
		}
		return v
	}
	hiMed, hiMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.HighSI }), trim)
	liMed, liMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.LowSI }), trim)
	dsiMed, dsiMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.DeltaSI }), trim)
	gsMed, gsMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.Gstat }), trim)
	edMed, edMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.ED }), trim)
	lodMed, lodMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.LOD }), trim)
	bblMed, bblMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.BBLogBF }), trim)

	// -----------------------------------------------------------------------
	// Pass 2 — group windows by chromosome and sort for consistent output.
	// -----------------------------------------------------------------------
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

	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics"

	robustZPage := components.NewPage()
	robustZPage.SetLayout(components.PageFlexLayout)
	robustZPage.PageTitle = "GoBSAseq — Robust Z-score Overlay"

	compositePage := components.NewPage()
	compositePage.SetLayout(components.PageFlexLayout)
	compositePage.PageTitle = "GoBSAseq — Composite Signal"

	for _, chrom := range chroms {
		stats := byChr[chrom]
		if len(stats) == 0 {
			continue
		}

		// Average permutation thresholds across windows — used to draw the
		// reference lines on the individual raw-value charts only.
		nf := float64(len(stats))
		var (
			sumHp99, sumHp95     float64
			sumHMp99, sumHMp95   float64
			sumLp99, sumLp95     float64
			sumLMp99, sumLMp95   float64
			sumDp99, sumDp95     float64
			sumDMp99, sumDMp95   float64
			sumGs99, sumGs95     float64
			sumEp99, sumEp95     float64
			sumLodp99, sumLodp95 float64
			sumBbp99, sumBbp95   float64
		)
		for _, s := range stats {
			t := s.thresholds
			sumHp99 += t.HighP99; sumHp95 += t.HighP95
			sumHMp99 += t.HighMp99; sumHMp95 += t.HighMp95
			sumLp99 += t.LowP99;  sumLp95 += t.LowP95
			sumLMp99 += t.LowMp99; sumLMp95 += t.LowMp95
			sumDp99 += t.DsiP99;  sumDp95 += t.DsiP95
			sumDMp99 += t.DsiMp99; sumDMp95 += t.DsiMp95
			sumGs99 += t.GsP99;  sumGs95 += t.GsP95
			sumEp99 += t.EdP99;  sumEp95 += t.EdP95
			sumLodp99 += t.LodP99; sumLodp95 += t.LodP95
			sumBbp99 += t.BbP99;  sumBbp95 += t.BbP95
		}
		avgHp99, avgHp95     := sumHp99/nf, sumHp95/nf
		avgHMp99, avgHMp95   := sumHMp99/nf, sumHMp95/nf
		avgLp99, avgLp95     := sumLp99/nf, sumLp95/nf
		avgLMp99, avgLMp95   := sumLMp99/nf, sumLMp95/nf
		avgDp99, avgDp95     := sumDp99/nf, sumDp95/nf
		avgDMp99, avgDMp95   := sumDMp99/nf, sumDMp95/nf
		avgGs99, avgGs95     := sumGs99/nf, sumGs95/nf
		avgEp99, avgEp95     := sumEp99/nf, sumEp95/nf
		avgLodp99, avgLodp95 := sumLodp99/nf, sumLodp95/nf
		avgBbp99, avgBbp95   := sumBbp99/nf, sumBbp95/nf

		// Extract per-window data and per-window p99/p95 threshold arrays.
		n := len(stats)
		x := make([]int64, n)
		hi := make([]float64, n)
		li := make([]float64, n)
		dsi := make([]float64, n)
		gs := make([]float64, n)
		ed := make([]float64, n)
		lod := make([]float64, n)
		bbl := make([]float64, n)
		hiT99, hiTM99 := make([]float64, n), make([]float64, n)
		hiT95, hiTM95 := make([]float64, n), make([]float64, n)
		liT99, liTM99 := make([]float64, n), make([]float64, n)
		liT95, liTM95 := make([]float64, n), make([]float64, n)
		dsiT99, dsiTM99 := make([]float64, n), make([]float64, n)
		dsiT95, dsiTM95 := make([]float64, n), make([]float64, n)
		gsT99, gsT95 := make([]float64, n), make([]float64, n)
		edT99, edT95 := make([]float64, n), make([]float64, n)
		lodT99, lodT95 := make([]float64, n), make([]float64, n)
		bblT99, bblT95 := make([]float64, n), make([]float64, n)

		for i, s := range stats {
			x[i] = s.POS
			hi[i], li[i], dsi[i] = s.HighSI, s.LowSI, s.DeltaSI
			gs[i], ed[i], lod[i], bbl[i] = s.Gstat, s.ED, s.LOD, s.BBLogBF
			t := s.thresholds
			hiT99[i], hiTM99[i] = t.HighP99, t.HighMp99
			hiT95[i], hiTM95[i] = t.HighP95, t.HighMp95
			liT99[i], liTM99[i] = t.LowP99, t.LowMp99
			liT95[i], liTM95[i] = t.LowP95, t.LowMp95
			dsiT99[i], dsiTM99[i] = t.DsiP99, t.DsiMp99
			dsiT95[i], dsiTM95[i] = t.DsiP95, t.DsiMp95
			gsT99[i], gsT95[i] = t.GsP99, t.GsP95
			edT99[i], edT95[i] = t.EdP99, t.EdP95
			lodT99[i], lodT95[i] = t.LodP99, t.LodP95
			bblT99[i], bblT95[i] = t.BbP99, t.BbP95
		}

		// ----------------------------------------------------------------
		// QTL detection — Method 1: locally-adaptive permutation thresholds.
		// Each window is compared to its own depth-specific p99/p95 threshold.
		// ----------------------------------------------------------------
		var chromQTLs []QTLRecord
		// p99
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, hi, hiT99, "HighSI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, hi, hiTM99, "HighSI", "99", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, li, liT99, "LowSI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, li, liTM99, "LowSI", "99", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, dsi, dsiT99, "DeltaSI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, dsi, dsiTM99, "DeltaSI", "99", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, gs, gsT99, "Gstat", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, ed, edT99, "ED4", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, lod, lodT99, "LOD", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, bbl, bblT99, "BBLogBF", "99", false, "Permutation")...)
		// p95
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, hi, hiT95, "HighSI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, hi, hiTM95, "HighSI", "95", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, li, liT95, "LowSI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, li, liTM95, "LowSI", "95", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, dsi, dsiT95, "DeltaSI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, dsi, dsiTM95, "DeltaSI", "95", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, gs, gsT95, "Gstat", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, ed, edT95, "ED4", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, lod, lodT95, "LOD", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, bbl, bblT95, "BBLogBF", "95", false, "Permutation")...)

		// ----------------------------------------------------------------
		// QTL detection — Method 2: genome-wide robust Z-score at z=±3.
		// Uses Pass 1 median/MAD so signals are normalised across the genome.
		// ----------------------------------------------------------------
		hiZ := robustZScore(hi, hiMed, hiMAD)
		liZ := robustZScore(li, liMed, liMAD)
		dsiZ := robustZScore(dsi, dsiMed, dsiMAD)
		gsZ := robustZScore(gs, gsMed, gsMAD)
		edZ := robustZScore(ed, edMed, edMAD)
		lodZ := robustZScore(lod, lodMed, lodMAD)
		bblZ := robustZScore(bbl, bblMed, bblMAD)

		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, hiZ, zSig, "HighSI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, liZ, zSig, "LowSI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, dsiZ, zSig, "DeltaSI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, dsiZ, -zSig, "DeltaSI_Z", "z3", true, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, gsZ, zSig, "Gstat_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, edZ, zSig, "ED4_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, lodZ, zSig, "LOD_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, bblZ, zSig, "BBLogBF_Z", "z3", false, "ZScore")...)

		// ----------------------------------------------------------------
		// QTL detection — Method 3: multi-statistic consensus.
		// Windows where ≥consensusMinStats statistics exceed their p99
		// thresholds simultaneously are merged into consensus intervals.
		// ----------------------------------------------------------------
		chromQTLs = append(chromQTLs, detectConsensusQTLs(chrom, stats)...)

		// ----------------------------------------------------------------
		// BRM blocks and Method 4: high-confidence intersection.
		// Permutation QTLs overlapping a BRM block are promoted to
		// HighConfidence, providing cross-method validation.
		// ----------------------------------------------------------------
		chromBRMBlocks := calculateBRMBlocks(chrom, stats, highBulkSize, lowBulkSize, popLevel, brmUAlpha)
		allBRMBlocks = append(allBRMBlocks, chromBRMBlocks...)
		chromQTLs = append(chromQTLs, intersectQTLsWithBRM(chromQTLs, chromBRMBlocks)...)
		allQTLs = append(allQTLs, chromQTLs...)

		robustZPage.AddCharts(createRobustZOverlayChart(chrom, x, hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))

		individualPage.AddCharts(
			createInteractiveLineChart(chrom+" HighSI", x, hi, avgHp99, avgHp95, avgHMp99, avgHMp95, true, chromBRMBlocks),
			createInteractiveLineChart(chrom+" LowSI", x, li, avgLp99, avgLp95, avgLMp99, avgLMp95, true, chromBRMBlocks),
			createInteractiveLineChart(chrom+" DeltaSI", x, dsi, avgDp99, avgDp95, avgDMp99, avgDMp95, true, chromBRMBlocks),
			createInteractiveLineChart(chrom+" Gstat", x, gs, avgGs99, avgGs95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" ED4", x, ed, avgEp99, avgEp95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" LOD", x, lod, avgLodp99, avgLodp95, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" BBLogBF", x, bbl, avgBbp99, avgBbp95, 0, 0, false, chromBRMBlocks),
		)

		compositePage.AddCharts(createCompositeSignalChart(chrom, x, hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))
	}

	// Write HTML files.
	if err := writeHTMLPage(individualPage, filepath.Join(outDir, "GoBSAseq_IndividualPlots.html")); err != nil {
		return err
	}
	if err := writeHTMLPage(robustZPage, filepath.Join(outDir, "GoBSAseq_RobustZScore.html")); err != nil {
		return err
	}
	if err := writeHTMLPage(compositePage, filepath.Join(outDir, "GoBSAseq_CompositeSignal.html")); err != nil {
		return err
	}

	// Write QTL TSV.
	fTsv, err := os.Create(qtlOutFile)
	if err != nil {
		return fmt.Errorf("create qtl file: %w", err)
	}
	fmt.Fprintf(fTsv, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\tSOURCE\n")
	for _, q := range allQTLs {
		fmt.Fprintf(fTsv, "%s\t%d\t%d\t%.6f\t%s\t%s\t%s\n",
			q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI, q.Source)
	}
	if err := fTsv.Close(); err != nil {
		return fmt.Errorf("close qtl file: %w", err)
	}

	// Write BRM blocks TSV.
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

// ---------------------------------------------------------------------------
// detectQTLsAdaptive — per-window threshold variant of detectQTLs
// ---------------------------------------------------------------------------

// detectQTLsAdaptive calls QTLs using a per-window threshold slice rather than
// a single chromosome-average threshold, making detection locally adaptive to
// sequencing depth.  thresholds[i] is the significance threshold for window x[i].
// When isValley is true, a window is considered significant if y[i] < thresholds[i].
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
			} else {
				if (!isValley && val > peak) || (isValley && val < peak) {
					peak = val
				}
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

	// Gap-bridging.
	merged := []run{runs[0]}
	for i := 1; i < len(runs); i++ {
		prev := &merged[len(merged)-1]
		prevStopIdx := sort.Search(len(x), func(j int) bool { return x[j] >= prev.stop })
		nextStartIdx := sort.Search(len(x), func(j int) bool { return x[j] >= runs[i].start })
		gap := nextStartIdx - prevStopIdx - 1
		if gap <= maxGapWindows {
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
		qtls = append(qtls, QTLRecord{
			Chrom:  chrom,
			Start:  r.start,
			Stop:   r.stop,
			Peak:   r.peak,
			Stat:   statName,
			CI:     ci,
			Source: source,
		})
	}
	return qtls
}

func RunTwoBulkTwoParents(cfg utils.AnalysisConfig, hfCfg utils.HardFilterConfig) {
	highParIdx := cfg.HighParentIdx
	//highParDP := cfg.HighParentDepth
	lowParIdx := cfg.LowParentIdx
	//lowParDP := cfg.LowParentDepth
	highBulkIdx := cfg.HighBulkIdx
	//highBulkDP := cfg.HighBulkDepth
	lowBulkIdx := cfg.LowBulkIdx
	//lowBulkDP := cfg.LowBulkDepth
	vcfRdr := cfg.Rdr
	outDir := cfg.OutputDir

	windowSize := int64(cfg.WindowSize)
	stepSize := int64(cfg.StepSize)
	rep := cfg.Rep
	pop := cfg.Population

	highSmAF := utils.SimulateAF(pop, float64(cfg.HighBulkSize), cfg.Rep)
	lowSmAF := utils.SimulateAF(pop, float64(cfg.LowBulkSize), cfg.Rep)

	overallStart := time.Now()

	// -----------------------------------------------------------------------
	// Stage 0 — hard filtering
	// -----------------------------------------------------------------------
	color.Cyan("============================ GATK Hard Filtering ============================\n\n")

	filteredVcfPath := filepath.Join(outDir, "GoBSAseq.hard_filtered.vcf.gz")
	badVcfPath := filepath.Join(outDir, "GoBSAseq.bad_variants.vcf.gz")

	passedVariants, original, hardFiltered, err := utils.HardFilterVcf(vcfRdr, filteredVcfPath, badVcfPath, cfg, hfCfg)
	if err != nil {
		color.Red("Hard filter error: %v", err)
		return
	}
	color.Green("Original variants: %v\nHard filtered variants: %v", original, hardFiltered)
	// -----------------------------------------------------------------------
	// Stage 1 — per-SNP statistics (concurrent)
	// -----------------------------------------------------------------------
	color.Cyan("============================ Calculating Statistics =============================\n\n")

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

	bar := progressbar.Default(int64(len(passedVariants)), "Processing variants")

	// Feed worker pool from the pre-filtered slice.
	variantChan := make(chan *vcfgo.Variant, 10000)
	go func() {
		for _, v := range passedVariants {
			variantChan <- v
		}
		close(variantChan)
	}()

	for i := 0; i < numWorkers; i++ {
		workerWG.Add(1)
		go func() {
			defer workerWG.Done()
			for variant := range variantChan {
				// Basic safety: Ensure samples exist for the configured indices
				maxIdx := highParIdx
				for _, idx := range []int{lowParIdx, highBulkIdx, lowBulkIdx} {
					if idx > maxIdx {
						maxIdx = idx
					}
				}
				if len(variant.Samples) <= maxIdx {
					_ = bar.Add(1)
					continue
				}

				hpS := variant.Samples[highParIdx]
				lpS := variant.Samples[lowParIdx]
				hbS := variant.Samples[highBulkIdx]
				lbS := variant.Samples[lowBulkIdx]

				// Identify the one 'real' ALT allele index.
				realAltIdx := -1
				for i, alt := range variant.Alt() {
					if !utils.IsSymbolic(alt) {
						realAltIdx = i
						break
					}
				}
				if realAltIdx == -1 {
					_ = bar.Add(1)
					continue
				}

				// Identify which allele is from the High Parent.
				// passesTwoBulkBSAseqFilters guarantees parents are homozygous
				// and carry either 0 (REF) or the one real ALT allele.
				hpAllele := hpS.GT[0]

				// Get ADs for each bulk.
				hbRefDep, _ := hbS.RefDepth()
				hbAltDeps, _ := hbS.AltDepths()
				lbRefDep, _ := lbS.RefDepth()
				lbAltDeps, _ := lbS.AltDepths()

				// Safety check: ensure the target ALT allele depth exists.
				if len(hbAltDeps) <= realAltIdx || len(lbAltDeps) <= realAltIdx {
					_ = bar.Add(1)
					continue
				}

				var hbH, hbL, lbH, lbL int
				// SI (Success Index) is defined as the frequency of the High Parent allele.
				if hpAllele == 0 {
					// High Parent allele is REF (0), Low Parent is ALT.
					hbH, hbL = hbRefDep, hbAltDeps[realAltIdx]
					lbH, lbL = lbRefDep, lbAltDeps[realAltIdx]
				} else {
					// High Parent allele is ALT, Low Parent is REF (0).
					hbH, hbL = hbAltDeps[realAltIdx], hbRefDep
					lbH, lbL = lbAltDeps[realAltIdx], lbRefDep
				}

				// Total depth for SI calculation is the sum of relevant alleles.
				// This is more robust than s.DP if there are symbolic/other alleles.
				hbTotal := hbH + hbL
				lbTotal := lbH + lbL

				if hbTotal == 0 || lbTotal == 0 {
					_ = bar.Add(1)
					continue
				}

				hSI := float64(hbH) / float64(hbTotal)
				lSI := float64(lbH) / float64(lbTotal)
				minDepth := hbTotal
				if lbTotal < minDepth {
					minDepth = lbTotal
				}

				s := BSAstats{
					CHROM:      variant.Chromosome,
					POS:        int64(variant.Pos),
					REF:        variant.Reference,
					ALT:        variant.Alt()[realAltIdx],
					HighParGT:  hpS.GT,
					LowParGT:   lpS.GT,
					HighBulkGT: hbS.GT,
					HighBulkAD: fmt.Sprintf("%d,%d", hbRefDep, hbAltDeps[realAltIdx]),
					LowBulkGT:  lbS.GT,
					LowBulkAD:  fmt.Sprintf("%d,%d", lbRefDep, lbAltDeps[realAltIdx]),
					HighBulkL:  hbL,
					HighBulkH:  hbH,
					LowBulkL:   lbL,
					LowBulkH:   lbH,
					HighSI:     math.Round(hSI*1e6) / 1e6,
					LowSI:      math.Round(lSI*1e6) / 1e6,
					DeltaSI:    math.Round((hSI-lSI)*1e6) / 1e6,
					Gstat:      math.Round(GStatistic(hbH, hbL, lbH, lbL)*1e6) / 1e6,
					ED:         math.Round(euclideanDistance4(hSI, lSI)*1e6) / 1e6,
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

	// -----------------------------------------------------------------------
	// Stage 2 — smoothing
	// -----------------------------------------------------------------------
	color.Cyan("\n============================ Smoothing Statistics =============================\n\n")

	var allSmoothed []SmoothedStats
	for chrom, stats := range chromStats {
		color.Yellow("Smoothing %s: %d SNPs", chrom, len(stats))
		smoothed := smoothChromosome(stats, windowSize, stepSize)
		allSmoothed = append(allSmoothed, smoothed...)
	}

	// -----------------------------------------------------------------------
	// Stage 3 — threshold simulation
	// -----------------------------------------------------------------------
	color.Cyan("\n============================ Calculating Thresholds (%d simulations per depth pair) ==============================\n\n", rep)
	calcAllThresholds(allSmoothed, highSmAF, lowSmAF, rep)
	// Attach per-window thresholds to the SmoothedStats slice so downstream
	// QTL detection can use locally adaptive thresholds rather than chromosome averages.
	for i := range allSmoothed {
		allSmoothed[i].thresholds = calcThresholdsCached(
			allSmoothed[i].MeanHighBulkDP, allSmoothed[i].MeanLowBulkDP, highSmAF, lowSmAF, rep)
	}
	color.Green("\nThreshold calculations complete.")

	// -----------------------------------------------------------------------
	// Stage 4 — smoothed TSV
	// -----------------------------------------------------------------------
	color.Cyan("\n=========================================== Writing Smoothed TSV =================================================\n\n")
	smoothTSV := filepath.Join(outDir, "GoBSAseq.smooth.tsv")
	if err := writeSmoothedTSV(smoothTSV, allSmoothed, highSmAF, lowSmAF, rep); err != nil {
		color.Red("Error writing smoothed TSV: %v", err)
	} else {
		color.Green("Wrote %d smoothed windows to %s", len(allSmoothed), smoothTSV)
	}
	color.Green("Raw stats written to %s", filepath.Join(outDir, "GoBSAseq.raw.tsv"))
	color.Green("\nTotal time: %s\n", time.Since(overallStart).Round(time.Second))

	// -----------------------------------------------------------------------
	// Stage 5 — plots and QTL detection
	// -----------------------------------------------------------------------
	color.Cyan("\n============================ Generating HTML Plots & QTLs ========================================\n\n")
	htmlFile := filepath.Join(outDir, "GoBSAseq_RobustZScore.html")
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
		color.Green("HTML plots written to %s", outDir)
		color.Green("QTL tabular results written to %s", qtlFile)
	}
}

// ---------------------------------------------------------------------------
// Threshold simulation and caching
// ---------------------------------------------------------------------------
