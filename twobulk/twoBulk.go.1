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
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
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
	ED      float64
	LOD     float64
	BBLogBF float64

	Depth int
}

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

var thresholdCache sync.Map

func calcThresholds(highBulkDP, lowBulkDP int, highSmAF, lowSmAF float64, rep int) Thresholds {
	src := rand.NewSource(time.Now().UnixNano())
	rng := rand.New(src)

	// gonum distuv.Binomial needs its own source
	distSrc := distuv.Binomial{N: float64(highBulkDP), P: highSmAF, Src: rng}
	distLow := distuv.Binomial{N: float64(lowBulkDP), P: lowSmAF, Src: rng}

	highSIArr := make([]float64, rep)
	lowSIArr := make([]float64, rep)
	dsiArr := make([]float64, rep)
	gsArr := make([]float64, rep)
	edArr := make([]float64, rep)
	lodArr := make([]float64, rep)
	bbArr := make([]float64, rep)

	for i := 0; i < rep; i++ {
		hAlt := distSrc.Rand()
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
	key := fmt.Sprintf("%d_%d_%.6f_%.6f", highBulkDP, lowBulkDP, highSmAF, lowSmAF)
	if v, ok := thresholdCache.Load(key); ok {
		return v.(Thresholds)
	}
	t := calcThresholds(highBulkDP, lowBulkDP, highSmAF, lowSmAF, rep)
	thresholdCache.Store(key, t)
	return t
}

func calcAllThresholds(allSmoothed []SmoothedStats, highSmAF, lowSmAF float64, rep int) {
	// Collect unique depth pairs
	type depthPair struct{ h, l int }
	seen := make(map[depthPair]bool)
	for _, sm := range allSmoothed {
		dp := depthPair{sm.MeanHighBulkDP, sm.MeanLowBulkDP}
		if sm.MeanHighBulkDP > 0 && sm.MeanLowBulkDP > 0 {
			seen[dp] = true
		}
	}

	pairs := make([]depthPair, 0, len(seen))
	for dp := range seen {
		pairs = append(pairs, dp)
	}

	total := len(pairs)
	if total == 0 {
		return
	}

	bar := progressbar.NewOptions(total,
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

	pairChan := make(chan depthPair, total)
	for _, dp := range pairs {
		pairChan <- dp
	}
	close(pairChan)

	numWorkers := runtime.NumCPU()
	var wg sync.WaitGroup
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for dp := range pairChan {
				calcThresholdsCached(dp.h, dp.l, highSmAF, lowSmAF, rep)
				bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println() // newline after progress bar
}

func gStatisticFloat(hAlt, hRef, lAlt, lRef float64) float64 {
	n1 := hAlt + hRef
	n2 := lAlt + lRef
	if n1 == 0 || n2 == 0 {
		return 0
	}
	totalAlt := hAlt + lAlt
	total := n1 + n2
	p := totalAlt / total

	g := 0.0
	if hAlt > 0 {
		g += hAlt * math.Log(hAlt/(n1*p))
	}
	if hRef > 0 {
		g += hRef * math.Log(hRef/(n1*(1-p)))
	}
	if lAlt > 0 {
		g += lAlt * math.Log(lAlt/(n2*p))
	}
	if lRef > 0 {
		g += lRef * math.Log(lRef/(n2*(1-p)))
	}
	return 2 * g
}

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
		return 0.0
	}

	eps := 1e-10
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

func GoodVariants(v *vcfgo.Variant, highPar, highParDP, lowPar, lowParDP, highBulk, highBulkDP, lowBulk, lowBulkDP int) bool {
	indices := []int{highPar, lowPar, highBulk, lowBulk}
	if len(v.Alt()) != 1 {
		return false
	}

	for _, idx := range indices {
		if idx >= len(v.Samples) {
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
	if hpDP < highParDP || lpDP < lowParDP || hbDP < highBulkDP || lbDP < lowBulkDP {
		return false
	}
	return true
}

func tricubeWeight(d, D float64) float64 {
	if d >= D {
		return 0
	}
	x := 1 - math.Pow(d/D, 3)
	return x * x * x
}

func smoothChromosome(stats []BSAstats, windowSize int64, step int64) []SmoothedStats {
	if len(stats) == 0 {
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
			wDeltaSI := w * depthWeight

			sumDeltaSI += s.DeltaSI * wDeltaSI
			sumWeightDeltaSI += wDeltaSI

			sumED += s.ED * wDeltaSI
			sumWeightED += wDeltaSI

			sumGstat += s.Gstat * wDeltaSI
			sumLOD += s.LOD * wDeltaSI
			sumBBLogBF += s.BBLogBF * wDeltaSI

			sumHighSI += s.HighSI * wDeltaSI
			sumLowSI += s.LowSI * wDeltaSI
			sumWeightSI += wDeltaSI

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
					bar.Add(1)
					continue
				}

				lpGT := variant.Samples[lowParIdx].GT
				hbGT := variant.Samples[highBulkIdx].GT
				lbGT := variant.Samples[lowBulkIdx].GT
				hbDP, lbDP := variant.Samples[highBulkIdx].DP, variant.Samples[lowBulkIdx].DP

				if hbDP == 0 || lbDP == 0 {
					bar.Add(1)
					continue
				}

				hbRefDep, _ := variant.Samples[highBulkIdx].RefDepth()
				hbAltDeps, _ := variant.Samples[highBulkIdx].AltDepths()
				lbRefDep, _ := variant.Samples[lowBulkIdx].RefDepth()
				lbAltDeps, _ := variant.Samples[lowBulkIdx].AltDepths()

				if len(hbAltDeps) == 0 || len(lbAltDeps) == 0 {
					bar.Add(1)
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

				//var hbL, hbH, lbL, lbH int
				//lpAllele := 0
				//if len(lpGT) > 0 {
				//	lpAllele = lpGT[0]
				//}
				//hbAllele := 0
				//if len(hbGT) > 0 {
				//	hbAllele = hbGT[0]
				//}
				//lbAllele := 0
				//if len(lbGT) > 0 {
				//	lbAllele = lbGT[0]
				//}
				//
				//if hbAllele == lpAllele {
				//	hbL, hbH = hbRefDep, hbAltDeps[0]
				//} else {
				//	hbL, hbH = hbAltDeps[0], hbRefDep
				//}
				//if lbAllele == lpAllele {
				//	lbL, lbH = lbRefDep, lbAltDeps[0]
				//} else {
				//	lbL, lbH = lbAltDeps[0], lbRefDep
				//}

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
					HighSI:     math.Round((hSI)*1e6) / 1e6,
					LowSI:      math.Round((lSI)*1e6) / 1e6,
					DeltaSI:    math.Round((hSI-lSI)*1e6) / 1e6,
					Gstat:      math.Round((GStatistic(hbH, hbL, lbH, lbL))*1e6) / 1e6,
					ED:         math.Round((math.Abs(hSI-lSI))*1e6) / 1e6,
					LOD:        math.Round((lod(hbL, hbH, lbL, lbH))*1e6) / 1e6,
					BBLogBF:    math.Round((betaBinomialLogBF(hbH, hbL, lbH, lbL))*1e6) / 1e6,
					Depth:      minDepth,
				}
				statsChan <- s
				rawWriteChan <- s
				bar.Add(1)
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
	bar.Finish()
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
	if err := writeSmoothedTSV(filepath.Join(outDir, "GoBSAseq.smooth.tsv"), allSmoothed, highSmAF, lowSmAF, rep); err != nil {
		color.Red("Error writing smoothed TSV: %v", err)
	} else {
		color.Green("Wrote %d smoothed windows to %s", len(allSmoothed), filepath.Join(outDir, "GoBSAseq.smooth.tsv"))
	}
	color.Green("Raw stats written to %s", filepath.Join(outDir, "GoBSAseq.raw.tsv"))
	color.Green("\nTotal time: %s\n", time.Since(overallStart).Round(time.Second))

	color.Cyan("\n============================ Generating HTML Plots & QTLs ========================================\n\n")
	htmlFile := filepath.Join(outDir, "GoBSAseq_InteractivePlots.html")
	qtlFile := filepath.Join(outDir, "GoBSAseq_QTL.tsv")

	if err := GenerateHtmlPlotsAndQTL(allSmoothed, highSmAF, lowSmAF, rep, htmlFile, qtlFile); err != nil {
		color.Red("Error generating Plots and QTLs: %v", err)
	} else {
		color.Green("Interactive HTML plots written to %s", htmlFile)
		color.Green("QTL tabular results written to %s", qtlFile)
	}
}
