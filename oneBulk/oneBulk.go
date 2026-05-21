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

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/genespace"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

const (
	minSNPsPerWindow  = 5
	maxGapWindows     = 3
	consensusMinStats = 3
	afpFloor          = 0.05
	defaultBRMAlpha   = 0.05
)

// OneBulkStats holds per-SNP statistics for single-bulk BSA-seq.
type OneBulkStats struct {
	CHROM     string
	POS       int64
	REF       string
	ALT       string
	HighParGT []int
	LowParGT  []int
	BulkGT    []int
	BulkAD    string
	SI        float64
	AbsSI     float64
	Gstat     float64
	LOD       float64
	BBLogBF   float64
	ED        float64
	Depth     int
}

// SmoothedStats holds window-averaged statistics for one genomic window.
type SmoothedStats struct {
	CHROM      string
	POS        int64
	SI         float64
	AbsSI      float64
	Gstat      float64
	ED         float64
	LOD        float64
	BBLogBF    float64
	NumSNPs    int
	MeanBulkDP int
	thresholds Thresholds
}

// Thresholds holds permutation-derived significance levels for each statistic.
type Thresholds struct {
	SiP99   float64
	SiP95   float64
	SiMp99  float64
	SiMp95  float64
	AbsP99  float64
	AbsP95  float64
	GsP99   float64
	GsP95   float64
	EdP99   float64
	EdP95   float64
	LodP99  float64
	LodP95  float64
	BbP99   float64
	BbP95   float64
}

// QTLRecord holds a detected QTL interval and its peak value.
type QTLRecord struct {
	Chrom  string
	Start  int64
	Stop   int64
	Peak   float64
	Stat   string
	CI     string
	Source string
}

// BRMBlock holds one BRM-style interval.
type BRMBlock struct {
	Chrom     string
	Start     int64
	Stop      int64
	PeakPos   int64
	Peak      float64
	Threshold float64
}

var thresholdCache sync.Map

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
	alt, ref := gt[0], gt[1]
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
	logNull := total * math.Log(0.5)
	logAlt := logBeta(a+1, r+1) - logBeta(1, 1)
	return logAlt - logNull
}

func buildOneBulkStats(v *vcfgo.Variant, highParIdx, lowParIdx, bulkIdx int, trackHighParent bool) (OneBulkStats, bool) {
	realAltIdx := -1
	for i, alt := range v.Alt() {
		if !(alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<')) {
			realAltIdx = i
			break
		}
	}
	if realAltIdx == -1 {
		return OneBulkStats{}, false
	}

	maxIdx := bulkIdx
	for _, idx := range []int{highParIdx, lowParIdx} {
		if idx > maxIdx {
			maxIdx = idx
		}
	}
	if len(v.Samples) <= maxIdx {
		return OneBulkStats{}, false
	}

	hpS := v.Samples[highParIdx]
	lpS := v.Samples[lowParIdx]
	bS := v.Samples[bulkIdx]

	bRefDep, _ := bS.RefDepth()
	bAltDeps, _ := bS.AltDepths()
	if len(bAltDeps) <= realAltIdx {
		return OneBulkStats{}, false
	}
	bAltDep := bAltDeps[realAltIdx]
	bTotal := bRefDep + bAltDep
	if bTotal == 0 {
		return OneBulkStats{}, false
	}

	hpAllele := hpS.GT[0]
	var bulkSus, bulkRes int
	if trackHighParent {
		if hpAllele == 0 {
			bulkSus, bulkRes = bRefDep, bAltDep
		} else {
			bulkSus, bulkRes = bAltDep, bRefDep
		}
	} else {
		if lpS.GT[0] == 0 {
			bulkSus, bulkRes = bRefDep, bAltDep
		} else {
			bulkSus, bulkRes = bAltDep, bRefDep
		}
	}

	si := float64(bulkSus) / float64(bTotal)
	absSI := math.Abs(si - 0.5)
	bulkGT := []int{bulkSus, bulkRes}

	return OneBulkStats{
		CHROM:     v.Chromosome,
		POS:       int64(v.Pos),
		REF:       v.Reference,
		ALT:       v.Alt()[realAltIdx],
		HighParGT: hpS.GT,
		LowParGT:  lpS.GT,
		BulkGT:    bulkGT,
		BulkAD:    fmt.Sprintf("%d,%d", bulkSus, bulkRes),
		SI:        math.Round(si*1e6) / 1e6,
		AbsSI:     math.Round(absSI*1e6) / 1e6,
		ED:        math.Round(math.Pow(absSI, 4)*1e6) / 1e6,
		Gstat:     math.Round(GStatisticOneBulk(bulkGT)*1e6) / 1e6,
		LOD:       math.Round(LodOneBulk(bulkGT)*1e6) / 1e6,
		BBLogBF:   math.Round(BetaBinomialOneBulk(bulkGT)*1e6) / 1e6,
		Depth:     bTotal,
	}, true
}

func calcThresholdsOneBulk(bulkDP int, smAF float64, rep int) Thresholds {
	if bulkDP <= 0 || rep <= 0 || smAF <= 0 || smAF >= 1 {
		return Thresholds{}
	}

	src := rand.NewSource(time.Now().UnixNano())
	rng := rand.New(src)
	dist := distuv.Binomial{N: float64(bulkDP), P: smAF, Src: rng}

	siArr := make([]float64, rep)
	absArr := make([]float64, rep)
	gsArr := make([]float64, rep)
	edArr := make([]float64, rep)
	lodArr := make([]float64, rep)
	bbArr := make([]float64, rep)

	for i := 0; i < rep; i++ {
		alt := int(dist.Rand())
		ref := bulkDP - alt
		gt := []int{alt, ref}
		si := float64(alt) / float64(bulkDP)
		absSI := math.Abs(si - 0.5)

		siArr[i] = math.Round(si*1e6) / 1e6
		absArr[i] = math.Round(absSI*1e6) / 1e6
		gsArr[i] = math.Round(GStatisticOneBulk(gt)*1e6) / 1e6
		edArr[i] = math.Round(math.Pow(absSI, 4)*1e6) / 1e6
		lodArr[i] = math.Round(LodOneBulk(gt)*1e6) / 1e6
		bbArr[i] = math.Round(BetaBinomialOneBulk(gt)*1e6) / 1e6
	}

	sort.Float64s(siArr)
	sort.Float64s(absArr)
	sort.Float64s(gsArr)
	sort.Float64s(edArr)
	sort.Float64s(lodArr)
	sort.Float64s(bbArr)

	return Thresholds{
		SiP99:  roundQuantile(siArr, 0.995),
		SiP95:  roundQuantile(siArr, 0.95),
		SiMp99: roundQuantile(siArr, 0.005),
		SiMp95: roundQuantile(siArr, 0.05),
		AbsP99: roundQuantile(absArr, 0.995),
		AbsP95: roundQuantile(absArr, 0.95),
		GsP99:  roundQuantile(gsArr, 0.995),
		GsP95:  roundQuantile(gsArr, 0.95),
		EdP99:  roundQuantile(edArr, 0.995),
		EdP95:  roundQuantile(edArr, 0.95),
		LodP99: roundQuantile(lodArr, 0.995),
		LodP95: roundQuantile(lodArr, 0.95),
		BbP99:  roundQuantile(bbArr, 0.995),
		BbP95:  roundQuantile(bbArr, 0.95),
	}
}

func roundQuantile(sorted []float64, p float64) float64 {
	return math.Round(stat.Quantile(p, stat.Empirical, sorted, nil)*1e6) / 1e6
}

func calcThresholdsCachedOneBulk(bulkDP int, smAF float64, rep int) Thresholds {
	key := fmt.Sprintf("1b_%d_%.6f_%d", bulkDP, smAF, rep)
	if v, ok := thresholdCache.Load(key); ok {
		return v.(Thresholds)
	}
	t := calcThresholdsOneBulk(bulkDP, smAF, rep)
	actual, _ := thresholdCache.LoadOrStore(key, t)
	return actual.(Thresholds)
}

func calcAllThresholdsOneBulk(allSmoothed []SmoothedStats, smAF float64, rep int) {
	seen := make(map[int]bool)
	for _, sm := range allSmoothed {
		if sm.MeanBulkDP > 0 {
			seen[sm.MeanBulkDP] = true
		}
	}
	deps := make([]int, 0, len(seen))
	for d := range seen {
		deps = append(deps, d)
	}
	if len(deps) == 0 {
		return
	}

	bar := progressbar.NewOptions(len(deps),
		progressbar.OptionSetDescription("Computing single-bulk thresholds"),
		progressbar.OptionSetWidth(40),
		progressbar.OptionShowCount(),
	)
	depChan := make(chan int, len(deps))
	for _, d := range deps {
		depChan <- d
	}
	close(depChan)

	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for dp := range depChan {
				calcThresholdsCachedOneBulk(dp, smAF, rep)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println()
}

func tricubeWeight(d, D float64) float64 {
	if D <= 0 || d >= D {
		return 0
	}
	x := 1 - math.Pow(d/D, 3)
	return x * x * x
}

func smoothChromosomeOneBulk(stats []OneBulkStats, windowSize, step int64) []SmoothedStats {
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
			sumSI, sumAbs, sumWeightStat float64
			sumGstat, sumWeightGs        float64
			sumLOD, sumWeightLod         float64
			sumBBLogBF, sumWeightBB      float64
			sumED, sumWeightED           float64
			sumDP                        float64
			nSNPs                        int
		)

		for _, s := range stats {
			if s.POS < windowStart || s.POS > windowEnd {
				continue
			}
			nSNPs++
			d := math.Abs(float64(s.POS - center))
			w := tricubeWeight(d, float64(windowSize)/2)
			wStat := w * math.Sqrt(float64(s.Depth))

			sumSI += s.SI * wStat
			sumAbs += s.AbsSI * wStat
			sumWeightStat += wStat
			sumGstat += s.Gstat * wStat
			sumWeightGs += wStat
			sumLOD += s.LOD * wStat
			sumWeightLod += wStat
			sumBBLogBF += s.BBLogBF * wStat
			sumWeightBB += wStat
			sumED += s.ED * wStat
			sumWeightED += wStat
			sumDP += float64(s.Depth)
		}

		if nSNPs < minSNPsPerWindow {
			continue
		}

		sm := SmoothedStats{CHROM: chrom, POS: center, NumSNPs: nSNPs, MeanBulkDP: int(sumDP / float64(nSNPs))}
		if sumWeightStat > 0 {
			sm.SI = sumSI / sumWeightStat
			sm.AbsSI = sumAbs / sumWeightStat
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
		smoothed = append(smoothed, sm)
	}
	return smoothed
}

func writeRawTSVOneBulk(filename string, statsChan <-chan OneBulkStats) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	fmt.Fprintln(w, "CHROM\tPOS\tREF\tALT\tHighParGT\tLowParGT\tBulkGT\tBulkAD\tSI\tAbsSI\tGstat\tED4\tLOD\tBBLogBF\tDepth")
	for s := range statsChan {
		fmt.Fprintf(w, "%s\t%d\t%s\t%s\t%v\t%v\t%v\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
			s.CHROM, s.POS, s.REF, s.ALT,
			s.HighParGT, s.LowParGT, s.BulkGT, s.BulkAD,
			s.SI, s.AbsSI, s.Gstat, s.ED, s.LOD, s.BBLogBF, s.Depth)
	}
	return nil
}

func writeSmoothedTSVOneBulk(filename string, data []SmoothedStats, smAF float64, rep int) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	header := "CHROM\tPOS\tSI\tAbsSI\tGstat\tED4\tLOD\tBBLogBF\tNumSNPs\tMeanBulkDP" +
		"\tSI_p99\tSI_p95\tSI_m_p99\tSI_m_p95" +
		"\tAbsSI_p99\tAbsSI_p95" +
		"\tGstat_p99\tGstat_p95" +
		"\tED4_p99\tED4_p95" +
		"\tLOD_p99\tLOD_p95" +
		"\tBBLogBF_p99\tBBLogBF_p95"
	fmt.Fprintln(w, header)

	for _, d := range data {
		t := calcThresholdsCachedOneBulk(d.MeanBulkDP, smAF, rep)
		fmt.Fprintf(w, "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d"+
			"\t%.6f\t%.6f\t%.6f\t%.6f"+
			"\t%.6f\t%.6f"+
			"\t%.6f\t%.6f"+
			"\t%.6f\t%.6f"+
			"\t%.6f\t%.6f"+
			"\t%.6f\t%.6f\n",
			d.CHROM, d.POS, d.SI, d.AbsSI, d.Gstat, d.ED, d.LOD, d.BBLogBF, d.NumSNPs, d.MeanBulkDP,
			t.SiP99, t.SiP95, t.SiMp99, t.SiMp95,
			t.AbsP99, t.AbsP95,
			t.GsP99, t.GsP95,
			t.EdP99, t.EdP95,
			t.LodP99, t.LodP95,
			t.BbP99, t.BbP95)
	}
	return nil
}

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
		if s.SI > t.SiP99 || s.SI < t.SiMp99 {
			fired = append(fired, "SI")
		}
		if s.AbsSI > t.AbsP99 {
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

func calculateBRMBlocksOneBulk(chrom string, stats []SmoothedStats, bulkSize, popLevel int, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || bulkSize <= 0 {
		return nil
	}
	n := float64(bulkSize)
	popScale := math.Pow(2, float64(popLevel))
	varianceScale := 1.0 / (popScale * n)

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
			Chrom: chrom, Start: start, Stop: stop,
			PeakPos: stats[peakIdx].POS, Peak: peak, Threshold: threshold,
		})
	}

	for i, s := range stats {
		afp := s.SI
		if afp < afpFloor {
			afp = afpFloor
		}
		if afp > 1-afpFloor {
			afp = 1 - afpFloor
		}
		threshold := uAlpha * math.Sqrt(varianceScale*afp*(1-afp))
		deviation := math.Abs(s.SI - 0.5)
		significant := threshold > 0 && deviation >= threshold

		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = deviation
				peakThreshold = threshold
				continue
			}
			if deviation > peak {
				peakIdx = i
				peak = deviation
				peakThreshold = threshold
			}
			continue
		}
		if inBlock {
			emitBlock(startIdx, i-1, peakIdx, peak, peakThreshold)
			inBlock = false
		}
	}
	if inBlock {
		emitBlock(startIdx, len(stats)-1, peakIdx, peak, peakThreshold)
	}
	return blocks
}

// RunTwoParentsLowBulk runs BSA-seq with two parents and a low (susceptible) bulk.
func RunTwoParentsLowBulk(cfg utils.AnalysisConfig, hfCfg utils.HardFilterConfig) error {
	return runTwoParentsOneBulk(cfg, hfCfg, 2, false, cfg.LowBulkSize)
}

// RunTwoParentsHighBulk runs BSA-seq with two parents and a high (resistant) bulk.
func RunTwoParentsHighBulk(cfg utils.AnalysisConfig, hfCfg utils.HardFilterConfig) error {
	return runTwoParentsOneBulk(cfg, hfCfg, 3, true, cfg.HighBulkSize)
}

func runTwoParentsOneBulk(cfg utils.AnalysisConfig, hfCfg utils.HardFilterConfig, bsaType int, trackHighParent bool, bulkSize int) error {
	highParIdx := cfg.HighParentIdx
	lowParIdx := cfg.LowParentIdx
	bulkIdx := cfg.LowBulkIdx
	if trackHighParent {
		bulkIdx = cfg.HighBulkIdx
	}

	vcfRdr := cfg.Rdr
	outDir := cfg.OutputDir
	windowSize := int64(cfg.WindowSize)
	stepSize := int64(cfg.StepSize)
	rep := cfg.Rep
	pop := cfg.Population

	smAF := utils.SimulateAF(pop, float64(bulkSize), rep)
	overallStart := time.Now()

	label := "low bulk"
	if trackHighParent {
		label = "high bulk"
	}
	color.Cyan("============================ GATK Hard Filtering (Two parents, %s) ============================\n\n", label)

	filteredVcfPath := filepath.Join(outDir, "GoBSAseq.hard_filtered.vcf.gz")
	badVcfPath := filepath.Join(outDir, "GoBSAseq.bad_variants.vcf.gz")

	passedVariants, original, hardFiltered, err := utils.HardFilterVcf(vcfRdr, filteredVcfPath, badVcfPath, cfg, hfCfg, bsaType)
	if err != nil {
		return fmt.Errorf("hard filter: %w", err)
	}

	cfg.HighParentIdx = 0
	cfg.LowParentIdx = 1
	if trackHighParent {
		cfg.HighBulkIdx = 2
	} else {
		cfg.LowBulkIdx = 2
	}
	highParIdx = 0
	lowParIdx = 1
	bulkIdx = 2

	color.Green("Original variants: %v\nHard filtered variants: %v", original, hardFiltered)
	color.Cyan("============================ Calculating Single-Bulk Statistics =============================\n\n")

	statsChan := make(chan OneBulkStats, 10000)
	rawWriteChan := make(chan OneBulkStats, 10000)
	numWorkers := runtime.NumCPU()
	var workerWG sync.WaitGroup

	var rawWG sync.WaitGroup
	rawWG.Add(1)
	go func() {
		defer rawWG.Done()
		if err := writeRawTSVOneBulk(filepath.Join(outDir, "GoBSAseq.raw.tsv"), rawWriteChan); err != nil {
			color.Red("Error writing raw TSV: %v", err)
		}
	}()

	bar := progressbar.Default(int64(len(passedVariants)), "Processing variants")
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
				s, ok := buildOneBulkStats(variant, highParIdx, lowParIdx, bulkIdx, trackHighParent)
				if !ok {
					_ = bar.Add(1)
					continue
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

	chromStats := make(map[string][]OneBulkStats)
	for s := range statsChan {
		chromStats[s.CHROM] = append(chromStats[s.CHROM], s)
	}
	_ = bar.Finish()
	rawWG.Wait()

	color.Cyan("\n============================ Smoothing Statistics =============================\n\n")
	var allSmoothed []SmoothedStats
	for chrom, stats := range chromStats {
		color.Yellow("Smoothing %s: %d SNPs", chrom, len(stats))
		allSmoothed = append(allSmoothed, smoothChromosomeOneBulk(stats, windowSize, stepSize)...)
	}

	color.Cyan("\n============================ Calculating Thresholds (%d simulations per depth) ==============================\n\n", rep)
	calcAllThresholdsOneBulk(allSmoothed, smAF, rep)
	for i := range allSmoothed {
		allSmoothed[i].thresholds = calcThresholdsCachedOneBulk(allSmoothed[i].MeanBulkDP, smAF, rep)
	}
	color.Green("\nThreshold calculations complete.")

	color.Cyan("\n=========================================== Writing Smoothed TSV =================================================\n\n")
	smoothTSV := filepath.Join(outDir, "GoBSAseq.smooth.tsv")
	if err := writeSmoothedTSVOneBulk(smoothTSV, allSmoothed, smAF, rep); err != nil {
		color.Red("Error writing smoothed TSV: %v", err)
	} else {
		color.Green("Wrote %d smoothed windows to %s", len(allSmoothed), smoothTSV)
	}
	color.Green("Raw stats written to %s", filepath.Join(outDir, "GoBSAseq.raw.tsv"))
	color.Green("\nStats pipeline time: %s\n", time.Since(overallStart).Round(time.Second))

	color.Cyan("\n============================ Generating HTML Plots & QTLs ========================================\n\n")
	htmlFile := filepath.Join(outDir, "GoBSAseq_RobustZScore.html")
	qtlFile := filepath.Join(outDir, "GoBSAseq_QTL.tsv")

	finalQTLs, err := GenerateHtmlPlotsAndQTLOneBulk(allSmoothed, bulkSize, pop, cfg.Alphas, htmlFile, qtlFile)
	if err != nil {
		return fmt.Errorf("plots and QTLs: %w", err)
	}
	color.Green("HTML plots written to %s", outDir)
	color.Green("QTL tabular results written to %s", qtlFile)

	geneSpaceEnabled := cfg.SnpEffDB != "" && cfg.GeneDesc != "" && cfg.Prg != "" && cfg.Gff != ""
	if !geneSpaceEnabled {
		color.Yellow("Skipping gene space analysis: required parameters were not provided.")
		return nil
	}

	color.Green("Step 1: Annotating genes with SnpEff (creating Super VCF table)")
	_, hasEFF := cfg.Rdr.Header.Infos["EFF"]
	if hasEFF {
		color.Yellow("EFF column is already present; skipping gene space analysis.")
		return nil
	}
	err, annotatedTsvFiles := annotation.CreateSuperVcf([]string{filteredVcfPath}, cfg.SnpEffDB, true, cfg.GeneDesc, cfg.Prg)
	if err != nil {
		return fmt.Errorf("snpeff annotation: %w", err)
	}
	if len(annotatedTsvFiles) == 0 {
		color.Yellow("Skipping gene space analysis: SnpEff did not produce an annotated TSV.")
		return nil
	}

	color.Cyan("Performing Gene space analysis for %d QTL intervals ...", len(finalQTLs))
	for _, qtl := range finalQTLs {
		_, err = genespace.GeneSpace(
			cfg.Gff, annotatedTsvFiles[0], qtl.Chrom,
			int(qtl.Start), int(qtl.Stop), []string{cfg.HighParentName},
			[]string{cfg.LowParentName}, cfg.GeneDesc, cfg.Prg, outDir)
		if err != nil {
			return fmt.Errorf("gene space: %w", err)
		}
	}
	color.Green("Gene space analysis complete.")
	return nil
}
