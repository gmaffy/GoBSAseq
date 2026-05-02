package twobulk

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
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
	CHROM   string
	POS     int64
	DeltaSI float64
	Gstat   float64
	ED      float64
	LOD     float64
	BBLogBF float64
	NumSNPs int
}

// Thresholds holds the genome-wide significance thresholds derived from
// max-statistic permutation testing. Each value is the (1-alpha) quantile
// of the null distribution of per-permutation genome-wide maxima.
type Thresholds struct {
	DeltaSI float64
	Gstat   float64
	ED      float64
	LOD     float64
	BBLogBF float64
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

func EuclideanDist(refBulk1, altBulk1, refBulk2, altBulk2 int) float64 {
	total1 := float64(refBulk1 + altBulk1)
	total2 := float64(refBulk2 + altBulk2)
	if total1 == 0 || total2 == 0 {
		return 0
	}

	pAlt := float64(altBulk1) / total1
	qAlt := float64(altBulk2) / total2

	return math.Abs(pAlt - qAlt)
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

		var sumDeltaSI, sumED float64
		var sumWeightDeltaSI, sumWeightED float64
		var sumGstat, sumLOD, sumBBLogBF float64
		var nSNPs int

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
		}

		if nSNPs == 0 {
			continue
		}

		sm := SmoothedStats{
			CHROM:   chrom,
			POS:     center,
			NumSNPs: nSNPs,
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

		smoothed = append(smoothed, sm)
	}

	return smoothed
}

// permuteSNPs returns a copy of stats with the high/low bulk read counts
// randomly swapped per-SNP, producing a null dataset. Positions and depths
// are preserved so the smoothing geometry is identical to the real data.
func permuteSNPs(stats []BSAstats, rng *rand.Rand) []BSAstats {
	perm := make([]BSAstats, len(stats))
	copy(perm, stats)
	for i := range perm {
		if rng.Intn(2) == 0 {
			perm[i].HighBulkL, perm[i].LowBulkL = perm[i].LowBulkL, perm[i].HighBulkL
			perm[i].HighBulkH, perm[i].LowBulkH = perm[i].LowBulkH, perm[i].HighBulkH
		}
		// recalculate all derived stats from the (possibly swapped) counts
		hbDP := perm[i].HighBulkL + perm[i].HighBulkH
		lbDP := perm[i].LowBulkL + perm[i].LowBulkH
		if hbDP == 0 || lbDP == 0 {
			perm[i].HighSI, perm[i].LowSI = 0, 0
			perm[i].DeltaSI, perm[i].Gstat, perm[i].ED, perm[i].LOD, perm[i].BBLogBF = 0, 0, 0, 0, 0
			continue
		}
		hSI := float64(perm[i].HighBulkH) / float64(hbDP)
		lSI := float64(perm[i].LowBulkH) / float64(lbDP)
		perm[i].HighSI = hSI
		perm[i].LowSI = lSI
		perm[i].DeltaSI = hSI - lSI
		perm[i].ED = math.Abs(hSI - lSI)
		perm[i].Gstat = GStatistic(perm[i].HighBulkH, perm[i].HighBulkL, perm[i].LowBulkH, perm[i].LowBulkL)
		perm[i].LOD = lod(perm[i].HighBulkL, perm[i].HighBulkH, perm[i].LowBulkL, perm[i].LowBulkH)
		perm[i].BBLogBF = betaBinomialLogBF(perm[i].HighBulkH, perm[i].HighBulkL, perm[i].LowBulkH, perm[i].LowBulkL)
	}
	return perm
}

// quantile returns the value at quantile q (0–1) of a sorted slice using
// linear interpolation.
func quantile(sorted []float64, q float64) float64 {
	if len(sorted) == 0 {
		return 0
	}
	idx := q * float64(len(sorted)-1)
	lo := int(math.Floor(idx))
	hi := int(math.Ceil(idx))
	if lo == hi {
		return sorted[lo]
	}
	frac := idx - float64(lo)
	return sorted[lo]*(1-frac) + sorted[hi]*frac
}

// calcThresholds runs nReps permutations in parallel. Each permutation records
// the genome-wide maximum smoothed value for each stat (max-statistic approach),
// which controls FWER across all windows without a separate multiple-testing
// correction step. The (1-alpha) quantile of those maxima is returned.
func calcThresholds(chromStats map[string][]BSAstats, windowSize, stepSize int64, nReps int, alpha float64) Thresholds {
	type maxima struct {
		deltaSI, gstat, ed, lod, bbLogBF float64
	}

	// Flatten all SNPs for permutation; chromosome labels are preserved inside
	// each BSAstats so smoothChromosome can re-partition them after the swap.
	var allStats []BSAstats
	for _, s := range chromStats {
		allStats = append(allStats, s...)
	}

	results := make([]maxima, nReps)
	repChan := make(chan int, nReps)
	for r := 0; r < nReps; r++ {
		repChan <- r
	}
	close(repChan)

	numWorkers := runtime.NumCPU()
	var wg sync.WaitGroup
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			rng := rand.New(rand.NewSource(time.Now().UnixNano()))
			for r := range repChan {
				permAll := permuteSNPs(allStats, rng)

				// re-partition permuted SNPs back into chromosomes
				permChrom := make(map[string][]BSAstats)
				for _, s := range permAll {
					permChrom[s.CHROM] = append(permChrom[s.CHROM], s)
				}

				var mx maxima
				for _, stats := range permChrom {
					for _, sm := range smoothChromosome(stats, windowSize, stepSize) {
						if math.Abs(sm.DeltaSI) > mx.deltaSI {
							mx.deltaSI = math.Abs(sm.DeltaSI)
						}
						if sm.Gstat > mx.gstat {
							mx.gstat = sm.Gstat
						}
						if sm.ED > mx.ed {
							mx.ed = sm.ED
						}
						if sm.LOD > mx.lod {
							mx.lod = sm.LOD
						}
						if sm.BBLogBF > mx.bbLogBF {
							mx.bbLogBF = sm.BBLogBF
						}
					}
				}
				results[r] = mx
			}
		}()
	}
	wg.Wait()

	// Collect each stat's maxima into sorted slices then take the quantile
	ds := make([]float64, nReps)
	gs := make([]float64, nReps)
	ed := make([]float64, nReps)
	ls := make([]float64, nReps)
	bb := make([]float64, nReps)
	for i, m := range results {
		ds[i] = m.deltaSI
		gs[i] = m.gstat
		ed[i] = m.ed
		ls[i] = m.lod
		bb[i] = m.bbLogBF
	}
	sort.Float64s(ds)
	sort.Float64s(gs)
	sort.Float64s(ed)
	sort.Float64s(ls)
	sort.Float64s(bb)

	q := 1.0 - alpha
	return Thresholds{
		DeltaSI: quantile(ds, q),
		Gstat:   quantile(gs, q),
		ED:      quantile(ed, q),
		LOD:     quantile(ls, q),
		BBLogBF: quantile(bb, q),
	}
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

// writeSmoothedTSV writes smoothed windows with thresholds as comment lines at
// the top of the file. Comment lines are prefixed with # so standard TSV parsers
// and R/Python tools that skip comments can read the data columns unchanged.
func writeSmoothedTSV(filename string, data []SmoothedStats, thresh Thresholds) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	fmt.Fprintf(w, "#Threshold_DeltaSI\t%.6f\n", thresh.DeltaSI)
	fmt.Fprintf(w, "#Threshold_Gstat\t%.6f\n", thresh.Gstat)
	fmt.Fprintf(w, "#Threshold_ED\t%.6f\n", thresh.ED)
	fmt.Fprintf(w, "#Threshold_LOD\t%.6f\n", thresh.LOD)
	fmt.Fprintf(w, "#Threshold_BBLogBF\t%.6f\n", thresh.BBLogBF)

	fmt.Fprintln(w, "CHROM\tPOS\tDeltaSI\tGstat\tED\tLOD\tBBLogBF\tNumSNPs")

	for _, d := range data {
		fmt.Fprintf(w, "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
			d.CHROM, d.POS, d.DeltaSI, d.Gstat, d.ED, d.LOD, d.BBLogBF, d.NumSNPs)
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
	outPrefix := cfg.OutputFile

	windowSize := int64(cfg.WindowSize)
	stepSize := int64(cfg.StepSize)
	nReps := cfg.Rep
	alpha := cfg.Alpha

	overallStart := time.Now()
	variantChan := make(chan *vcfgo.Variant, 10000)
	statsChan := make(chan BSAstats, 10000)
	rawWriteChan := make(chan BSAstats, 10000)

	numWorkers := runtime.NumCPU()
	var workerWG sync.WaitGroup

	// Start raw TSV writer goroutine
	var rawWG sync.WaitGroup
	rawWG.Add(1)
	go func() {
		defer rawWG.Done()
		if err := writeRawTSV(outPrefix+".raw.tsv", rawWriteChan); err != nil {
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
				lpAllele := 0
				if len(lpGT) > 0 {
					lpAllele = lpGT[0]
				}
				hbAllele := 0
				if len(hbGT) > 0 {
					hbAllele = hbGT[0]
				}
				lbAllele := 0
				if len(lbGT) > 0 {
					lbAllele = lbGT[0]
				}

				if hbAllele == lpAllele {
					hbL, hbH = hbRefDep, hbAltDeps[0]
				} else {
					hbL, hbH = hbAltDeps[0], hbRefDep
				}
				if lbAllele == lpAllele {
					lbL, lbH = lbRefDep, lbAltDeps[0]
				} else {
					lbL, lbH = lbAltDeps[0], lbRefDep
				}

				hSI, lSI := float64(hbH)/float64(hbDP), float64(lbH)/float64(lbDP)
				minDepth := hbDP
				if lbDP < minDepth {
					minDepth = lbDP
				}

				s := BSAstats{
					CHROM: variant.Chromosome, POS: int64(variant.Pos), REF: variant.Reference, ALT: variant.Alt()[0],
					HighParGT: variant.Samples[highParIdx].GT, LowParGT: variant.Samples[lowParIdx].GT,
					HighBulkGT: variant.Samples[highBulkIdx].GT, HighBulkAD: fmt.Sprintf("%d,%d", hbRefDep, hbAltDeps[0]),
					LowBulkGT: variant.Samples[lowBulkIdx].GT, LowBulkAD: fmt.Sprintf("%d,%d", lbRefDep, lbAltDeps[0]),
					HighBulkL: hbL, HighBulkH: hbH, LowBulkL: lbL, LowBulkH: lbH,
					HighSI: hSI, LowSI: lSI, DeltaSI: hSI - lSI,
					Gstat: GStatistic(hbH, hbL, lbH, lbL),
					ED:    math.Abs(hSI - lSI), //EuclideanDist(hbL, hbH, lbL, lbH),
					LOD:   lod(hbL, hbH, lbL, lbH), BBLogBF: betaBinomialLogBF(hbH, hbL, lbH, lbL),
					Depth: minDepth,
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

	color.Cyan("\n============================ Calculating Thresholds (%d permutations, alpha=%.3f) =============================\n\n", nReps, alpha)
	thresh := calcThresholds(chromStats, windowSize, stepSize, nReps, alpha)

	color.Green("  DeltaSI threshold : %.6f", thresh.DeltaSI)
	color.Green("  Gstat   threshold : %.6f", thresh.Gstat)
	color.Green("  ED      threshold : %.6f", thresh.ED)
	color.Green("  LOD     threshold : %.6f", thresh.LOD)
	color.Green("  BBLogBF threshold : %.6f", thresh.BBLogBF)

	color.Cyan("\n============================ Writing Smoothed TSV =============================\n\n")
	if err := writeSmoothedTSV(outPrefix+".smooth.tsv", allSmoothed, thresh); err != nil {
		color.Red("Error writing smoothed TSV: %v", err)
	} else {
		color.Green("Wrote %d smoothed windows to %s.smooth.tsv", len(allSmoothed), outPrefix)
	}
	color.Green("Raw stats written to %s.raw.tsv", outPrefix)

	color.Green("\nTotal time: %s\n", time.Since(overallStart).Round(time.Second))
}
