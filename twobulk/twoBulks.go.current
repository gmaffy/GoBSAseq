package twobulk

import (
	"fmt"
	"html"
	"math"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
)

// --- Types ---

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

	DeltaSIK float64
	Gprime   float64
	EDK      float64

	BBLogBF float64
	BBK     float64
	EWMA    float64
	CUSUM   float64
	ScanLLR float64
	SegMean float64
}

type Block struct {
	Chr   string
	Start int64
	End   int64
	AFD   float64
	Z     float64
	Sig   bool
}

type ScanWindow struct {
	Chr      string
	Center   int64
	Start    int64
	End      int64
	SNPs     int
	LLR      float64
	HighProp float64
	LowProp  float64
}

type ChangePoint struct {
	Chr       string
	Pos       int64
	LeftMean  float64
	RightMean float64
	Score     float64
}

type QTL struct {
	Chr       string
	Start     int64
	End       int64
	PeakPos   int64
	PeakScore float64
	Method    string
}

type plotSeries struct {
	Label  string
	Color  string
	Values []float64
}

// --- Analysis Engine ---

type Analyzer struct {
	WindowSize int
	HighPar    int
	HighParDP  int
	LowPar     int
	LowParDP   int
	HighBulk   int
	HighBulkDP int
	LowBulk    int
	LowBulkDP  int

	Stats      []BSAstats
	ByChr      map[string][]BSAstats
	Scans      []ScanWindow
	Changes    []ChangePoint
	Blocks     []Block
	QTLs       []QTL
}

func NewAnalyzer(windowSize int) *Analyzer {
	return &Analyzer{
		WindowSize: windowSize,
		ByChr:      make(map[string][]BSAstats),
	}
}

func (a *Analyzer) Run(vcfRdr *vcfgo.Reader) {
	overallStart := time.Now()
	color.New(color.FgHiWhite, color.Bold).Println("\n==================== TWO-BULK BSA ANALYSIS ====================")
	fmt.Printf("Workers: %d | Kernel window: %d bp\n", runtime.NumCPU(), a.WindowSize)

	a.CollectStats(vcfRdr)
	a.GroupAndSort()
	a.ProcessChromosomes()
	a.ComputeBlocks()
	a.IdentifyQTLs()
	a.WriteOutputs()

	color.New(color.FgHiWhite, color.Bold).Printf("\nAnalysis complete in %s.\n", time.Since(overallStart).Round(time.Millisecond))
}

func (a *Analyzer) CollectStats(vcfRdr *vcfgo.Reader) {
	stage := stageStart("VCF scan and per-SNP statistics")
	var total, passed int
	var mu sync.Mutex
	variantChan := make(chan *vcfgo.Variant, 1000)
	var wg sync.WaitGroup

	for i := 0; i < runtime.NumCPU(); i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for v := range variantChan {
				if !a.isInformative(v) {
					continue
				}
				stat := a.computeSNPStats(v)
				if stat != nil {
					mu.Lock()
					a.Stats = append(a.Stats, *stat)
					passed++
					mu.Unlock()
				}
			}
		}()
	}

	for v := vcfRdr.Read(); v != nil; v = vcfRdr.Read() {
		total++
		variantChan <- v
	}
	close(variantChan)
	wg.Wait()
	stageDone("VCF scan and per-SNP statistics", stage, "Scanned %d variants, retained %d SNPs.", total, passed)
}

func (a *Analyzer) isInformative(v *vcfgo.Variant) bool {
	if len(v.Alt()) != 1 {
		return false
	}
	indices := []int{a.HighPar, a.LowPar, a.HighBulk, a.LowBulk}
	for _, idx := range indices {
		if idx >= len(v.Samples) || len(v.Samples[idx].GT) == 0 {
			return false
		}
		for _, allele := range v.Samples[idx].GT {
			if allele < 0 {
				return false
			}
		}
	}

	hp := v.Samples[a.HighPar]
	lp := v.Samples[a.LowPar]
	hb := v.Samples[a.HighBulk]
	lb := v.Samples[a.LowBulk]

	if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) || hp.GT[0] == lp.GT[0] {
		return false
	}
	if hp.DP < a.HighParDP || lp.DP < a.LowParDP || hb.DP < a.HighBulkDP || lb.DP < a.LowBulkDP {
		return false
	}
	return true
}

func (a *Analyzer) computeSNPStats(v *vcfgo.Variant) *BSAstats {
	hb := v.Samples[a.HighBulk]
	lb := v.Samples[a.LowBulk]
	if hb.DP == 0 || lb.DP == 0 {
		return nil
	}

	hbRef, _ := hb.RefDepth()
	hbAlts, _ := hb.AltDepths()
	lbRef, _ := lb.RefDepth()
	lbAlts, _ := lb.AltDepths()
	if len(hbAlts) == 0 || len(lbAlts) == 0 {
		return nil
	}

	lpGT := v.Samples[a.LowPar].GT
	var hbL, hbH, lbL, lbH int
	if hb.GT[0] == lpGT[0] {
		hbL, hbH = hbRef, hbAlts[0]
	} else {
		hbL, hbH = hbAlts[0], hbRef
	}
	if lb.GT[0] == lpGT[0] {
		lbL, lbH = lbRef, lbAlts[0]
	} else {
		lbL, lbH = lbAlts[0], lbRef
	}

	hSI := safeRatio(hbH, hb.DP)
	lSI := safeRatio(lbH, lb.DP)
	deltaSI := hSI - lSI

	return &BSAstats{
		CHROM:      v.Chromosome,
		POS:        int64(v.Pos),
		REF:        v.Reference,
		ALT:        v.Alt()[0],
		HighParGT:  v.Samples[a.HighPar].GT,
		LowParGT:   v.Samples[a.LowPar].GT,
		HighBulkGT: hb.GT,
		HighBulkAD: fmt.Sprintf("%v,%v", hbRef, hbAlts[0]),
		LowBulkGT:  lb.GT,
		LowBulkAD:  fmt.Sprintf("%v,%v", lbRef, lbAlts[0]),
		HighBulkL:  hbL,
		HighBulkH:  hbH,
		LowBulkL:   lbL,
		LowBulkH:   lbH,
		HighSI:     hSI,
		LowSI:      lSI,
		DeltaSI:    deltaSI,
		Gstat:      GStatistic(hbAlts[0], hbRef, lbAlts[0], lbRef),
		ED:         math.Pow(deltaSI, 4), // ED^4 update
		BBLogBF:    BetaBinomialLogBayesFactor(hbH, hbL, lbH, lbL, 0.5, 0.5),
	}
}

func (a *Analyzer) GroupAndSort() {
	for _, s := range a.Stats {
		a.ByChr[s.CHROM] = append(a.ByChr[s.CHROM], s)
	}
	for chr := range a.ByChr {
		sort.Slice(a.ByChr[chr], func(i, j int) bool {
			return a.ByChr[chr][i].POS < a.ByChr[chr][j].POS
		})
	}
}

func (a *Analyzer) ProcessChromosomes() {
	stage := stageStart("Smoothing, scan statistics, and changepoints")
	winSizes := uniqueWindowSizes(a.WindowSize)

	for chr, snps := range a.ByChr {
		smoothChromosome(snps, a.WindowSize)
		computeEWMA(snps, 0.2)
		computeCUSUM(snps, 0.5)
		a.Scans = append(a.Scans, scanChromosome(snps, winSizes)...)
		a.Changes = append(a.Changes, detectChangePoints(snps, 12)...)
		fmt.Printf("  %-18s SNPs=%6d | scan peaks=%4d | changepoints=%3d\n",
			chr, len(snps), len(filterScansByChr(a.Scans, chr)), len(filterChangesByChr(a.Changes, chr)))
	}
	stageDone("Smoothing, scan statistics, and changepoints", stage, "Processed %d chromosomes.", len(a.ByChr))
}

func (a *Analyzer) ComputeBlocks() {
	a.Blocks = brmBlocks(a.Stats, 100000, 3.0)
}

func (a *Analyzer) IdentifyQTLs() {
	// Simple QTL identification based on Gprime and EDK thresholds
	// For Gprime, use a threshold of 5.0 (approx p < 0.01)
	// For EDK, we use a top percentile or fixed threshold.
	const gThresh = 5.0
	const edThresh = 0.05 // Example threshold for ED^4

	for chr, snps := range a.ByChr {
		// Identify contiguous regions above Gprime threshold
		inQTL := false
		var currentQTL QTL
		for _, s := range snps {
			if s.Gprime >= gThresh {
				if !inQTL {
					inQTL = true
					currentQTL = QTL{Chr: chr, Start: s.POS, Method: "Gprime"}
					currentQTL.PeakScore = s.Gprime
					currentQTL.PeakPos = s.POS
				} else {
					if s.Gprime > currentQTL.PeakScore {
						currentQTL.PeakScore = s.Gprime
						currentQTL.PeakPos = s.POS
					}
				}
				currentQTL.End = s.POS
			} else {
				if inQTL {
					a.QTLs = append(a.QTLs, currentQTL)
					inQTL = false
				}
			}
		}
		if inQTL {
			a.QTLs = append(a.QTLs, currentQTL)
		}

		// Identify regions above ED^4 threshold
		inQTL = false
		for _, s := range snps {
			if s.EDK >= edThresh {
				if !inQTL {
					inQTL = true
					currentQTL = QTL{Chr: chr, Start: s.POS, Method: "ED4"}
					currentQTL.PeakScore = s.EDK
					currentQTL.PeakPos = s.POS
				} else {
					if s.EDK > currentQTL.PeakScore {
						currentQTL.PeakScore = s.EDK
						currentQTL.PeakPos = s.POS
					}
				}
				currentQTL.End = s.POS
			} else {
				if inQTL {
					a.QTLs = append(a.QTLs, currentQTL)
					inQTL = false
				}
			}
		}
		if inQTL {
			a.QTLs = append(a.QTLs, currentQTL)
		}
	}
}

func (a *Analyzer) WriteOutputs() {
	// snps.tsv
	writeSNPs("snps.tsv", a.Stats)
	// brm_blocks.tsv
	writeBlocks("brm_blocks.tsv", a.Blocks)
	// scan_windows.tsv
	writeScanWindows("scan_windows.tsv", a.Scans)
	// changepoints.tsv
	writeChangePoints("changepoints.tsv", a.Changes)
	// qtls.tsv
	writeQTLs("qtls.tsv", a.QTLs)
	// plots.html
	writePlots("plots.html", a.ByChr, a.Blocks, a.Scans, a.Changes)
}

// --- Helper Functions (Consolidated/Reduced) ---

func stageStart(label string) time.Time {
	color.New(color.FgCyan, color.Bold).Printf("\n[%s] Starting...\n", label)
	return time.Now()
}

func stageDone(label string, start time.Time, format string, args ...interface{}) {
	elapsed := time.Since(start).Round(time.Millisecond)
	detail := ""
	if format != "" {
		detail = fmt.Sprintf(format, args...)
	}
	if detail != "" {
		color.New(color.FgGreen).Printf("[%s] Done in %s. %s\n", label, elapsed, detail)
		return
	}
	color.New(color.FgGreen).Printf("[%s] Done in %s.\n", label, elapsed)
}

func isHomozygous(gt []int) bool {
	for _, a := range gt[1:] {
		if a != gt[0] {
			return false
		}
	}
	return true
}

func safeRatio(num, den int) float64 {
	if den == 0 {
		return 0
	}
	return float64(num) / float64(den)
}

func GStatistic(hAlt, hRef, lAlt, lRef int) float64 {
	hTot := float64(hAlt + hRef)
	lTot := float64(lAlt + lRef)
	tot := hTot + lTot
	if hTot == 0 || lTot == 0 {
		return 0
	}

	expHAlt := hTot * float64(hAlt+lAlt) / tot
	expHRef := hTot * float64(hRef+lRef) / tot
	expLAlt := lTot * float64(hAlt+lAlt) / tot
	expLRef := lTot * float64(hRef+lRef) / tot

	g := 0.0
	term := func(obs int, exp float64) {
		if obs > 0 && exp > 0 {
			g += float64(obs) * math.Log(float64(obs)/exp)
		}
	}
	term(hAlt, expHAlt)
	term(hRef, expHRef)
	term(lAlt, expLAlt)
	term(lRef, expLRef)
	return 2 * g
}

func BetaBinomialLogBayesFactor(hS, hF, lS, lF int, alpha, beta float64) float64 {
	lBeta := func(a, b float64) float64 {
		ga, _ := math.Lgamma(a)
		gb, _ := math.Lgamma(b)
		gab, _ := math.Lgamma(a + b)
		return ga + gb - gab
	}
	logPrior := lBeta(alpha, beta)
	logAlt := (lBeta(alpha+float64(hS), beta+float64(hF)) - logPrior) + (lBeta(alpha+float64(lS), beta+float64(lF)) - logPrior)
	logNull := lBeta(alpha+float64(hS+lS), beta+float64(hF+lF)) - logPrior
	return logAlt - logNull
}

func smoothChromosome(snps []BSAstats, bandwidth int) {
	if bandwidth <= 0 {
		bandwidth = 1
	}
	biweight := func(x float64) float64 {
		if x >= 1 {
			return 0
		}
		t := 1 - x*x
		return t * t
	}
	tricube := func(x float64) float64 {
		if x >= 1 {
			return 0
		}
		return math.Pow(1-math.Pow(x, 3), 3)
	}

	for i := range snps {
		center := snps[i].POS
		var nD, nE, nG, nBB, dD, dE, dG, dBB float64

		// Search backwards
		for j := i; j >= 0; j-- {
			dist := float64(center - snps[j].POS)
			if dist > float64(bandwidth) {
				break
			}
			x := dist / float64(bandwidth)
			wB := biweight(x)
			wT := tricube(x)

			nD += wB * snps[j].DeltaSI
			nE += wB * snps[j].ED
			nBB += wB * snps[j].BBLogBF
			nG += wT * snps[j].Gstat
			dD += wB
			dE += wB
			dBB += wB
			dG += wT
		}
		// Search forwards
		for j := i + 1; j < len(snps); j++ {
			dist := float64(snps[j].POS - center)
			if dist > float64(bandwidth) {
				break
			}
			x := dist / float64(bandwidth)
			wB := biweight(x)
			wT := tricube(x)

			nD += wB * snps[j].DeltaSI
			nE += wB * snps[j].ED
			nBB += wB * snps[j].BBLogBF
			nG += wT * snps[j].Gstat
			dD += wB
			dE += wB
			dBB += wB
			dG += wT
		}

		if dD > 0 {
			snps[i].DeltaSIK = nD / dD
			snps[i].EDK = nE / dE
		}
		if dBB > 0 {
			snps[i].BBK = nBB / dBB
		}
		if dG > 0 {
			snps[i].Gprime = nG / dG
		}
	}
}

func computeEWMA(snps []BSAstats, alpha float64) {
	if len(snps) == 0 {
		return
	}
	snps[0].EWMA = snps[0].DeltaSI
	for i := 1; i < len(snps); i++ {
		snps[i].EWMA = alpha*snps[i].DeltaSI + (1-alpha)*snps[i-1].EWMA
	}
}

func computeCUSUM(snps []BSAstats, reference float64) {
	if len(snps) == 0 {
		return
	}
	var sum float64
	for _, s := range snps {
		sum += s.DeltaSI
	}
	mean := sum / float64(len(snps))
	var ss float64
	for _, s := range snps {
		d := s.DeltaSI - mean
		ss += d * d
	}
	sd := math.Sqrt(ss / float64(len(snps)))
	if sd == 0 {
		return
	}

	pos, neg := 0.0, 0.0
	for i := range snps {
		x := (snps[i].DeltaSI - mean) / sd
		pos = math.Max(0, pos+x-reference)
		neg = math.Min(0, neg+x+reference)
		if pos >= -neg {
			snps[i].CUSUM = pos
		} else {
			snps[i].CUSUM = neg
		}
	}
}

func brmBlocks(snps []BSAstats, blockSize int64, zThresh float64) []Block {
	byChr := make(map[string][]BSAstats)
	for _, s := range snps {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}

	var blocks []Block
	for chr, list := range byChr {
		sort.Slice(list, func(i, j int) bool { return list[i].POS < list[j].POS })
		for i := 0; i < len(list); {
			start := list[i].POS
			end := start + blockSize
			var sum float64
			var n int
			blockStartIdx := i
			for i < len(list) && list[i].POS < end {
				sum += list[i].DeltaSI
				n++
				i++
			}
			if n == 0 {
				continue
			}
			afd := sum / float64(n)
			var variance float64
			for k := blockStartIdx; k < i; k++ {
				diff := list[k].DeltaSI - afd
				variance += diff * diff
			}
			z := 0.0
			if n > 1 {
				std := math.Sqrt(variance / float64(n-1))
				if std > 0 {
					z = afd / (std / math.Sqrt(float64(n)))
				}
			}
			blocks = append(blocks, Block{Chr: chr, Start: start, End: end, AFD: afd, Z: z, Sig: math.Abs(z) >= zThresh})
		}
	}
	return blocks
}

func scanChromosome(snps []BSAstats, windowSizes []int64) []ScanWindow {
	n := len(snps)
	if n < 2 {
		return nil
	}
	pos := make([]int64, n)
	hS, hT, lS, lT := make([]int, n+1), make([]int, n+1), make([]int, n+1), make([]int, n+1)
	for i, s := range snps {
		pos[i] = s.POS
		hS[i+1] = hS[i] + s.HighBulkH
		hT[i+1] = hT[i] + s.HighBulkH + s.HighBulkL
		lS[i+1] = lS[i] + s.LowBulkH
		lT[i+1] = lT[i] + s.LowBulkH + s.LowBulkL
	}

	bestWin := make([]ScanWindow, n)
	for i, s := range snps {
		best := ScanWindow{Chr: s.CHROM, Center: s.POS}
		for _, size := range windowSizes {
			lPos, rPos := s.POS-size/2, s.POS+size/2
			lIdx := sort.Search(n, func(j int) bool { return pos[j] >= lPos })
			rIdx := sort.Search(n, func(j int) bool { return pos[j] > rPos })
			if rIdx-lIdx < 2 {
				continue
			}
			hSucc, hTot := hS[rIdx]-hS[lIdx], hT[rIdx]-hT[lIdx]
			lSucc, lTot := lS[rIdx]-lS[lIdx], lT[rIdx]-lT[lIdx]
			llr := binomialLLR(hSucc, hTot, lSucc, lTot)
			if llr > best.LLR {
				best = ScanWindow{Chr: s.CHROM, Center: s.POS, Start: pos[lIdx], End: pos[rIdx-1], SNPs: rIdx - lIdx, LLR: llr, HighProp: safeRatio(hSucc, hTot), LowProp: safeRatio(lSucc, lTot)}
			}
		}
		snps[i].ScanLLR = best.LLR
		bestWin[i] = best
	}

	var windows []ScanWindow
	for i, w := range bestWin {
		if w.LLR <= 0 {
			continue
		}
		if (i == 0 || w.LLR >= bestWin[i-1].LLR) && (i == n-1 || w.LLR >= bestWin[i+1].LLR) {
			windows = append(windows, w)
		}
	}
	return windows
}

func binomialLLR(hS, hT, lS, lT int) float64 {
	if hT == 0 || lT == 0 {
		return 0
	}
	logL := func(s, t int) float64 {
		if t == 0 {
			return 0
		}
		p := float64(s) / float64(t)
		res := 0.0
		if s > 0 && p > 0 {
			res += float64(s) * math.Log(p)
		}
		if t-s > 0 && 1-p > 0 {
			res += float64(t-s) * math.Log(1-p)
		}
		return res
	}
	llAlt := logL(hS, hT) + logL(lS, lT)
	llNull := logL(hS+lS, hT+lT)
	return math.Max(0, 2*(llAlt-llNull))
}

func detectChangePoints(snps []BSAstats, minSize int) []ChangePoint {
	n := len(snps)
	if n < 2*minSize {
		return nil
	}
	vals := make([]float64, n)
	for i, s := range snps {
		vals[i] = s.DeltaSIK
	}
	sum, sumSq := make([]float64, n+1), make([]float64, n+1)
	for i, v := range vals {
		sum[i+1], sumSq[i+1] = sum[i]+v, sumSq[i]+v*v
	}
	sse := func(lo, hi int) float64 {
		if hi-lo <= 1 {
			return 0
		}
		s, s2 := sum[hi]-sum[lo], sumSq[hi]-sumSq[lo]
		return s2 - (s*s)/float64(hi-lo)
	}

	var ssTotal float64
	m := sum[n] / float64(n)
	for _, v := range vals {
		ssTotal += (v - m) * (v - m)
	}
	sd := math.Sqrt(ssTotal / float64(n))
	penalty := math.Max(0.02, 2*sd*sd*math.Log(float64(n)+1))

	var cps []int
	var split func(int, int)
	split = func(lo, hi int) {
		if hi-lo < 2*minSize {
			return
		}
		base, bestIdx, bestGain := sse(lo, hi), -1, 0.0
		for k := lo + minSize; k <= hi-minSize; k++ {
			gain := base - sse(lo, k) - sse(k, hi)
			if gain > bestGain {
				bestGain, bestIdx = gain, k
			}
		}
		if bestIdx != -1 && bestGain > penalty {
			cps = append(cps, bestIdx)
			split(lo, bestIdx)
			split(bestIdx, hi)
		}
	}
	split(0, n)
	sort.Ints(cps)

	bounds := append([]int{0}, cps...)
	bounds = append(bounds, n)
	var changes []ChangePoint
	for i := 0; i < len(bounds)-1; i++ {
		lo, hi := bounds[i], bounds[i+1]
		mean := (sum[hi] - sum[lo]) / float64(hi-lo)
		for j := lo; j < hi; j++ {
			snps[j].SegMean = mean
		}
		if i < len(bounds)-2 {
			cpIdx := bounds[i+1]
			nextMean := (sum[bounds[i+2]] - sum[bounds[i+1]]) / float64(bounds[i+2]-bounds[i+1])
			changes = append(changes, ChangePoint{Chr: snps[cpIdx].CHROM, Pos: snps[cpIdx].POS, LeftMean: mean, RightMean: nextMean, Score: math.Abs(nextMean - mean)})
		}
	}
	return changes
}

func uniqueWindowSizes(base int) []int64 {
	sizes := []int64{int64(base / 4), int64(base / 2), int64(base)}
	m := make(map[int64]bool)
	var res []int64
	for _, s := range sizes {
		if s > 0 && !m[s] {
			m[s], res = true, append(res, s)
		}
	}
	sort.Slice(res, func(i, j int) bool { return res[i] < res[j] })
	return res
}

func filterScansByChr(scans []ScanWindow, chr string) []ScanWindow {
	var res []ScanWindow
	for _, s := range scans {
		if s.Chr == chr {
			res = append(res, s)
		}
	}
	return res
}

func filterChangesByChr(changes []ChangePoint, chr string) []ChangePoint {
	var res []ChangePoint
	for _, c := range changes {
		if c.Chr == chr {
			res = append(res, c)
		}
	}
	return res
}

func filterBlocksByChr(blocks []Block, chr string) []Block {
	var res []Block
	for _, b := range blocks {
		if b.Chr == chr {
			res = append(res, b)
		}
	}
	return res
}

// --- IO Functions ---

func writeSNPs(path string, stats []BSAstats) {
	f, _ := os.Create(path)
	defer f.Close()
	fmt.Fprintln(f, "CHR\tPOS\tSNPiH\tSNPiL\tDELTA\tDELTAK\tED\tEDK\tG\tGPRIME\tBBLOGBF\tBBK\tEWMA\tCUSUM\tSCANLLR\tSEGMEAN")
	for _, s := range stats {
		fmt.Fprintf(f, "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
			s.CHROM, s.POS, s.HighSI, s.LowSI, s.DeltaSI, s.DeltaSIK, s.ED, s.EDK, s.Gstat, s.Gprime, s.BBLogBF, s.BBK, s.EWMA, s.CUSUM, s.ScanLLR, s.SegMean)
	}
}

func writeBlocks(path string, blocks []Block) {
	f, _ := os.Create(path)
	defer f.Close()
	fmt.Fprintln(f, "CHR\tSTART\tEND\tAFD\tZ\tSIGNIFICANT")
	for _, b := range blocks {
		fmt.Fprintf(f, "%s\t%d\t%d\t%.4f\t%.2f\t%v\n", b.Chr, b.Start, b.End, b.AFD, b.Z, b.Sig)
	}
}

func writeScanWindows(path string, scans []ScanWindow) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "CHR\tCENTER\tSTART\tEND\tSNPS\tLLR\tHIGHPROP\tLOWPROP")
	for _, w := range scans {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n", w.Chr, w.Center, w.Start, w.End, w.SNPs, w.LLR, w.HighProp, w.LowProp)
	}
	return nil
}

func writeChangePoints(path string, changes []ChangePoint) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "CHR\tPOS\tLEFTMEAN\tRIGHTMEAN\tSHIFT")
	for _, c := range changes {
		fmt.Fprintf(f, "%s\t%d\t%.4f\t%.4f\t%.4f\n", c.Chr, c.Pos, c.LeftMean, c.RightMean, c.Score)
	}
	return nil
}

func writeQTLs(path string, qtls []QTL) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "CHR\tSTART\tEND\tPEAK_POS\tPEAK_SCORE\tMETHOD")
	for _, q := range qtls {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%.4f\t%s\n", q.Chr, q.Start, q.End, q.PeakPos, q.PeakScore, q.Method)
	}
	return nil
}

func buildSeries(snps []BSAstats, label, color string, pick func(BSAstats) float64) plotSeries {
	vals := make([]float64, len(snps))
	for i, s := range snps {
		vals[i] = pick(s)
	}
	return plotSeries{Label: label, Color: color, Values: vals}
}

func writePlots(path string, byChr map[string][]BSAstats, blocks []Block, scans []ScanWindow, changes []ChangePoint) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	chroms := make([]string, 0, len(byChr))
	for chr := range byChr {
		chroms = append(chroms, chr)
	}
	sort.Strings(chroms)

	fmt.Fprintln(f, "<!doctype html><html><head><meta charset=\"utf-8\"><title>GoBSAseq report</title><style>body{font-family:Arial,sans-serif;background:#f8fafc;color:#0f172a;margin:24px;} .section{margin:24px 0 40px;} .plot{background:#fff;border:1px solid #e2e8f0;border-radius:8px;padding:12px;box-shadow:0 1px 2px rgba(15,23,42,.05);margin-bottom:12px;} .plot-title{font-weight:600;margin-bottom:8px;} table{border-collapse:collapse;width:100%;background:#fff;} th,td{border:1px solid #e2e8f0;padding:6px 8px;text-align:left;font-size:12px;} th{background:#f1f5f9;}</style></head><body><h1>GoBSAseq two-bulk report</h1>")

	for _, chr := range chroms {
		snps := byChr[chr]
		fmt.Fprintf(f, "<div class=\"section\"><h2>%s</h2>", html.EscapeString(chr))

		deltaPlot := []plotSeries{
			buildSeries(snps, "Delta", "#94a3b8", func(s BSAstats) float64 { return s.DeltaSI }),
			buildSeries(snps, "Smoothed Delta", "#2563eb", func(s BSAstats) float64 { return s.DeltaSIK }),
			buildSeries(snps, "Segment Mean", "#7c3aed", func(s BSAstats) float64 { return s.SegMean }),
		}
		fmt.Fprintln(f, drawPlot("Delta-oriented signals", snps, deltaPlot, filterBlocksByChr(blocks, chr), nil, filterChangesByChr(changes, chr)))

		evidencePlot := []plotSeries{
			buildSeries(snps, "G-prime", "#dc2626", func(s BSAstats) float64 { return s.Gprime }),
			buildSeries(snps, "BB smooth", "#0f766e", func(s BSAstats) float64 { return s.BBK }),
			buildSeries(snps, "Scan LLR", "#d97706", func(s BSAstats) float64 { return s.ScanLLR }),
		}
		fmt.Fprintln(f, drawPlot("Evidence tracks", snps, evidencePlot, filterBlocksByChr(blocks, chr), filterScansByChr(scans, chr), filterChangesByChr(changes, chr)))

		controlPlot := []plotSeries{
			buildSeries(snps, "ED4 smooth", "#059669", func(s BSAstats) float64 { return s.EDK }),
			buildSeries(snps, "EWMA", "#0891b2", func(s BSAstats) float64 { return s.EWMA }),
			buildSeries(snps, "CUSUM", "#ea580c", func(s BSAstats) float64 { return s.CUSUM }),
		}
		fmt.Fprintln(f, drawPlot("Control and distance tracks", snps, controlPlot, filterBlocksByChr(blocks, chr), nil, filterChangesByChr(changes, chr)))
		fmt.Fprintln(f, "</div>")
	}
	fmt.Fprintln(f, "</body></html>")
	return nil
}

func drawPlot(title string, snps []BSAstats, series []plotSeries, blocks []Block, windows []ScanWindow, changes []ChangePoint) string {
	const (
		width  = 1100
		height = 280
		margin = 48.0
	)
	if len(snps) == 0 {
		return ""
	}

	minX, maxX := snps[0].POS, snps[len(snps)-1].POS
	minY, maxY := math.Inf(1), math.Inf(-1)
	for _, s := range series {
		for _, v := range s.Values {
			if !math.IsNaN(v) && !math.IsInf(v, 0) {
				if v < minY {
					minY = v
				}
				if v > maxY {
					maxY = v
				}
			}
		}
	}
	if math.IsInf(minY, 1) {
		minY, maxY = -1, 1
	}
	pad := (maxY - minY) * 0.1
	if pad == 0 {
		pad = 1
	}
	minY, maxY = minY-pad, maxY+pad

	scaleX := func(p int64) float64 { return margin + float64(p-minX)/float64(maxX-minX)*(width-2*margin) }
	scaleY := func(v float64) float64 { return height - margin - (v-minY)/(maxY-minY)*(height-2*margin) }

	var b strings.Builder
	fmt.Fprintf(&b, "<div class=\"plot\"><div class=\"plot-title\">%s</div><svg viewBox=\"0 0 %d %d\">", html.EscapeString(title), width, height)
	fmt.Fprintf(&b, "<rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#fff\"/>", width, height)

	for _, blk := range blocks {
		if blk.Sig {
			x := scaleX(blk.Start)
			w := scaleX(blk.End) - x
			fmt.Fprintf(&b, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#dcfce7\" opacity=\"0.4\"/>", x, margin, math.Max(1, w), height-2*margin)
		}
	}
	for i, w := range windows {
		if i < 5 {
			x := scaleX(w.Start)
			fmt.Fprintf(&b, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"#fde68a\" opacity=\"0.2\"/>", x, margin, math.Max(1, scaleX(w.End)-x), height-2*margin)
		}
	}
	for _, cp := range changes {
		x := scaleX(cp.Pos)
		fmt.Fprintf(&b, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"#7c3aed\" stroke-dasharray=\"5 4\" opacity=\"0.6\"/>", x, margin, x, height-margin)
	}

	for _, s := range series {
		var pts []string
		for i, v := range s.Values {
			if !math.IsNaN(v) && !math.IsInf(v, 0) {
				pts = append(pts, fmt.Sprintf("%.2f,%.2f", scaleX(snps[i].POS), scaleY(v)))
			}
		}
		fmt.Fprintf(&b, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.5\" points=\"%s\"/>", s.Color, strings.Join(pts, " "))
	}
	fmt.Fprintf(&b, "</svg></div>")
	return b.String()
}

// RunTwoBulkTwoParents is the entry point, now using the Analyzer struct
func RunTwoBulkTwoParents(vcfRdr *vcfgo.Reader, highPar int, highParDP int, lowPar int, lowParDP int, highBulk int, highBulkDP int, lowBulk int, lowBulkDP int, windowSize int) {
	a := NewAnalyzer(windowSize)
	a.HighPar, a.HighParDP = highPar, highParDP
	a.LowPar, a.LowParDP = lowPar, lowParDP
	a.HighBulk, a.HighBulkDP = highBulk, highBulkDP
	a.LowBulk, a.LowBulkDP = lowBulk, lowBulkDP
	a.Run(vcfRdr)
}
