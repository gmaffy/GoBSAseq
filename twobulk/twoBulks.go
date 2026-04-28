package twobulk

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
)

type BSAstats struct {
	CHROM     string
	POS       int64
	REF       string
	ALT       string
	HighParGT []int
	//HighParAD string
	LowParGT []int
	//LowParAD   string
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

	DeltaSIK float64
	Gprime   float64
	EDK      float64
	LODK     float64
	BBLogBFK float64
}

type PopulationType string

const (
	F2  PopulationType = "F2"
	F3  PopulationType = "F3"
	BC  PopulationType = "BC"
	RIL PopulationType = "RIL"
)

func isHomozygous(gt []int) bool {
	if len(gt) == 0 {
		return false
	}
	for _, a := range gt[1:] {
		if a != gt[0] {
			return false
		}
	}
	return true
}

func GStatistic(highBulkAlt, highBulkRef, lowBulkAlt, lowBulkRef int) float64 {
	// G-test of independence
	highBulkTotal := float64(highBulkAlt + highBulkRef)
	lowBulkTotal := float64(lowBulkAlt + lowBulkRef)
	total := highBulkTotal + lowBulkTotal

	if highBulkTotal == 0 || lowBulkTotal == 0 || total == 0 {
		return 0
	}

	// Expected values under null hypothesis
	expHighAlt := highBulkTotal * float64(highBulkAlt+lowBulkAlt) / total
	expHighRef := highBulkTotal * float64(highBulkRef+lowBulkRef) / total
	expLowAlt := lowBulkTotal * float64(highBulkAlt+lowBulkAlt) / total
	expLowRef := lowBulkTotal * float64(highBulkRef+lowBulkRef) / total

	g := 0.0
	if highBulkAlt > 0 && expHighAlt > 0 {
		g += float64(highBulkAlt) * math.Log(float64(highBulkAlt)/expHighAlt)
	}
	if highBulkRef > 0 && expHighRef > 0 {
		g += float64(highBulkRef) * math.Log(float64(highBulkRef)/expHighRef)
	}
	if lowBulkAlt > 0 && expLowAlt > 0 {
		g += float64(lowBulkAlt) * math.Log(float64(lowBulkAlt)/expLowAlt)
	}
	if lowBulkRef > 0 && expLowRef > 0 {
		g += float64(lowBulkRef) * math.Log(float64(lowBulkRef)/expLowRef)
	}

	return 2 * g // Multiply by 2 for G-statistic
}

func EuclideanDist(refBulk1, altBulk1, refBulk2, altBulk2 int) float64 {
	// allele frequencies
	total1 := float64(refBulk1 + altBulk1)
	total2 := float64(refBulk2 + altBulk2)

	pRef := float64(refBulk1) / total1
	pAlt := float64(altBulk1) / total1
	qRef := float64(refBulk2) / total2
	qAlt := float64(altBulk2) / total2

	// Euclidean distance
	return math.Sqrt(math.Pow(pRef-qRef, 2) + math.Pow(pAlt-qAlt, 2))
}

func logBeta(a, b float64) float64 {
	la, _ := math.Lgamma(a)
	lb, _ := math.Lgamma(b)
	lab, _ := math.Lgamma(a + b)
	return la + lb - lab
}

// LOD: log10 likelihood ratio (bulks differ vs pooled)
func lod(ref1, alt1, ref2, alt2 int) float64 {
	n1 := float64(ref1 + alt1)
	n2 := float64(ref2 + alt2)
	total := n1 + n2

	p1 := float64(alt1) / n1
	p2 := float64(alt2) / n2
	p0 := float64(alt1+alt2) / total

	// Likelihoods
	L0 := math.Pow(p0, float64(alt1+alt2)) *
		math.Pow(1-p0, float64(ref1+ref2))
	L1 := math.Pow(p1, float64(alt1)) *
		math.Pow(1-p1, float64(ref1)) *
		math.Pow(p2, float64(alt2)) *
		math.Pow(1-p2, float64(ref2))

	if L0 == 0 || L1 == 0 {
		return 0.0
	}
	return math.Log10(L1 / L0)
}

func betaBinomialLogBF(highSucc, highFail, lowSucc, lowFail, N_high, N_low int) float64 {
	if N_high <= 0 {
		N_high = 1
	}
	if N_low <= 0 {
		N_low = 1
	}

	// Jeffreys prior scaled by pool size
	alphaH := 0.5 * float64(N_high)
	betaH := 0.5 * float64(N_high)
	alphaL := 0.5 * float64(N_low)
	betaL := 0.5 * float64(N_low)

	// H1: bulks differ
	logAlt := logBeta(alphaH+float64(highSucc), betaH+float64(highFail)) - logBeta(alphaH, betaH)
	logAlt += logBeta(alphaL+float64(lowSucc), betaL+float64(lowFail)) - logBeta(alphaL, betaL)

	// H0: bulks same
	alpha0 := 0.5 * float64(N_high+N_low) / 2.0
	beta0 := 0.5 * float64(N_high+N_low) / 2.0
	logNull := logBeta(alpha0+float64(highSucc+lowSucc), beta0+float64(highFail+lowFail)) - logBeta(alpha0, beta0)

	return logAlt - logNull
}

func biweight(x float64) float64 {
	if x >= 1 {
		return 0
	}
	t := 1 - x*x
	return t * t
}

func tricube(x float64) float64 {
	if x >= 1 {
		return 0
	}
	return math.Pow(1-math.Pow(x, 3), 3)
}

func smoothChromosomeBSA(stats []BSAstats, bandwidth int) {
	n := len(stats)

	for i := 0; i < n; i++ {
		center := stats[i].POS

		var numDelta, numED, numLOD, numBB, numG float64
		var denB, denG float64

		// look left
		for j := i; j >= 0; j-- {
			d := float64(center - stats[j].POS)
			if d > float64(bandwidth) {
				break
			}
			x := d / float64(bandwidth)

			wB := biweight(x)
			wG := tricube(x)

			numDelta += wB * stats[j].DeltaSI
			numED += wB * stats[j].ED
			numLOD += wB * stats[j].LOD
			numBB += wB * stats[j].BBLogBF
			numG += wG * stats[j].Gstat

			denB += wB
			denG += wG
		}

		// look right
		for j := i + 1; j < n; j++ {
			d := float64(stats[j].POS - center)
			if d > float64(bandwidth) {
				break
			}
			x := d / float64(bandwidth)

			wB := biweight(x)
			wG := tricube(x)

			numDelta += wB * stats[j].DeltaSI
			numED += wB * stats[j].ED
			numLOD += wB * stats[j].LOD
			numBB += wB * stats[j].BBLogBF
			numG += wG * stats[j].Gstat

			denB += wB
			denG += wG
		}

		if denB > 0 {
			stats[i].DeltaSIK = numDelta / denB
			stats[i].EDK = numED / denB
			stats[i].LODK = numLOD / denB
			stats[i].BBLogBFK = numBB / denB
		}

		if denG > 0 {
			stats[i].Gprime = numG / denG
		}
	}
}

func expectedAF(pop PopulationType, bcAltIsRecurrent bool) float64 {
	switch pop {
	case F2, F3, RIL:
		return 0.5
	case BC:
		if bcAltIsRecurrent {
			return 0.25
		}
		return 0.75
	default:
		return 0.5
	}
}

func simulateAltCount(depth int, p float64) int {
	c := 0
	for i := 0; i < depth; i++ {
		if rand.Float64() < p {
			c++
		}
	}
	return c
}

type SimMax struct {
	MaxDelta float64
	MaxED    float64
	MaxG     float64
	MaxLOD   float64
	MaxBB    float64
}

type Thresholds struct {
	Chrom string
	Delta float64
	ED    float64
	G     float64
	LOD   float64
	BB    float64
}

func simulateChromosomeMax(
	obs []BSAstats,
	pop PopulationType,
	bcAltIsRecurrent bool,
	windowSize int,
) SimMax {

	p0 := expectedAF(pop, bcAltIsRecurrent)

	sim := make([]BSAstats, len(obs))

	for i, s := range obs {
		hbDP := s.HighBulkL + s.HighBulkH
		lbDP := s.LowBulkL + s.LowBulkH

		altH := simulateAltCount(hbDP, p0)
		altL := simulateAltCount(lbDP, p0)

		refH := hbDP - altH
		refL := lbDP - altL

		hSI := float64(altH) / float64(hbDP)
		lSI := float64(altL) / float64(lbDP)

		sim[i] = BSAstats{
			CHROM:   s.CHROM,
			POS:     s.POS,
			DeltaSI: hSI - lSI,
			ED:      EuclideanDist(refH, altH, refL, altL),
			Gstat:   GStatistic(altH, refH, altL, refL),
			LOD:     lod(refH, altH, refL, altL),
			BBLogBF: betaBinomialLogBF(altH, refH, altL, refL, hbDP, lbDP),
		}
	}

	// Apply identical smoothing
	smoothChromosomeBSA(sim, windowSize)

	// Extract maxima
	var max SimMax
	for _, s := range sim {
		max.MaxDelta = math.Max(max.MaxDelta, math.Abs(s.DeltaSIK))
		max.MaxED = math.Max(max.MaxED, s.EDK)
		max.MaxG = math.Max(max.MaxG, s.Gprime)
		max.MaxLOD = math.Max(max.MaxLOD, s.LODK)
		max.MaxBB = math.Max(max.MaxBB, s.BBLogBFK)
	}
	return max
}

func estimateThresholds(
	obs []BSAstats,
	pop PopulationType,
	bcAltIsRecurrent bool,
	windowSize int,
	nSim int,
	alpha float64,
) Thresholds {

	var deltas, eds, gs, lods, bbs []float64

	for i := 0; i < nSim; i++ {
		m := simulateChromosomeMax(obs, pop, bcAltIsRecurrent, windowSize)
		deltas = append(deltas, m.MaxDelta)
		eds = append(eds, m.MaxED)
		gs = append(gs, m.MaxG)
		lods = append(lods, m.MaxLOD)
		bbs = append(bbs, m.MaxBB)
	}

	sort.Float64s(deltas)
	sort.Float64s(eds)
	sort.Float64s(gs)
	sort.Float64s(lods)
	sort.Float64s(bbs)

	q := int(float64(nSim) * (1.0 - alpha))

	return Thresholds{
		Chrom: obs[0].CHROM,
		Delta: deltas[q],
		ED:    eds[q],
		G:     gs[q],
		LOD:   lods[q],
		BB:    bbs[q],
	}
}

func GoodVariants(v *vcfgo.Variant, highPar int, highParDP int, lowPar int, lowParDP int, highBulk int, highBulkDP int, lowBulk int, lowBulkDP int) bool {
	indices := []int{highPar, lowPar, highBulk, lowBulk}
	if len(v.Alt()) != 1 {
		return false
	}

	for _, idx := range indices {
		s := v.Samples[idx]

		// ── Filter 1: no missing data (GT alleles must all be >= 0) ───────────
		if len(s.GT) == 0 {
			return false
		}
		for _, allele := range s.GT {
			if allele < 0 {
				return false
			}
		}
	}

	hpGT := v.Samples[highPar].GT
	lpGT := v.Samples[lowPar].GT

	hpDP := v.Samples[highPar].DP
	lpDP := v.Samples[lowPar].DP

	hbDP := v.Samples[highBulk].DP
	lbDP := v.Samples[lowBulk].DP

	if !isHomozygous(hpGT) || !isHomozygous(lpGT) {
		return false
	}
	// Homozygous means all alleles are the same, so compare the first allele
	if hpGT[0] == lpGT[0] {
		return false
	}

	if hpDP < highParDP || lpDP < lowParDP || hbDP < highBulkDP || lbDP <= lowBulkDP {
		return false
	}

	return true
}

func RunTwoBulkTwoParents(vcfRdr *vcfgo.Reader, highPar int, highParDP int, lowPar int, lowParDP int, highBulk int, highBulkDP int, lowBulk int, lowBulkDP int, windowSize int) {
	overallStart := time.Now()
	variantChan := make(chan *vcfgo.Variant, 1000)
	statsChan := make(chan BSAstats, 1000)

	// Worker pool
	numWorkers := runtime.NumCPU()
	var workerWG sync.WaitGroup
	var stats []BSAstats
	color.Cyan("============================ Calculating ∆SI, ED^4, LOD, and BBLogBF =============================\n\n")
	color.Cyan("Running %d workers .....\n", numWorkers)
	for i := 0; i < numWorkers; i++ {
		workerWG.Add(1)

		go func() {
			defer workerWG.Done()

			for variant := range variantChan {

				if !GoodVariants(
					variant,
					highPar, highParDP,
					lowPar, lowParDP,
					highBulk, highBulkDP,
					lowBulk, lowBulkDP,
				) {
					continue
				}

				lpGT := variant.Samples[lowPar].GT
				hbGT := variant.Samples[highBulk].GT
				lbGT := variant.Samples[lowBulk].GT

				hbDP := variant.Samples[highBulk].DP
				lbDP := variant.Samples[lowBulk].DP

				// Avoid divide-by-zero
				if hbDP == 0 || lbDP == 0 {
					continue
				}

				hbRefDep, _ := variant.Samples[highBulk].RefDepth()
				hbAltDeps, _ := variant.Samples[highBulk].AltDepths()

				lbRefDep, _ := variant.Samples[lowBulk].RefDepth()
				lbAltDeps, _ := variant.Samples[lowBulk].AltDepths()

				// Skip malformed multiallelic edge cases
				if len(hbAltDeps) == 0 || len(lbAltDeps) == 0 {
					continue
				}

				var hbL, hbH, lbL, lbH int

				if hbGT[0] == lpGT[0] {
					hbL = hbRefDep
					hbH = hbAltDeps[0]
				} else {
					hbL = hbAltDeps[0]
					hbH = hbRefDep
				}

				if lbGT[0] == lpGT[0] {
					lbL = lbRefDep
					lbH = lbAltDeps[0]
				} else {
					lbL = lbAltDeps[0]
					lbH = lbRefDep
				}

				hSI := float64(hbH) / float64(hbDP)
				lSI := float64(lbH) / float64(lbDP)

				deltaSI := hSI - lSI

				gstat := GStatistic(hbH, hbL, lbH, lbL)

				ed := EuclideanDist(hbL, hbH, lbL, lbH)

				lodVal := lod(hbL, hbH, lbL, lbH)

				bbLogBF := betaBinomialLogBF(hbH, hbL, lbH, lbL, hbDP, lbDP)

				stat := BSAstats{
					CHROM:     variant.Chromosome,
					POS:       int64(variant.Pos),
					REF:       variant.Reference,
					ALT:       variant.Alt()[0],
					HighParGT: variant.Samples[highPar].GT,
					//HighParAD: fmt.Sprintf("%s,%s", hbRefDep, ),
					LowParGT:   variant.Samples[lowPar].GT,
					HighBulkGT: variant.Samples[highBulk].GT,
					HighBulkAD: fmt.Sprintf("%v,%v", hbRefDep, hbAltDeps[0]),
					LowBulkGT:  variant.Samples[lowBulk].GT,
					LowBulkAD:  fmt.Sprintf("%v,%v", hbRefDep, hbAltDeps[0]),

					HighBulkL: hbL,
					HighBulkH: hbH,
					LowBulkL:  lbL,
					LowBulkH:  lbH,
					HighSI:    hSI,
					LowSI:     lSI,
					DeltaSI:   deltaSI,
					Gstat:     gstat,
					ED:        ed,
					LOD:       lodVal,
					BBLogBF:   bbLogBF,
				}
				//stats = append(stats, stat)
				statsChan <- stat
			}
		}()
	}
	color.New(color.FgHiWhite, color.Bold).Printf("\nBSAseq stats tool %s.\n", time.Since(overallStart).Round(time.Millisecond))

	smoothingStart := time.Now()
	color.Cyan("============================ Smoothing BSAseq stats =============================\n\n")
	// -------------------------------------------------- feeder -----------------------------------------------------
	go func() {
		for {
			v := vcfRdr.Read()
			if v == nil {
				break
			}
			variantChan <- v
		}
		close(variantChan)
	}()

	// 	--------------------------------------------------- collector ------------------------------------------------
	go func() {
		workerWG.Wait()
		close(statsChan)
	}()

	//var stats []BSAstats
	for s := range statsChan {
		stats = append(stats, s)
	}

	// ── SMOOTHING STAGE (NEW) ─────────────────
	byChr := make(map[string][]BSAstats)
	for _, s := range stats {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}

	for chr := range byChr {
		sort.Slice(byChr[chr], func(i, j int) bool {
			return byChr[chr][i].POS < byChr[chr][j].POS
		})
		smoothChromosomeBSA(byChr[chr], windowSize)
	}

	stats = stats[:0]
	for _, s := range byChr {
		stats = append(stats, s...)
	}

	color.New(color.FgHiWhite, color.Bold).Printf("\nBSAseq stats smoothing tool %s.\n", time.Since(smoothingStart).Round(time.Millisecond))

	thresholds := make(map[string]Thresholds)

	for chr, chrStats := range byChr {
		th := estimateThresholds(
			chrStats,
			F2,    // or F3, BC, RIL
			false, // BC config
			windowSize,
			5000, // simulations
			0.05, // alpha
		)
		thresholds[chr] = th
	}

}
