package twobulk

import (
	"math"
	"runtime"
	"sync"

	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/common"
)

// ============================================================================
// STATISTICAL FUNCTIONS
// ============================================================================

func gStatistic(alt1, ref1, alt2, ref2 int) float64 {
	n1 := float64(alt1 + ref1)
	n2 := float64(alt2 + ref2)
	total := n1 + n2

	if n1 == 0 || n2 == 0 || total == 0 {
		return 0
	}

	expAlt1 := n1 * float64(alt1+alt2) / total
	expRef1 := n1 * float64(ref1+ref2) / total
	expAlt2 := n2 * float64(alt1+alt2) / total
	expRef2 := n2 * float64(ref1+ref2) / total

	g := 0.0
	if alt1 > 0 && expAlt1 > 0 {
		g += float64(alt1) * math.Log(float64(alt1)/expAlt1)
	}
	if ref1 > 0 && expRef1 > 0 {
		g += float64(ref1) * math.Log(float64(ref1)/expRef1)
	}
	if alt2 > 0 && expAlt2 > 0 {
		g += float64(alt2) * math.Log(float64(alt2)/expAlt2)
	}
	if ref2 > 0 && expRef2 > 0 {
		g += float64(ref2) * math.Log(float64(ref2)/expRef2)
	}
	return 2 * g
}

func euclideanDistance(alt1, ref1, alt2, ref2 int) float64 {
	total1 := float64(alt1 + ref1)
	total2 := float64(alt2 + ref2)
	if total1 == 0 || total2 == 0 {
		return 0
	}
	p1 := float64(alt1) / total1
	p2 := float64(alt2) / total2
	return math.Sqrt(2 * math.Pow(p1-p2, 2))
}

func lodScore(ref1, alt1, ref2, alt2 int) float64 {
	n1 := float64(ref1 + alt1)
	n2 := float64(ref2 + alt2)
	if n1 == 0 || n2 == 0 {
		return 0
	}

	p1 := float64(alt1) / n1
	p2 := float64(alt2) / n2
	p0 := float64(alt1+alt2) / (n1 + n2)

	logLik := func(k, n float64, p float64) float64 {
		if n == 0 {
			return 0
		}
		if p <= 0 || p >= 1 {
			if (p <= 0 && k == 0) || (p >= 1 && k == n) {
				return 0
			}
			return -1e10
		}
		return k*math.Log(p) + (n-k)*math.Log(1-p)
	}

	l1 := logLik(float64(alt1), n1, p1) + logLik(float64(alt2), n2, p2)
	l0 := logLik(float64(alt1), n1, p0) + logLik(float64(alt2), n2, p0)
	return (l1 - l0) / math.Log(10)
}

func betaBinomialLogBF(succ1, fail1, succ2, fail2 int) float64 {
	alphaH, betaH := 1.0, 1.0
	alphaL, betaL := 1.0, 1.0

	logBeta := func(a, b float64) float64 {
		la, _ := math.Lgamma(a)
		lb, _ := math.Lgamma(b)
		lab, _ := math.Lgamma(a + b)
		return la + lb - lab
	}

	logAlt := logBeta(alphaH+float64(succ1), betaH+float64(fail1)) - logBeta(alphaH, betaH)
	logAlt += logBeta(alphaL+float64(succ2), betaL+float64(fail2)) - logBeta(alphaL, betaL)

	logNull := logBeta(1.0+float64(succ1+succ2), 1.0+float64(fail1+fail2)) - logBeta(1.0, 1.0)

	return logAlt - logNull
}

type BSAdata struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string
	Type  string // SNP or INDEL

	// Parent genotypes
	HighParentGT []int
	LowParentGT  []int

	// Bulk data
	HighBulkGT  []int
	LowBulkGT   []int
	HighBulkRef int
	HighBulkAlt int
	LowBulkRef  int
	LowBulkAlt  int
	HighBulkDP  int
	LowBulkDP   int

	// Statistics
	HighSI  float64 // High bulk SNP index
	LowSI   float64 // Low bulk SNP index
	DeltaSI float64 // Difference in SNP index
	Gstat   float64 // G-statistic
	ED      float64 // Euclidean distance
	LOD     float64 // LOD score
	BBLogBF float64 // Beta-binomial Bayes factor

}

func determineVariantType(v *vcfgo.Variant) string {
	// Check if it's an INDEL (based on length difference)
	if len(v.Reference) != len(v.Alt()[0]) {
		return "INDEL"
	}
	return "SNP"
}

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

func RunTwoBulkTwoParents(cfg common.AnalysisConfig) error {
	//var mu sync.Mutex

	numWorkers := runtime.NumCPU()
	variantChan := make(chan *vcfgo.Variant, 10000)
	var wg sync.WaitGroup

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for v := range variantChan {
				if !GoodVariants(v, cfg.HighParentIdx, cfg.HighParentDepth, cfg.LowParentIdx, cfg.LowParentDepth, cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.LowBulkIdx, cfg.LowBulkDepth) {
					continue
				}

				// Parse bulk depths
				highDP := v.Samples[cfg.HighBulkIdx].DP
				lowDP := v.Samples[cfg.LowBulkIdx].DP

				highRef, _ := v.Samples[cfg.HighBulkIdx].RefDepth()
				highAltDeps, _ := v.Samples[cfg.HighBulkIdx].AltDepths()
				lowRef, _ := v.Samples[cfg.LowBulkIdx].RefDepth()
				lowAltDeps, _ := v.Samples[cfg.LowBulkIdx].AltDepths()

				if len(highAltDeps) == 0 || len(lowAltDeps) == 0 {
					continue
				}

				// Determine which allele corresponds to low parent (reference)
				lpGT := v.Samples[cfg.LowParentIdx].GT[0]
				hpGT := v.Samples[cfg.HighParentIdx].GT[0]

				var highAlt, highRefCnt, lowAlt, lowRefCnt int
				if hpGT == lpGT {
					// High parent has same allele as low parent
					highRefCnt, highAlt = highRef, highAltDeps[0]
					lowRefCnt, lowAlt = lowRef, lowAltDeps[0]
				} else {
					// High parent has alternate allele
					highRefCnt, highAlt = highAltDeps[0], highRef
					lowRefCnt, lowAlt = lowAltDeps[0], lowRef
				}

				variant := &BSAdata{
					Chrom:        v.Chromosome,
					Pos:          int64(v.Pos),
					Ref:          v.Reference,
					Alt:          v.Alt()[0],
					Type:         determineVariantType(v),
					HighParentGT: v.Samples[cfg.HighParentIdx].GT,
					LowParentGT:  v.Samples[cfg.LowParentIdx].GT,
					HighBulkGT:   v.Samples[cfg.HighBulkIdx].GT,
					LowBulkGT:    v.Samples[cfg.LowBulkIdx].GT,
					HighBulkRef:  highRefCnt,
					HighBulkAlt:  highAlt,
					LowBulkRef:   lowRefCnt,
					LowBulkAlt:   lowAlt,
					HighBulkDP:   highDP,
					LowBulkDP:    lowDP,
				}

				// Calculate statistics
				variant.HighSI = float64(highAlt) / float64(highDP)
				variant.LowSI = float64(lowAlt) / float64(lowDP)
				variant.DeltaSI = variant.HighSI - variant.LowSI
				variant.Gstat = gStatistic(highAlt, highRefCnt, lowAlt, lowRefCnt)
				variant.ED = euclideanDistance(highAlt, highRefCnt, lowAlt, lowRefCnt)
				variant.LOD = lodScore(highRefCnt, highAlt, lowRefCnt, lowAlt)
				variant.BBLogBF = betaBinomialLogBF(highAlt, highRefCnt, lowAlt, lowRefCnt)

			}
		}()
	}
	return nil
}
