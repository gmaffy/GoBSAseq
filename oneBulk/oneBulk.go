package oneBulk

import (
	"fmt"
	"math"

	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/utils"
)

type OneBulkStats struct {
	CHROM     string
	POS       int64
	REF       string
	ALT       string
	HighParGT []int
	LowParGT  []int
	BulkGT    []int
	BulkAD    string
	SI        float64 // ALT / (ALT+REF)
	AbsSI     float64 // |SI - 0.5|
	Gstat     float64 // one-bulk G vs. uniform
	LOD       float64 // one-bulk LOD vs. p=0.5
	BBLogBF   float64 // BF vs. p=0.5
	ED        float64 // |SI-0.5|^4
	Depth     int
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

func RunTwoParentsLowBulk(cfg utils.AnalysisConfig, hfcfg utils.HardFilterConfig) error {
	rdr := cfg.Rdr
	highParIdx := cfg.HighParentIdx
	lowParIdx := cfg.LowParentIdx
	lowBulkIdx := cfg.LowBulkIdx

	//-------------------------------------- Remove problematic fields ---------------------------------------------- //
	for _, id := range []string{"PGT", "PID"} {
		delete(rdr.Header.SampleFormats, id)
	}

	//-------------------------------------- Header for writing ----------------------------------------------------- //
	fmt.Println(rdr.Header.SampleNames[highParIdx], rdr.Header.SampleNames[lowParIdx], rdr.Header.SampleNames[lowBulkIdx], highParIdx, lowParIdx, lowBulkIdx)
	sampleNames := []string{rdr.Header.SampleNames[highParIdx], rdr.Header.SampleNames[lowParIdx], rdr.Header.SampleNames[lowBulkIdx]}

	writerHeader := *rdr.Header
	writerHeader.SampleNames = sampleNames

	// ------------------------------------------------- Run -------------------------------------------------------- //

	badVariant := 0
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
		//if err := rdr.Error(); err != nil {
		//	if strings.Contains(err.Error(), "bad sample string") {
		//		//fmt.Printf("Bad sample string: %s\n", v.String())
		//		//rdr.Clear()
		//		continue
		//	} else {
		//		//fmt.Printf("VCF parse error at line %d: %w", v.LineNumber, err)
		//		continue
		//	}
		//}

		passed := utils.PassesHardFilter(v, hfcfg) && v.Samples[cfg.HighParentIdx].DP >= cfg.HighParentDepth && v.Samples[cfg.LowParentIdx].DP >= cfg.LowParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
		if passed {
			highPar := v.Samples[cfg.HighParentIdx]
			lowPar := v.Samples[cfg.LowParentIdx]
			lowBulk := v.Samples[cfg.LowBulkIdx]

			bulkRefDep, _ := lowBulk.RefDepth()
			bulkAltDeps, _ := lowBulk.AltDepths()
			//fmt.Println(highPar.GT, lowPar.GT, lowBulk.GT, highPar.DP, lowPar.DP, lowBulk.DP)
			var bulkSusAlleleCount int
			//var bulkResAlleleCount int
			if lowBulk.GT[0] == lowPar.GT[0] {
				bulkSusAlleleCount = bulkRefDep
				//bulkResAlleleCount = bulkAltDeps[0]
			} else {
				bulkSusAlleleCount = bulkAltDeps[0]
				//bulkResAlleleCount = bulkRefDep
			}

			SI := float64(bulkSusAlleleCount) / float64(lowBulk.DP)
			fmt.Println(SI)
			s := OneBulkStats{
				CHROM:     v.Chromosome,
				POS:       int64(v.Pos),
				REF:       v.Reference,
				ALT:       v.Alt()[0],
				HighParGT: highPar.GT,
				LowParGT:  lowPar.GT,
				BulkGT:    lowBulk.GT,
				BulkAD:    fmt.Sprintf("%d,%d", bulkRefDep, bulkAltDeps[0]),
				SI:        SI,
				AbsSI:     math.Abs(SI - 0.5),
				ED:        math.Pow(math.Abs(SI-0.5), 4),

				Gstat:   math.Round(GStatisticOneBulk(lowBulk.GT)*1e6) / 1e6,
				LOD:     math.Round(LodOneBulk(lowBulk.GT)*1e6) / 1e6,
				BBLogBF: math.Round(BetaBinomialOneBulk(lowBulk.GT)*1e6) / 1e6,
				Depth:   cfg.LowBulkDepth,
			}
			v.Samples = []*vcfgo.SampleGenotype{highPar, lowPar, lowBulk}
		}
	}

	return nil

}
