package twobulk

import (
	"encoding/csv"
	"log"
	"math"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/brentp/vcfgo"
)

type BSAstats struct {
	CHROM      string
	POS        int64
	REF        string
	ALT        string
	HighParGT  []int
	LowParGT   []int
	HighBulkGT []int
	LowBulkGT  []int
	HighBulkL  int
	HighBulkH  int
	LowBulkL   int
	LowBulkH   int
	HighSI     float64
	LowSI      float64
	DeltaSI    float64
	Gstat      float64
	ED         float64
}

func GTToString(gt []int) string {
	if len(gt) == 0 {
		return "./."
	}

	parts := make([]string, len(gt))
	for i, allele := range gt {
		if allele < 0 {
			parts[i] = "."
		} else {
			parts[i] = strconv.Itoa(allele)
		}
	}

	return strings.Join(parts, "/")
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

func EuclideanDist(highSNP, lowSNP float64) float64 {
	// Euclidean distance between two bulk SNP indices
	diff := highSNP - lowSNP
	return math.Sqrt(diff * diff) // Could be extended to multi-dimensional
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

//func RunTwoBulkTwoParents(vcfRdr *vcfgo.Reader, highPar int, highParDP int, lowPar int, lowParDP int, highBulk int, highBulkDP int, lowBulk int, lowBulkDP int) {
//	statsCsv, err := os.Create("stats.csv")
//	if err != nil {
//		log.Fatal(err)
//	}
//	var stats []BSAstats
//	for variant := vcfRdr.Read(); variant != nil; variant = vcfRdr.Read() {
//		//fmt.Printf("Chromosome: %s, Pos:%v, Samples%s\n", variant.Chromosome, variant.Pos, variant.Samples[])
//		if GoodVariants(variant, highPar, highParDP, lowPar, lowParDP, highBulk, highBulkDP, lowBulk, lowBulkDP) {
//			//hpGT := variant.Samples[highPar].GT
//			lpGT := variant.Samples[lowPar].GT
//			hbGT := variant.Samples[highBulk].GT
//			lbGT := variant.Samples[lowBulk].GT
//			hbDP := variant.Samples[highBulk].DP
//			lbDP := variant.Samples[lowBulk].DP
//			hbRefDep, _ := variant.Samples[highBulk].RefDepth()
//			hbAltDeps, _ := variant.Samples[highBulk].AltDepths()
//			lbRefDep, _ := variant.Samples[lowBulk].RefDepth()
//			lbAltDeps, _ := variant.Samples[lowBulk].AltDepths()
//
//			var hbL int
//			var hbH int
//			var lbL int
//			var lbH int
//
//			if hbGT[0] == lpGT[0] {
//				hbL = hbRefDep
//				hbH = hbAltDeps[0]
//			} else {
//				hbL = hbAltDeps[0]
//				hbH = hbRefDep
//			}
//
//			if lbGT[0] == lpGT[0] {
//				lbL = lbRefDep
//				lbH = lbAltDeps[0]
//			} else {
//				lbL = lbAltDeps[0]
//				lbH = lbRefDep
//			}
//
//			hSI := float64(hbH / hbDP)
//			lSI := float64(lbH / lbDP)
//			deltaSI := hSI - lSI
//
//			gstat := GStatistic(hbAltDeps[0], hbRefDep, lbAltDeps[0], lbRefDep)
//			ed := EuclideanDist(float64(hbL), float64(lbL))
//
//			stat := BSAstats{
//				CHROM:      variant.Chromosome,
//				POS:        int64(variant.Pos),
//				REF:        variant.Reference,
//				ALT:        variant.Alt()[0],
//				HighParGT:  variant.Samples[highPar].GT,
//				LowParGT:   variant.Samples[lowPar].GT,
//				HighBulkGT: variant.Samples[highBulk].GT,
//				LowBulkGT:  variant.Samples[lowBulk].GT,
//				HighBulkL:  hbL,
//				HighBulkH:  hbH,
//				LowBulkL:   lbL,
//				LowBulkH:   lbH,
//				HighSI:     hSI,
//				LowSI:      lSI,
//				DeltaSI:    deltaSI,
//				Gstat:      gstat,
//				ED:         ed,
//			}
//
//
//			stats = append(stats, stat)
//
//		}
//
//	}
//}

func RunTwoBulkTwoParents(
	vcfRdr *vcfgo.Reader,
	highPar int, highParDP int,
	lowPar int, lowParDP int,
	highBulk int, highBulkDP int,
	lowBulk int, lowBulkDP int,
) {

	statsCsv, err := os.Create("stats.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer statsCsv.Close()

	csvWriter := csv.NewWriter(statsCsv)
	defer csvWriter.Flush()

	// CSV Header
	header := []string{
		"CHROM", "POS", "REF", "ALT",
		"HighParGT", "LowParGT", "HighBulkGT", "LowBulkGT",
		"HighBulkL", "HighBulkH", "LowBulkL", "LowBulkH",
		"HighSI", "LowSI", "DeltaSI", "Gstat", "ED",
	}
	if err := csvWriter.Write(header); err != nil {
		log.Fatal(err)
	}

	// Channels
	variantChan := make(chan *vcfgo.Variant, 1000)
	statsChan := make(chan BSAstats, 1000)

	// Writer goroutine
	var writerWG sync.WaitGroup
	writerWG.Add(1)

	go func() {
		defer writerWG.Done()

		for stat := range statsChan {
			record := []string{
				stat.CHROM,
				strconv.FormatInt(stat.POS, 10),
				stat.REF,
				stat.ALT,
				GTToString(stat.HighParGT),
				GTToString(stat.LowParGT),
				GTToString(stat.HighBulkGT),
				GTToString(stat.LowBulkGT),
				strconv.Itoa(stat.HighBulkL),
				strconv.Itoa(stat.HighBulkH),
				strconv.Itoa(stat.LowBulkL),
				strconv.Itoa(stat.LowBulkH),
				strconv.FormatFloat(stat.HighSI, 'f', 6, 64),
				strconv.FormatFloat(stat.LowSI, 'f', 6, 64),
				strconv.FormatFloat(stat.DeltaSI, 'f', 6, 64),
				strconv.FormatFloat(stat.Gstat, 'f', 6, 64),
				strconv.FormatFloat(stat.ED, 'f', 6, 64),
			}

			if err := csvWriter.Write(record); err != nil {
				log.Printf("CSV write error: %v", err)
			}
		}
	}()

	// Worker pool
	numWorkers := runtime.NumCPU()
	var workerWG sync.WaitGroup

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

				// IMPORTANT:
				// Your original code used integer division before float conversion.
				// This fixes that bug.
				hSI := float64(hbH) / float64(hbDP)
				lSI := float64(lbH) / float64(lbDP)

				deltaSI := hSI - lSI

				gstat := GStatistic(
					hbAltDeps[0],
					hbRefDep,
					lbAltDeps[0],
					lbRefDep,
				)

				ed := EuclideanDist(
					float64(hbL),
					float64(lbL),
				)

				stat := BSAstats{
					CHROM:      variant.Chromosome,
					POS:        int64(variant.Pos),
					REF:        variant.Reference,
					ALT:        variant.Alt()[0],
					HighParGT:  variant.Samples[highPar].GT,
					LowParGT:   variant.Samples[lowPar].GT,
					HighBulkGT: variant.Samples[highBulk].GT,
					LowBulkGT:  variant.Samples[lowBulk].GT,
					HighBulkL:  hbL,
					HighBulkH:  hbH,
					LowBulkL:   lbL,
					LowBulkH:   lbH,
					HighSI:     hSI,
					LowSI:      lSI,
					DeltaSI:    deltaSI,
					Gstat:      gstat,
					ED:         ed,
				}

				statsChan <- stat
			}
		}()
	}

	// Reader loop (single-threaded VCF stream)
	for variant := vcfRdr.Read(); variant != nil; variant = vcfRdr.Read() {
		variantChan <- variant
	}

	close(variantChan)

	workerWG.Wait()
	close(statsChan)

	writerWG.Wait()
}
