package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sync"
	"sync/atomic"

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

func RawStats(cfg utils.AnalysisConfig, bsaType string, idxs []int, passedVariants []*vcfgo.Variant) ([]BSAstats, error) {
	if len(idxs) == 0 {
		return nil, fmt.Errorf("no sample indices supplied for %s raw stats", bsaType)
	}

	highParentIdx, lowParentIdx, highBulkIdx, lowBulkIdx := -1, -1, -1, -1
	switch bsaType {
	case "2p2b":
		highParentIdx, lowParentIdx, highBulkIdx, lowBulkIdx = 0, 1, 2, 3
	case "2phb":
		highParentIdx, lowParentIdx, highBulkIdx = 0, 1, 2
	case "2plb":
		highParentIdx, lowParentIdx, lowBulkIdx = 0, 1, 2
	case "hp2b":
		highParentIdx, highBulkIdx, lowBulkIdx = 0, 1, 2
	case "lp2b":
		lowParentIdx, highBulkIdx, lowBulkIdx = 0, 1, 2
	case "hphb":
		highParentIdx, highBulkIdx = 0, 1
	case "hplb":
		highParentIdx, lowBulkIdx = 0, 1
	case "lphb":
		lowParentIdx, highBulkIdx = 0, 1
	case "lplb":
		lowParentIdx, lowBulkIdx = 0, 1
	case "2b":
		highBulkIdx, lowBulkIdx = 0, 1
	default:
		return nil, fmt.Errorf("unsupported bsaseq type %q", bsaType)
	}

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}

	color.Cyan("============================ Calculating Raw Statistics (%s) ============================\n\n", bsaType)

	results := make([]BSAstats, len(passedVariants))
	keep := make([]bool, len(passedVariants))
	bar := progressbar.Default(int64(len(passedVariants)), "Processing variants")

	var next atomic.Int64
	var wg sync.WaitGroup
	workers := runtime.GOMAXPROCS(0)
	if workers < 1 {
		workers = 1
	}

	wg.Add(workers)
	for range workers {
		go func() {
			defer wg.Done()
			for {
				i := int(next.Add(1)) - 1
				if i >= len(passedVariants) {
					return
				}

				v := passedVariants[i]
				realAltIdx := -1
				realAltCount := 0
				for altIdx, alt := range v.Alt() {
					if alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<') {
						continue
					}
					realAltIdx = altIdx
					realAltCount++
				}
				if realAltCount != 1 {
					_ = bar.Add(1)
					continue
				}

				s := BSAstats{
					CHROM: v.Chromosome,
					POS:   int64(v.Pos),
					REF:   v.Reference,
					ALT:   v.Alt()[realAltIdx],
				}

				if highParentIdx >= 0 && highParentIdx < len(v.Samples) && v.Samples[highParentIdx] != nil {
					s.HighParGT = v.Samples[highParentIdx].GT
				}
				if lowParentIdx >= 0 && lowParentIdx < len(v.Samples) && v.Samples[lowParentIdx] != nil {
					s.LowParGT = v.Samples[lowParentIdx].GT
				}
				if highBulkIdx >= 0 && highBulkIdx < len(v.Samples) && v.Samples[highBulkIdx] != nil {
					s.HighBulkGT = v.Samples[highBulkIdx].GT
				}
				if lowBulkIdx >= 0 && lowBulkIdx < len(v.Samples) && v.Samples[lowBulkIdx] != nil {
					s.LowBulkGT = v.Samples[lowBulkIdx].GT
				}

				highAllele := 0
				switch {
				case len(s.HighParGT) > 0 && s.HighParGT[0] >= 0:
					highAllele = s.HighParGT[0]
				case len(s.LowParGT) > 0 && s.LowParGT[0] > 0:
					highAllele = 0
				case len(s.LowParGT) > 0 && s.LowParGT[0] == 0:
					highAllele = realAltIdx + 1
				}

				var highTotal, lowTotal int
				if highBulkIdx >= 0 && highBulkIdx < len(v.Samples) && v.Samples[highBulkIdx] != nil {
					refDepth, _ := v.Samples[highBulkIdx].RefDepth()
					altDepths, _ := v.Samples[highBulkIdx].AltDepths()
					if len(altDepths) <= realAltIdx {
						_ = bar.Add(1)
						continue
					}
					altDepth := altDepths[realAltIdx]
					if highAllele == 0 {
						s.HighBulkH, s.HighBulkL = refDepth, altDepth
					} else {
						s.HighBulkH, s.HighBulkL = altDepth, refDepth
					}
					s.HighBulkAD = fmt.Sprintf("%d,%d", refDepth, altDepth)
					highTotal = s.HighBulkH + s.HighBulkL
					if highTotal > 0 {
						s.HighSI = math.Round((float64(s.HighBulkH)/float64(highTotal))*1e6) / 1e6
					}
				}

				if lowBulkIdx >= 0 && lowBulkIdx < len(v.Samples) && v.Samples[lowBulkIdx] != nil {
					refDepth, _ := v.Samples[lowBulkIdx].RefDepth()
					altDepths, _ := v.Samples[lowBulkIdx].AltDepths()
					if len(altDepths) <= realAltIdx {
						_ = bar.Add(1)
						continue
					}
					altDepth := altDepths[realAltIdx]
					if highAllele == 0 {
						s.LowBulkH, s.LowBulkL = refDepth, altDepth
					} else {
						s.LowBulkH, s.LowBulkL = altDepth, refDepth
					}
					s.LowBulkAD = fmt.Sprintf("%d,%d", refDepth, altDepth)
					lowTotal = s.LowBulkH + s.LowBulkL
					if lowTotal > 0 {
						s.LowSI = math.Round((float64(s.LowBulkH)/float64(lowTotal))*1e6) / 1e6
					}
				}

				if highTotal == 0 && lowTotal == 0 {
					_ = bar.Add(1)
					continue
				}

				switch {
				case highTotal > 0 && (lowTotal == 0 || highTotal < lowTotal):
					s.Depth = highTotal
				case lowTotal > 0:
					s.Depth = lowTotal
				}

				if highTotal > 0 && lowTotal > 0 {
					s.DeltaSI = math.Round((s.HighSI-s.LowSI)*1e6) / 1e6
					s.Gstat = math.Round(GStatistic(s.HighBulkH, s.HighBulkL, s.LowBulkH, s.LowBulkL)*1e6) / 1e6
					s.ED = math.Round(euclideanDistance4(s.HighSI, s.LowSI)*1e6) / 1e6
					s.LOD = math.Round(lod(s.HighBulkL, s.HighBulkH, s.LowBulkL, s.LowBulkH)*1e6) / 1e6
					s.BBLogBF = math.Round(betaBinomialLogBF(s.HighBulkH, s.HighBulkL, s.LowBulkH, s.LowBulkL)*1e6) / 1e6
				}

				if s.Depth > 0 {
					results[i] = s
					keep[i] = true
				}
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	_ = bar.Finish()

	stats := make([]BSAstats, 0, len(results))
	for i, s := range results {
		if keep[i] {
			stats = append(stats, s)
		}
	}

	rawPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.raw.tsv", bsaType))
	if err := writeRawTSV(rawPath, stats); err != nil {
		return nil, err
	}

	color.Green("Raw stats written to %s (%d variants)", rawPath, len(stats))
	return stats, nil
}

func writeRawTSV(filename string, stats []BSAstats) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	fmt.Fprintln(w, "CHROM\tPOS\tREF\tALT\tHighParGT\tLowParGT\tHighBulkGT\tHighBulkAD\tLowBulkGT\tLowBulkAD\tHighBulkL\tHighBulkH\tLowBulkL\tLowBulkH\tHighSI\tLowSI\tDeltaSI\tGstat\tED4\tLOD\tBBLogBF\tDepth")
	for _, s := range stats {
		fmt.Fprintf(w, "%s\t%d\t%s\t%s\t%v\t%v\t%v\t%s\t%v\t%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n",
			s.CHROM, s.POS, s.REF, s.ALT,
			s.HighParGT, s.LowParGT, s.HighBulkGT, s.HighBulkAD, s.LowBulkGT, s.LowBulkAD,
			s.HighBulkL, s.HighBulkH, s.LowBulkL, s.LowBulkH,
			s.HighSI, s.LowSI, s.DeltaSI, s.Gstat, s.ED, s.LOD, s.BBLogBF, s.Depth)
	}
	return nil
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
