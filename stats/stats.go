package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
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

	DeltaSI        float64
	Gstat          float64
	ED             float64
	LOD            float64
	BBLogBF        float64
	OneBulkP0      float64
	OneBulkAFDev   float64
	OneBulkGstat   float64
	OneBulkLOD     float64
	OneBulkBBLogBF float64

	Depth int
}

func safeGetDepths(s *vcfgo.SampleGenotype) (refDepth int, altDepths []int) {
	if s == nil || s.Fields == nil {
		return 0, nil
	}
	adStr, ok := s.Fields["AD"]
	if !ok || adStr == "" || adStr == "." {
		return 0, nil
	}
	parts := strings.Split(adStr, ",")
	if len(parts) > 0 {
		refDepth, _ = strconv.Atoi(parts[0])
		for i := 1; i < len(parts); i++ {
			d, _ := strconv.Atoi(parts[i])
			altDepths = append(altDepths, d)
		}
	}
	return
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

	hasOneBulk := (highBulkIdx >= 0) != (lowBulkIdx >= 0)
	oneBulkP0 := 0.0
	if hasOneBulk {
		p0, err := ExpectedAF(cfg.Population)
		if err != nil {
			return nil, err
		}
		oneBulkP0 = p0
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

				// Find the single real alt index (0-based). Variants reaching here
				// have already passed HardFilterVcf, but we still need the index to
				// look up allele depths and to record the ALT allele string.
				altIdx := -1
				for j, alt := range v.Alt() {
					if alt != "." && alt != "*" && !(len(alt) > 0 && alt[0] == '<') {
						altIdx = j
						break
					}
				}

				s := BSAstats{
					CHROM: v.Chromosome,
					POS:   int64(v.Pos),
					REF:   v.Reference,
					ALT:   v.Alt()[altIdx],
				}

				// Determine which allele index (0=ref, 1=alt) corresponds to the
				// high-parent phenotype. Defaults to ref (0) when no parent is present.
				highAllele := 0
				switch {
				case highParentIdx >= 0:
					highAllele = v.Samples[highParentIdx].GT[0]
				case lowParentIdx >= 0 && v.Samples[lowParentIdx].GT[0] == 0:
					highAllele = altIdx + 1
				}

				if highParentIdx >= 0 {
					s.HighParGT = v.Samples[highParentIdx].GT
				}
				if lowParentIdx >= 0 {
					s.LowParGT = v.Samples[lowParentIdx].GT
				}

				if highBulkIdx >= 0 {
					refDepth, altDepths := safeGetDepths(v.Samples[highBulkIdx])
					altDepth := 0
					if altIdx >= 0 && altIdx < len(altDepths) {
						altDepth = altDepths[altIdx]
					}
					s.HighBulkGT = v.Samples[highBulkIdx].GT
					s.HighBulkAD = fmt.Sprintf("%d,%d", refDepth, altDepth)
					if highAllele == 0 {
						s.HighBulkH, s.HighBulkL = refDepth, altDepth
					} else {
						s.HighBulkH, s.HighBulkL = altDepth, refDepth
					}
					highTotal := s.HighBulkH + s.HighBulkL
					if highTotal > 0 {
						s.HighSI = math.Round((float64(s.HighBulkH)/float64(highTotal))*1e6) / 1e6
					}
					s.Depth = highTotal
				}

				if lowBulkIdx >= 0 {
					refDepth, altDepths := safeGetDepths(v.Samples[lowBulkIdx])
					altDepth := 0
					if altIdx >= 0 && altIdx < len(altDepths) {
						altDepth = altDepths[altIdx]
					}
					s.LowBulkGT = v.Samples[lowBulkIdx].GT
					s.LowBulkAD = fmt.Sprintf("%d,%d", refDepth, altDepth)
					if highAllele == 0 {
						s.LowBulkH, s.LowBulkL = refDepth, altDepth
					} else {
						s.LowBulkH, s.LowBulkL = altDepth, refDepth
					}
					lowTotal := s.LowBulkH + s.LowBulkL
					if lowTotal > 0 {
						s.LowSI = math.Round((float64(s.LowBulkH)/float64(lowTotal))*1e6) / 1e6
					}
					if s.Depth == 0 || lowTotal < s.Depth {
						s.Depth = lowTotal
					}
				}

				if highBulkIdx >= 0 && lowBulkIdx >= 0 {
					highTotal := s.HighBulkH + s.HighBulkL
					lowTotal := s.LowBulkH + s.LowBulkL
					s.DeltaSI = math.Round((s.HighSI-s.LowSI)*1e6) / 1e6
					s.Gstat = math.Round(GStatistic(s.HighBulkH, s.HighBulkL, s.LowBulkH, s.LowBulkL)*1e6) / 1e6
					s.ED = math.Round(euclideanDistance4(s.HighSI, s.LowSI)*1e6) / 1e6
					s.LOD = math.Round(lod(s.HighBulkL, s.HighBulkH, s.LowBulkL, s.LowBulkH)*1e6) / 1e6
					s.BBLogBF = math.Round(betaBinomialLogBF(s.HighBulkH, s.HighBulkL, s.LowBulkH, s.LowBulkL)*1e6) / 1e6
					_ = highTotal
					_ = lowTotal
				}

				if hasOneBulk {
					s.OneBulkP0 = oneBulkP0
					h, l := s.HighBulkH, s.HighBulkL
					si := s.HighSI
					if lowBulkIdx >= 0 {
						h, l, si = s.LowBulkH, s.LowBulkL, s.LowSI
					}
					s.OneBulkAFDev = math.Round((si-oneBulkP0)*1e6) / 1e6
					s.OneBulkGstat = math.Round(oneBulkGStatistic(h, l, oneBulkP0)*1e6) / 1e6
					s.OneBulkLOD = math.Round(oneBulkLOD(h, l, oneBulkP0)*1e6) / 1e6
					s.OneBulkBBLogBF = math.Round(oneBulkBetaBinomialLogBF(h, l, oneBulkP0)*1e6) / 1e6
				}

				results[i] = s
				keep[i] = true
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
	if err := writeRawTSV(rawPath, stats, bsaType); err != nil {
		return nil, err
	}

	color.Green("Raw stats written to %s (%d variants)", rawPath, len(stats))
	return stats, nil
}

func writeRawTSV(filename string, stats []BSAstats, bsaType string) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	hasHighParent := strings.Contains(bsaType, "hp") || strings.Contains(bsaType, "2p")
	hasLowParent := strings.Contains(bsaType, "lp") || strings.Contains(bsaType, "2p")
	hasHighBulk := strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLowBulk := strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	hasBothBulks := hasHighBulk && hasLowBulk
	hasOneBulk := hasHighBulk != hasLowBulk

	header := []string{"CHROM", "POS", "REF", "ALT"}
	if hasHighParent {
		header = append(header, "HighParGT")
	}
	if hasLowParent {
		header = append(header, "LowParGT")
	}
	if hasHighBulk {
		header = append(header, "HighBulkGT", "HighBulkAD", "HighBulkL", "HighBulkH", "HighSI")
	}
	if hasLowBulk {
		header = append(header, "LowBulkGT", "LowBulkAD", "LowBulkL", "LowBulkH", "LowSI")
	}
	if hasBothBulks {
		header = append(header, "DeltaSI", "Gstat", "ED4", "LOD", "BBLogBF")
	}
	if hasOneBulk {
		header = append(header, "P0", "AFDev", "Gstat1", "LOD1", "BBLogBF1")
	}
	header = append(header, "Depth")

	fmt.Fprintln(w, strings.Join(header, "\t"))

	for _, s := range stats {
		row := []string{
			s.CHROM,
			fmt.Sprintf("%d", s.POS),
			s.REF,
			s.ALT,
		}
		if hasHighParent {
			row = append(row, fmt.Sprintf("%v", s.HighParGT))
		}
		if hasLowParent {
			row = append(row, fmt.Sprintf("%v", s.LowParGT))
		}
		if hasHighBulk {
			row = append(row, fmt.Sprintf("%v", s.HighBulkGT), s.HighBulkAD, fmt.Sprintf("%d", s.HighBulkL), fmt.Sprintf("%d", s.HighBulkH), fmt.Sprintf("%.6f", s.HighSI))
		}
		if hasLowBulk {
			row = append(row, fmt.Sprintf("%v", s.LowBulkGT), s.LowBulkAD, fmt.Sprintf("%d", s.LowBulkL), fmt.Sprintf("%d", s.LowBulkH), fmt.Sprintf("%.6f", s.LowSI))
		}
		if hasBothBulks {
			row = append(row, fmt.Sprintf("%.6f", s.DeltaSI), fmt.Sprintf("%.6f", s.Gstat), fmt.Sprintf("%.6f", s.ED), fmt.Sprintf("%.6f", s.LOD), fmt.Sprintf("%.6f", s.BBLogBF))
		}
		if hasOneBulk {
			row = append(row, fmt.Sprintf("%.6f", s.OneBulkP0), fmt.Sprintf("%.6f", s.OneBulkAFDev), fmt.Sprintf("%.6f", s.OneBulkGstat), fmt.Sprintf("%.6f", s.OneBulkLOD), fmt.Sprintf("%.6f", s.OneBulkBBLogBF))
		}
		row = append(row, fmt.Sprintf("%d", s.Depth))
		fmt.Fprintln(w, strings.Join(row, "\t"))
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

func ExpectedAF(popStruc string) (float64, error) {
	p := strings.ToUpper(strings.TrimSpace(popStruc))
	
	switch p {
	case "F2", "F3", "F4", "RIL":
		return 0.5, nil
	}

	if strings.HasPrefix(p, "BC") {
		suffix := p[2:]
		var n int
		var err error
		var isHigh bool
		
		if strings.HasSuffix(suffix, "F2H") || strings.HasSuffix(suffix, "F2L") {
			isHigh = strings.HasSuffix(suffix, "F2H")
			n, err = strconv.Atoi(suffix[:len(suffix)-3])
		} else if strings.HasSuffix(suffix, "H") || strings.HasSuffix(suffix, "L") {
			isHigh = strings.HasSuffix(suffix, "H")
			n, err = strconv.Atoi(suffix[:len(suffix)-1])
		}

		if err == nil && n >= 1 {
			p0 := 1.0 - math.Pow(0.5, float64(n+1))
			if isHigh {
				return p0, nil
			}
			return 1.0 - p0, nil
		}
	}
	
	return 0, fmt.Errorf("unknown population structure: %s", popStruc)
}

func oneBulkGStatistic(success, fail int, p0 float64) float64 {
	n := float64(success + fail)
	if n == 0 || p0 <= 0 || p0 >= 1 {
		return 0
	}

	obsSuccess := float64(success)
	obsFail := float64(fail)
	expSuccess := n * p0
	expFail := n * (1 - p0)

	g := 0.0
	if obsSuccess > 0 {
		g += obsSuccess * math.Log(obsSuccess/expSuccess)
	}
	if obsFail > 0 {
		g += obsFail * math.Log(obsFail/expFail)
	}
	return 2 * g
}

func oneBulkLOD(success, fail int, p0 float64) float64 {
	n := float64(success + fail)
	if n == 0 || p0 <= 0 || p0 >= 1 {
		return 0
	}

	phat := float64(success) / n
	logAlt := binomialLogLikelihood(float64(success), n, phat)
	logNull := binomialLogLikelihood(float64(success), n, p0)
	return (logAlt - logNull) / math.Log(10)
}

func oneBulkBetaBinomialLogBF(success, fail int, p0 float64) float64 {
	if success+fail == 0 || p0 <= 0 || p0 >= 1 {
		return 0
	}

	logAlt := logBeta(0.5+float64(success), 0.5+float64(fail)) - logBeta(0.5, 0.5)
	logNull := float64(success)*math.Log(p0) + float64(fail)*math.Log(1-p0)
	return (logAlt - logNull) / math.Log(10)
}

func binomialLogLikelihood(k, n, p float64) float64 {
	const eps = 1e-10
	if p <= 0 {
		p = eps
	}
	if p >= 1 {
		p = 1 - eps
	}
	return k*math.Log(p) + (n-k)*math.Log(1-p)
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
