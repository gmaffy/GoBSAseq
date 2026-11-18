package stats

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// BRMBlock is a simple block record for BRM output.
type BRMBlock struct {
	Chrom     string
	Start     int64
	Stop      int64
	PeakPos   int64
	Peak      float64
	Threshold float64
	Support   SupportInterval
}

// BRM significance
//
// BRM tests the smoothed allele-frequency difference between bulks (ΔSI; for
// one bulk, the deviation from the expected frequency) against u·σ, where σ
// is the standard deviation of that difference caused by which plants ended
// up in each bulk (pool sampling, "drift"):
//
//	σ² = af(1−af)·[1/(V·n1) + 1/(V·n2)]   two bulks
//	σ² = p0(1−p0)/(V·n)                   one bulk
//
// with V the population variance scale (PopulationVarianceScale). Read
// sampling noise is not part of σ: the tested value is a Gaussian-smoothed
// average over many markers, which removes it almost entirely, and using a
// single marker's depth made the threshold jump from marker to marker.
//
// u is genome-wide: the null-genome simulation (nullsim.go) records, for each
// simulated genome, the largest |ΔSI|/σ anywhere, and u is its (1−α)
// quantile (α = --brm-alpha). So α is the probability that a genome with no
// QTL produces any BRM block, at the cross, bulk sizes, marker density and
// depths of the real data; whatever read noise survives smoothing is
// included through the simulation.

// brmSigmaTwoBulk is σ for two bulks at pooled allele frequency af.
func brmSigmaTwoBulk(af float64, n1, n2 int, popScale float64) float64 {
	af = math.Min(math.Max(af, 0.05), 0.95)
	v := af * (1 - af) / popScale
	return math.Sqrt(v/float64(n1) + v/float64(n2))
}

// brmSigmaOneBulk is σ for one bulk of n individuals with expected frequency p0.
func brmSigmaOneBulk(p0 float64, n int, popScale float64) float64 {
	if math.IsNaN(p0) || p0 <= 0 || p0 >= 1 {
		p0 = 0.5
	}
	p0 = math.Min(math.Max(p0, 0.05), 0.95)
	return math.Sqrt(p0 * (1 - p0) / (popScale * float64(n)))
}

// bulkSizes returns the bulk sizes used for a run: high and low bulk, and the
// single bulk for one-bulk designs (OneBulkSize, else the size of the bulk
// that is present). Shared by the null simulation and BRM so both use the
// same numbers.
func bulkSizes(cfg utils.AnalysisConfig, bsaType string) (high, low, one int) {
	hasHigh, _, _, hasOne := BulkFlags(bsaType)
	high, low, one = cfg.HighBulkSize, cfg.LowBulkSize, cfg.OneBulkSize
	if hasOne && one <= 0 {
		if hasHigh {
			one = high
		} else {
			one = low
		}
	}
	return high, low, one
}

// PopulationVarianceScale returns the exact scaling factor for segregation variance
// based on the mapping population structure. This factor adjusts the binomial
// variance p(1-p)/n to match the true variance Var(p) of the pooled individuals.
func PopulationVarianceScale(pop string) float64 {
	p := strings.ToUpper(strings.TrimSpace(pop))
	switch p {
	case "F2":
		return 2.0
	case "F3":
		return 4.0 / 3.0 // 1.333...
	case "F4":
		return 8.0 / 7.0 // 1.142...
	case "RIL":
		return 1.0
	}

	if strings.HasPrefix(p, "BC") {
		suffix := p[2:]
		var n int
		var err error

		if strings.HasSuffix(suffix, "F2H") || strings.HasSuffix(suffix, "F2L") {
			n, err = strconv.Atoi(suffix[:len(suffix)-3])
			if err == nil && n >= 1 {
				q2 := math.Pow(0.5, float64(n-1))
				return (4.0 - q2) / (3.0 - q2)
			}
		} else if strings.HasSuffix(suffix, "H") || strings.HasSuffix(suffix, "L") {
			n, err = strconv.Atoi(suffix[:len(suffix)-1])
			if err == nil && n >= 1 {
				q := math.Pow(0.5, float64(n))
				return 2.0 * (1.0 - 0.5*q) / (1.0 - q)
			}
		}
	}

	return 1.0 // Default to RIL-like (most conservative)
}

// calculateBRMBlocksTwoBulk applies BRM block detection for two-bulk smoothed
// windows: a marker is significant when |ΔSI| ≥ u·σ (see "BRM significance"),
// with σ = brmSigmaTwoBulk at the smoothed pooled frequency and u the
// genome-wide multiplier from the null simulation.
func calculateBRMBlocksTwoBulk(chrom string, stats []SmoothedStats, highBulkSize, lowBulkSize int, popScale, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || highBulkSize <= 0 || lowBulkSize <= 0 || uAlpha <= 0 {
		return nil
	}

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
		blocks = append(blocks, BRMBlock{Chrom: chrom, Start: start, Stop: stop, PeakPos: stats[peakIdx].POS, Peak: peak, Threshold: threshold})
	}

	for i, s := range stats {
		af := (s.SmHighSI + s.SmLowSI) / 2
		threshold := uAlpha * brmSigmaTwoBulk(af, highBulkSize, lowBulkSize, popScale)
		significant := threshold > 0 && math.Abs(s.SmDeltaSI) >= threshold

		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = s.SmDeltaSI
				peakThreshold = threshold
				continue
			}
			if math.Abs(s.SmDeltaSI) > math.Abs(peak) {
				peakIdx = i
				peak = s.SmDeltaSI
				peakThreshold = threshold
			}
			continue
		}

		if inBlock {
			stopIdx := i - 1
			emitBlock(startIdx, stopIdx, peakIdx, peak, peakThreshold)
			inBlock = false
		}
	}

	if inBlock {
		stopIdx := len(stats) - 1
		emitBlock(startIdx, stopIdx, peakIdx, peak, peakThreshold)
	}
	return blocks
}

// calculateBRMBlocksOneBulk is the one-bulk version: a marker is significant
// when |SI − p0| ≥ u·σ with σ = brmSigmaOneBulk (see "BRM significance").
func calculateBRMBlocksOneBulk(chrom string, stats []SmoothedStats, bulkSize int, popScale, expectedSI, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || bulkSize <= 0 || uAlpha <= 0 {
		return nil
	}
	threshold := uAlpha * brmSigmaOneBulk(expectedSI, bulkSize, popScale)

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
		blocks = append(blocks, BRMBlock{Chrom: chrom, Start: start, Stop: stop, PeakPos: stats[peakIdx].POS, Peak: peak, Threshold: threshold})
	}

	for i, s := range stats {
		if !(threshold > 0) {
			continue
		}
		// In one-bulk case SmoothedStats.SmAFDev holds AF deviation; but sm.SmHighSI or SmAFDev may be used.
		// Use SmAFDev when available (one-bulk pipeline sets it), else compute deviation from expected using SmHighSI.
		var deviation float64
		if s.SmAFDev != 0 {
			deviation = s.SmAFDev
		} else {
			deviation = s.SmHighSI - expectedSI
		}
		significant := math.Abs(deviation) >= threshold
		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = deviation
				peakThreshold = threshold
				continue
			}
			if math.Abs(deviation) > math.Abs(peak) {
				peakIdx = i
				peak = deviation
				peakThreshold = threshold
			}
			continue
		}
		if inBlock {
			stopIdx := i - 1
			emitBlock(startIdx, stopIdx, peakIdx, peak, peakThreshold)
			inBlock = false
		}
	}
	if inBlock {
		stopIdx := len(stats) - 1
		emitBlock(startIdx, stopIdx, peakIdx, peak, peakThreshold)
	}
	return blocks
}

// RunBRM runs BRM block detection across chromosomes and writes a TSV of
// blocks. u is the genome-wide multiplier from CalculateThresholds
// (Thresholds.BRM.U).
func RunBRM(cfg utils.AnalysisConfig, bsaType string, sm []SmoothedStats, u float64) ([]BRMBlock, error) {

	color.Cyan("===================================== Running BRM block detection =========================================")
	if len(sm) == 0 {
		return nil, nil
	}

	_, _, hasBoth, hasOne := BulkFlags(bsaType)
	popScale := PopulationVarianceScale(cfg.Population)
	if !(u > 0) {
		return nil, fmt.Errorf("BRM: genome-wide threshold multiplier not available (u=%v)", u)
	}
	uAlpha := u
	highSize, lowSize, oneSize := bulkSizes(cfg, bsaType)

	// group by chrom
	byChrom := make(map[string][]SmoothedStats)
	for _, s := range sm {
		byChrom[s.CHROM] = append(byChrom[s.CHROM], s)
	}

	lodFn := supportLOD(bsaType)
	var allBlocks []BRMBlock
	for chrom, stats := range byChrom {
		sort.Slice(stats, func(i, j int) bool { return stats[i].POS < stats[j].POS })
		var blocks []BRMBlock
		if hasBoth {
			blocks = calculateBRMBlocksTwoBulk(chrom, stats, highSize, lowSize, popScale, uAlpha)
		} else if hasOne {
			// expectedSI: use expectedHighAlleleP0 if available, else 0.5
			expectedSI, _ := ExpectedAF(cfg.Population)
			if expectedSI == 0 {
				expectedSI = 0.5
			}
			blocks = calculateBRMBlocksOneBulk(chrom, stats, oneSize, popScale, expectedSI, uAlpha)
		}
		if len(blocks) > 0 {
			merged := mergeBRMBlocks(blocks, cfg.EffectiveMergeDistance())
			for i := range merged {
				merged[i].Support = supportInChrom(stats, lodFn, float64(merged[i].Start), float64(merged[i].Stop), SupportLODDrop)
			}
			allBlocks = append(allBlocks, merged...)
		}
	}

	// write TSV
	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.brm_blocks.tsv", bsaType))
	f, err := os.Create(outPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	// SUPPORT_* (appended so existing column positions are unchanged) is the
	// LOD-drop interval that localises the QTL; START/STOP is the whole block.
	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK\tTHRESHOLD\tSUPPORT_START\tSUPPORT_END\tSUPPORT_PEAK_POS\tSUPPORT_PEAK_LOD")
	for _, b := range allBlocks {
		s := b.Support
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d\t%d\t%d\t%.6f\n", b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.Threshold, s.Start, s.End, s.PeakPos, s.PeakLOD)
	}
	color.Green("BRM blocks written to %s", outPath)
	return allBlocks, nil
}

func mergeBRMBlocks(blocks []BRMBlock, gap int64) []BRMBlock {
	if len(blocks) <= 1 {
		return blocks
	}
	sorted := append([]BRMBlock(nil), blocks...)
	sort.Slice(sorted, func(i, j int) bool {
		if sorted[i].Chrom != sorted[j].Chrom {
			return sorted[i].Chrom < sorted[j].Chrom
		}
		return sorted[i].Start < sorted[j].Start
	})

	out := make([]BRMBlock, 0, len(sorted))
	current := sorted[0]
	for i := 1; i < len(sorted); i++ {
		next := sorted[i]
		sameSign := current.Peak == 0 || next.Peak == 0 || (current.Peak > 0) == (next.Peak > 0)
		if next.Chrom == current.Chrom && sameSign && next.Start <= current.Stop+gap {
			if next.Stop > current.Stop {
				current.Stop = next.Stop
			}
			if math.Abs(next.Peak) > math.Abs(current.Peak) {
				current.Peak = next.Peak
				current.PeakPos = next.PeakPos
				current.Threshold = next.Threshold
			}
			continue
		}
		out = append(out, current)
		current = next
	}
	out = append(out, current)
	return out
}
