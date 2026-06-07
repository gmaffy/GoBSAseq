package stats

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
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
}

// inverseNormalCDF implements an approximation of the inverse standard normal CDF
// using the algorithm by Peter J. Acklam. Accurate enough for typical alpha values.
func inverseNormalCDF(p float64) float64 {
	if p <= 0 || p >= 1 {
		if p <= 0 {
			return math.Inf(-1)
		}
		return math.Inf(1)
	}
	// Coefficients in rational approximations
	// Implementation adapted from https://web.archive.org/web/20150910023224/http://home.online.no/~pjacklam/notes/invnorm/
	var a = []float64{-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
		1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00}
	var b = []float64{-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
		6.680131188771972e+01, -1.328068155288572e+01}
	var c = []float64{-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
		-2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00}
	var d = []float64{7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00}

	// Define break-points.
	pLow := 0.02425
	pHigh := 1 - pLow
	var q, r, x float64
	if p < pLow {
		q = math.Sqrt(-2 * math.Log(p))
		x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q + c[5]) / (((d[0]*q+d[1])*q+d[2])*q + d[3])
		return -x
	}
	if p > pHigh {
		q = math.Sqrt(-2 * math.Log(1-p))
		x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q + c[5]) / (((d[0]*q+d[1])*q+d[2])*q + d[3])
		return x
	}
	q = p - 0.5
	r = q * q
	x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r + a[5]) * q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r + 1)
	return x
}

// popLevelFromPopulation returns the popLevel used by BRM variance scaling.
// Very small heuristic: F2 -> 1, else 0.
func popLevelFromPopulation(pop string) int {
	p := strings.ToUpper(strings.TrimSpace(pop))
	if p == "F2" {
		return 1
	}
	return 0
}

// calculateBRMBlocksTwoBulk applies BRM block detection for two-bulk smoothed windows.
// uAlpha should be the standard-normal quantile corresponding to 1 - alpha/2.
func calculateBRMBlocksTwoBulk(chrom string, stats []SmoothedStats, highBulkSize, lowBulkSize, popLevel int, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || highBulkSize <= 0 || lowBulkSize <= 0 || uAlpha <= 0 {
		return nil
	}

	n1 := float64(highBulkSize)
	n2 := float64(lowBulkSize)
	popScale := math.Pow(2, float64(popLevel))
	varianceScale := (n1 + n2) / (popScale * n1 * n2)

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
		if af < 0.05 {
			af = 0.05
		}
		if af > 0.95 {
			af = 0.95
		}
		threshold := uAlpha * math.Sqrt(varianceScale*af*(1-af))
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

// calculateBRMBlocksOneBulk mirrors the one-bulk BRM logic (testing SI vs expectedSI)
func calculateBRMBlocksOneBulk(chrom string, stats []SmoothedStats, bulkSize, popLevel int, expectedSI, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || bulkSize <= 0 || uAlpha <= 0 {
		return nil
	}
	p0 := expectedSI
	if math.IsNaN(p0) || p0 <= 0 || p0 >= 1 {
		p0 = 0.5
	}
	if p0 < 0.05 {
		p0 = 0.05
	}
	if p0 > 0.95 {
		p0 = 0.95
	}

	n := float64(bulkSize)
	popScale := math.Pow(2, float64(popLevel))
	threshold := uAlpha * math.Sqrt((p0*(1-p0))/(popScale*n))
	if !(threshold > 0) {
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

// RunBRM runs BRM block detection across chromosomes and writes a TSV of blocks.
func RunBRM(cfg utils.AnalysisConfig, bsaType string, sm []SmoothedStats) ([]BRMBlock, error) {
	if len(sm) == 0 {
		return nil, nil
	}

	_, _, hasBoth, hasOne := bulkFlags(bsaType)
	popLevel := popLevelFromPopulation(cfg.Population)
	uAlpha := inverseNormalCDF(1 - cfg.BrmAlpha/2)

	// group by chrom
	byChrom := make(map[string][]SmoothedStats)
	for _, s := range sm {
		byChrom[s.CHROM] = append(byChrom[s.CHROM], s)
	}

	var allBlocks []BRMBlock
	for chrom, stats := range byChrom {
		sort.Slice(stats, func(i, j int) bool { return stats[i].POS < stats[j].POS })
		var blocks []BRMBlock
		if hasBoth {
			blocks = calculateBRMBlocksTwoBulk(chrom, stats, cfg.HighBulkSize, cfg.LowBulkSize, popLevel, uAlpha)
		} else if hasOne {
			// expectedSI: use expectedHighAlleleP0 if available, else 0.5
			expectedSI, _ := expectedHighAlleleP0(cfg.Population)
			if expectedSI == 0 {
				expectedSI = 0.5
			}
			blocks = calculateBRMBlocksOneBulk(chrom, stats, cfg.OneBulkSize, popLevel, expectedSI, uAlpha)
		}
		if len(blocks) > 0 {
			allBlocks = append(allBlocks, blocks...)
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
	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK\tTHRESHOLD")
	for _, b := range allBlocks {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%.6f\t%.6f\n", b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.Threshold)
	}
	color.Green("BRM blocks written to %s", outPath)
	return allBlocks, nil
}
