package stats

import (
	"math"
)

// BRMBlock represents a BRM-style contiguous interval where a statistic is
// consistently beyond a variance-derived threshold.
type BRMBlock struct {
	Chrom     string
	Start     int64
	Stop      int64
	PeakPos   int64
	Peak      float64
	Threshold float64
	// ExpectedSI is used for one-bulk BRM blocks (can be zero for two-bulk)
	ExpectedSI float64
}

const afpFloor = 0.05

// CalculateBRMBlocks implements the two-bulk BRM block detection algorithm.
// - stats: smoothed windows (BSASmoothStats)
// - highBulkSize, lowBulkSize: integer sample sizes for variance scaling
// - popLevel: population level (e.g., 0 for F2, 1 for backcross generation scaling as 2^popLevel)
// - uAlpha: z-quantile (e.g., NormalQuantile(1 - alpha/2))
func CalculateBRMBlocks(chrom string, stats []BSASmoothStats, highBulkSize, lowBulkSize, popLevel int, uAlpha float64) []BRMBlock {
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
		start := stats[startIdx].WindowStart
		if startIdx > 0 {
			// midpoint between previous window end and this window start
			start = (stats[startIdx-1].WindowEnd + stats[startIdx].WindowStart) / 2
		}
		stop := stats[stopIdx].WindowEnd
		if stopIdx < len(stats)-1 {
			stop = (stats[stopIdx].WindowEnd + stats[stopIdx+1].WindowStart) / 2
		}
		if stop < start {
			stop = start
		}
		blocks = append(blocks, BRMBlock{
			Chrom:     chrom,
			Start:     start,
			Stop:      stop,
			PeakPos:   stats[peakIdx].WindowCenter,
			Peak:      peak,
			Threshold: threshold,
		})
	}

	for i, s := range stats {
		afp := (s.HighSI + s.LowSI) / 2
		if afp < afpFloor {
			afp = afpFloor
		}
		if afp > 1-afpFloor {
			afp = 1 - afpFloor
		}
		threshold := uAlpha * math.Sqrt(varianceScale*afp*(1-afp))
		significant := threshold > 0 && math.Abs(s.DeltaSI) >= threshold

		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = s.DeltaSI
				peakThreshold = threshold
				continue
			}
			if math.Abs(s.DeltaSI) > math.Abs(peak) {
				peakIdx = i
				peak = s.DeltaSI
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

// CalculateBRMBlocksOneBulk applies the one-bulk BRM threshold. expectedSI is
// the expected allele frequency under the null (e.g., 0.5 for F2).
func CalculateBRMBlocksOneBulk(chrom string, stats []BSASmoothStats, bulkSize, popLevel int, expectedSI, uAlpha float64) []BRMBlock {
	if len(stats) == 0 || bulkSize <= 0 || uAlpha <= 0 {
		return nil
	}

	p0 := expectedSI
	if math.IsNaN(p0) || math.IsInf(p0, 0) || p0 <= 0 || p0 >= 1 {
		p0 = 0.5
	}
	if p0 < afpFloor {
		p0 = afpFloor
	}
	if p0 > 1-afpFloor {
		p0 = 1 - afpFloor
	}

	n := float64(bulkSize)
	popScale := math.Pow(2, float64(popLevel))
	threshold := uAlpha * math.Sqrt((p0*(1-p0))/(popScale*n))
	if threshold <= 0 || math.IsNaN(threshold) || math.IsInf(threshold, 0) {
		return nil
	}

	var blocks []BRMBlock
	inBlock := false
	startIdx := 0
	peakIdx := 0
	peak := 0.0

	emitBlock := func(startIdx, stopIdx, peakIdx int, peak float64) {
		start := stats[startIdx].WindowStart
		if startIdx > 0 {
			start = (stats[startIdx-1].WindowEnd + stats[startIdx].WindowStart) / 2
		}
		stop := stats[stopIdx].WindowEnd
		if stopIdx < len(stats)-1 {
			stop = (stats[stopIdx].WindowEnd + stats[stopIdx+1].WindowStart) / 2
		}
		if stop < start {
			stop = start
		}
		blocks = append(blocks, BRMBlock{
			Chrom:      chrom,
			Start:      start,
			Stop:       stop,
			PeakPos:    stats[peakIdx].WindowCenter,
			Peak:       peak,
			ExpectedSI: p0,
			Threshold:  threshold,
		})
	}

	for i, s := range stats {
		// One-bulk smoothed deviation already stored in OneBulkAFDev (smoothed deviation)
		deviation := s.OneBulkAFDev
		significant := math.Abs(deviation) >= threshold

		if significant {
			if !inBlock {
				inBlock = true
				startIdx = i
				peakIdx = i
				peak = deviation
				continue
			}
			if math.Abs(deviation) > math.Abs(peak) {
				peakIdx = i
				peak = deviation
			}
			continue
		}

		if inBlock {
			emitBlock(startIdx, i-1, peakIdx, peak)
			inBlock = false
		}
	}

	if inBlock {
		emitBlock(startIdx, len(stats)-1, peakIdx, peak)
	}
	return blocks
}
