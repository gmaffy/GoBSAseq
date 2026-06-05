package stats

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

const minSNPsPerWindow = 5

type BSASmoothStats struct {
	CHROM        string
	WindowStart  int64
	WindowEnd    int64
	WindowCenter int64
	NSites       int
	MeanDepth    float64

	HighSI  float64
	LowSI   float64
	DeltaSI float64
	Gprime  float64
	ED      float64
	LOD     float64
	BBLogBF float64

	OneBulkP0      float64
	OneBulkAFDev   float64
	OneBulkGprime  float64
	OneBulkLOD     float64
	OneBulkBBLogBF float64
}

func SmoothStats(cfg utils.AnalysisConfig, bsaType string, rawStats []BSAstats) ([]BSASmoothStats, error) {
	if cfg.WindowSize <= 0 {
		return nil, fmt.Errorf("window size must be greater than 0")
	}
	if cfg.StepSize <= 0 {
		return nil, fmt.Errorf("step size must be greater than 0")
	}

	hasHighBulk := strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLowBulk := strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	hasBothBulks := hasHighBulk && hasLowBulk
	hasOneBulk := hasHighBulk != hasLowBulk

	byChrom := make(map[string][]BSAstats)
	chroms := make([]string, 0)
	seenChrom := make(map[string]bool)
	for _, s := range rawStats {
		byChrom[s.CHROM] = append(byChrom[s.CHROM], s)
		if !seenChrom[s.CHROM] {
			chroms = append(chroms, s.CHROM)
			seenChrom[s.CHROM] = true
		}
	}

	for _, chrom := range chroms {
		sort.Slice(byChrom[chrom], func(i, j int) bool {
			return byChrom[chrom][i].POS < byChrom[chrom][j].POS
		})
	}

	color.Cyan("============================ Smoothing Statistics (%s) ============================\n\n", bsaType)

	smoothed := make([]BSASmoothStats, 0)
	windowSize := int64(cfg.WindowSize)
	stepSize := int64(cfg.StepSize)
	halfWindow := float64(windowSize) / 2

	for _, chrom := range chroms {
		chromStats := byChrom[chrom]
		if len(chromStats) == 0 {
			continue
		}

		left := 0
		minPos := chromStats[0].POS
		maxPos := chromStats[len(chromStats)-1].POS
		for start := minPos; start <= maxPos; start += stepSize {
			end := start + windowSize - 1
			center := start + windowSize/2
			for left < len(chromStats) && chromStats[left].POS < start {
				left++
			}

			row := BSASmoothStats{
				CHROM:        chrom,
				WindowStart:  start,
				WindowEnd:    end,
				WindowCenter: center,
			}

			var depthWeightSum, tricubeWeightSum float64
			var highSISum, lowSISum, deltaSISum, edSum, afDevSum float64
			var depthSum, gprimeSum, lodSum, bbLogBFSum float64
			var oneBulkGprimeSum, oneBulkLODSum, oneBulkBBLogBFSum float64

			for i := left; i < len(chromStats) && chromStats[i].POS <= end; i++ {
				s := chromStats[i]
				distance := math.Abs(float64(s.POS - center))
				scaledDistance := distance / halfWindow
				if scaledDistance > 1 {
					continue
				}

				tricubeWeight := math.Pow(1-math.Pow(scaledDistance, 3), 3)
				depthWeight := tricubeWeight * math.Sqrt(float64(s.Depth))
				if tricubeWeight <= 0 || depthWeight <= 0 {
					continue
				}

				row.NSites++
				tricubeWeightSum += tricubeWeight
				depthWeightSum += depthWeight
				depthSum += tricubeWeight * float64(s.Depth)

				if hasHighBulk {
					highSISum += depthWeight * s.HighSI
				}
				if hasLowBulk {
					lowSISum += depthWeight * s.LowSI
				}
				if hasBothBulks {
					deltaSISum += depthWeight * s.DeltaSI
					edSum += depthWeight * s.ED
					gprimeSum += depthWeight * s.Gstat
					lodSum += tricubeWeight * s.LOD
					bbLogBFSum += tricubeWeight * s.BBLogBF
				}
				if hasOneBulk {
					row.OneBulkP0 = s.OneBulkP0
					afDevSum += depthWeight * s.OneBulkAFDev
					oneBulkGprimeSum += depthWeight * s.OneBulkGstat
					oneBulkLODSum += tricubeWeight * s.OneBulkLOD
					oneBulkBBLogBFSum += tricubeWeight * s.OneBulkBBLogBF
				}
			}

			if row.NSites < minSNPsPerWindow || tricubeWeightSum == 0 || depthWeightSum == 0 {
				continue
			}

			row.MeanDepth = math.Round((depthSum/tricubeWeightSum)*1e6) / 1e6
			if hasHighBulk {
				row.HighSI = math.Round((highSISum/depthWeightSum)*1e6) / 1e6
			}
			if hasLowBulk {
				row.LowSI = math.Round((lowSISum/depthWeightSum)*1e6) / 1e6
			}
			if hasBothBulks {
				row.DeltaSI = math.Round((deltaSISum/depthWeightSum)*1e6) / 1e6
				row.ED = math.Round((edSum/depthWeightSum)*1e6) / 1e6
				row.Gprime = math.Round((gprimeSum/depthWeightSum)*1e6) / 1e6
				row.LOD = math.Round(lodSum*1e6) / 1e6
				row.BBLogBF = math.Round(bbLogBFSum*1e6) / 1e6
			}
			if hasOneBulk {
				row.OneBulkAFDev = math.Round((afDevSum/depthWeightSum)*1e6) / 1e6
				row.OneBulkGprime = math.Round((oneBulkGprimeSum/depthWeightSum)*1e6) / 1e6
				row.OneBulkLOD = math.Round(oneBulkLODSum*1e6) / 1e6
				row.OneBulkBBLogBF = math.Round(oneBulkBBLogBFSum*1e6) / 1e6
			}
			smoothed = append(smoothed, row)
		}
	}

	color.Green("Smoothed %d windows", len(smoothed))
	return smoothed, nil
}
