package stats

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"

	"github.com/gmaffy/GoBSAseq/utils"
)

// ThresholdPoint describes a point on the smoothed-stat and threshold tracks.
type ThresholdPoint struct {
	POS       float64
	Stat      float64
	Threshold float64
}

// ThresholdPeak is one threshold-bounded high peak or low trough.
type ThresholdPeak struct {
	Chrom          string
	Start          ThresholdPoint
	Peak           ThresholdPoint
	Stop           ThresholdPoint
	AboveThreshold bool
}

// IndividualStatQTL is a merged interval supported by one or more individual
// smoothed statistics. Non-contributing peak fields are NaN.
type IndividualStatQTL struct {
	Chrom       string
	QTLStart    float64
	QTLStop     float64
	DeltaSIPeak float64
	AFDevPeak   float64
	GstatPeak   float64
	EDPeak      float64
	LODPeak     float64
	BBLogBFPeak float64
}

// FindIndividualStatQTLs finds per-stat threshold peaks, then merges any
// overlapping intervals on the same chromosome. thresholdTier accepts "P99" or
// "P95"; an empty value defaults to "P99".
func FindIndividualStatQTLs(
	smoothed []SmoothedStats,
	thresholds []Thresholds,
	bsaType string,
	thresholdTier string,
) ([]IndividualStatQTL, error) {
	if len(smoothed) != len(thresholds) {
		return nil, fmt.Errorf("smoothed and threshold slices must have the same length")
	}

	useP95 := false
	switch thresholdTier {
	case "", "P99", "p99", "99":
	case "P95", "p95", "95":
		useP95 = true
	default:
		return nil, fmt.Errorf("unsupported threshold tier %q", thresholdTier)
	}

	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	type statTrack struct {
		name             string
		significantAbove bool
		statValue        func(SmoothedStats) float64
		thresholdValue   func(Thresholds) float64
	}
	type statInterval struct {
		chrom       string
		start, stop float64
		stat        string
		peak        float64
	}

	var tracks []statTrack
	if hasBothBulks {
		tracks = append(tracks,
			statTrack{
				name:             "DeltaSI",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmDeltaSI },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.TwoBulk.DeltaSIP95
					}
					return t.TwoBulk.DeltaSIP99
				},
			},
			statTrack{
				name:             "DeltaSI",
				significantAbove: false,
				statValue:        func(s SmoothedStats) float64 { return s.SmDeltaSI },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.TwoBulk.DeltaSIMp95
					}
					return t.TwoBulk.DeltaSIMp99
				},
			},
			statTrack{
				name:             "Gstat",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmGstat },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.TwoBulk.GstatP95
					}
					return t.TwoBulk.GstatP99
				},
			},
			statTrack{
				name:             "ED",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmED },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.TwoBulk.ED4P95
					}
					return t.TwoBulk.ED4P99
				},
			},
			statTrack{
				name:             "LOD",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmLOD },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.TwoBulk.LODP95
					}
					return t.TwoBulk.LODP99
				},
			},
			statTrack{
				name:             "BBLogBF",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmBBLogBF },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.TwoBulk.BBLogBFP95
					}
					return t.TwoBulk.BBLogBFP99
				},
			},
		)
	}
	if hasOneBulk {
		tracks = append(tracks,
			statTrack{
				name:             "AFDev",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmAFDev },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.OneBulk.AFDevP95
					}
					return t.OneBulk.AFDevP99
				},
			},
			statTrack{
				name:             "AFDev",
				significantAbove: false,
				statValue:        func(s SmoothedStats) float64 { return s.SmAFDev },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.OneBulk.AFDevMp95
					}
					return t.OneBulk.AFDevMp99
				},
			},
			statTrack{
				name:             "Gstat",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmOneBulkG },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.OneBulk.OneBulkGstatP95
					}
					return t.OneBulk.OneBulkGstatP99
				},
			},
			statTrack{
				name:             "LOD",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmOneBulkLOD },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.OneBulk.OneBulkLODP95
					}
					return t.OneBulk.OneBulkLODP99
				},
			},
			statTrack{
				name:             "BBLogBF",
				significantAbove: true,
				statValue:        func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF },
				thresholdValue: func(t Thresholds) float64 {
					if useP95 {
						return t.OneBulk.OneBulkBBLogBFP95
					}
					return t.OneBulk.OneBulkBBLogBFP99
				},
			},
		)
	}
	if len(tracks) == 0 {
		return nil, fmt.Errorf("unsupported bsaseq type %q", bsaType)
	}

	var intervals []statInterval
	for _, track := range tracks {
		peaks, err := FindThresholdPeaks(smoothed, thresholds, track.statValue, track.thresholdValue)
		if err != nil {
			return nil, err
		}
		for _, peak := range peaks {
			if peak.AboveThreshold != track.significantAbove {
				continue
			}
			intervals = append(intervals, statInterval{
				chrom: peak.Chrom,
				start: peak.Start.POS,
				stop:  peak.Stop.POS,
				stat:  track.name,
				peak:  peak.Peak.Stat,
			})
		}
	}
	if len(intervals) == 0 {
		return nil, nil
	}

	sort.Slice(intervals, func(i, j int) bool {
		if intervals[i].chrom != intervals[j].chrom {
			return intervals[i].chrom < intervals[j].chrom
		}
		if intervals[i].start == intervals[j].start {
			return intervals[i].stop < intervals[j].stop
		}
		return intervals[i].start < intervals[j].start
	})

	applyPeak := func(qtl *IndividualStatQTL, interval statInterval) {
		setPeak := func(current *float64) {
			if math.IsNaN(*current) || math.Abs(interval.peak) > math.Abs(*current) {
				*current = interval.peak
			}
		}
		switch interval.stat {
		case "DeltaSI":
			setPeak(&qtl.DeltaSIPeak)
		case "AFDev":
			setPeak(&qtl.AFDevPeak)
		case "Gstat":
			setPeak(&qtl.GstatPeak)
		case "ED":
			setPeak(&qtl.EDPeak)
		case "LOD":
			setPeak(&qtl.LODPeak)
		case "BBLogBF":
			setPeak(&qtl.BBLogBFPeak)
		}
	}
	newQTL := func(interval statInterval) IndividualStatQTL {
		qtl := IndividualStatQTL{
			Chrom:       interval.chrom,
			QTLStart:    interval.start,
			QTLStop:     interval.stop,
			DeltaSIPeak: math.NaN(),
			AFDevPeak:   math.NaN(),
			GstatPeak:   math.NaN(),
			EDPeak:      math.NaN(),
			LODPeak:     math.NaN(),
			BBLogBFPeak: math.NaN(),
		}
		applyPeak(&qtl, interval)
		return qtl
	}

	qtls := []IndividualStatQTL{newQTL(intervals[0])}
	for _, interval := range intervals[1:] {
		current := &qtls[len(qtls)-1]
		if interval.chrom == current.Chrom && interval.start <= current.QTLStop {
			if interval.stop > current.QTLStop {
				current.QTLStop = interval.stop
			}
			applyPeak(current, interval)
			continue
		}
		qtls = append(qtls, newQTL(interval))
	}

	return qtls, nil
}

// DetectIndividualStatQTLs finds merged individual-stat QTLs and writes them to
// the standard individual-statistics QTL TSV.
func DetectIndividualStatQTLs(
	smoothed []SmoothedStats,
	thresholds []Thresholds,
	bsaType string,
	cfg *utils.AnalysisConfig,
) ([]IndividualStatQTL, error) {
	if cfg == nil {
		return nil, fmt.Errorf("analysis config must not be nil")
	}

	qtls, err := FindIndividualStatQTLs(smoothed, thresholds, bsaType, "P99")
	if err != nil {
		return nil, err
	}

	outDir := filepath.Join(cfg.OutputDir, "qtls")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, fmt.Errorf("DetectIndividualStatQTLs: mkdir %s: %w", outDir, err)
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.individual_stats_qtls.tsv", bsaType))

	f, err := os.Create(outPath)
	if err != nil {
		return nil, fmt.Errorf("DetectIndividualStatQTLs: create %s: %w", outPath, err)
	}
	defer f.Close()

	formatPeak := func(v float64) string {
		if math.IsNaN(v) {
			return "NA"
		}
		return fmt.Sprintf("%.6f", v)
	}

	fmt.Fprintln(f, "CHROM\tQTL_START\tQTL_STOP\tDELTA_SI_PEAK\tAFDEV_PEAK\tGSTAT_PEAK\tED_PEAK\tLOD_PEAK\tBBLOGBF_PEAK")
	for _, qtl := range qtls {
		fmt.Fprintf(f, "%s\t%.6f\t%.6f\t%s\t%s\t%s\t%s\t%s\t%s\n",
			qtl.Chrom,
			qtl.QTLStart,
			qtl.QTLStop,
			formatPeak(qtl.DeltaSIPeak),
			formatPeak(qtl.AFDevPeak),
			formatPeak(qtl.GstatPeak),
			formatPeak(qtl.EDPeak),
			formatPeak(qtl.LODPeak),
			formatPeak(qtl.BBLogBFPeak),
		)
	}

	return qtls, nil
}

// FindThresholdPeaks walks each chromosome from left to right and returns every
// threshold-bounded peak or trough between consecutive threshold crossings.
func FindThresholdPeaks(
	smoothed []SmoothedStats,
	thresholds []Thresholds,
	statValue func(SmoothedStats) float64,
	thresholdValue func(Thresholds) float64,
) ([]ThresholdPeak, error) {
	if len(smoothed) != len(thresholds) {
		return nil, fmt.Errorf("smoothed and threshold slices must have the same length")
	}
	if statValue == nil {
		return nil, fmt.Errorf("stat value function must not be nil")
	}
	if thresholdValue == nil {
		return nil, fmt.Errorf("threshold value function must not be nil")
	}
	if len(smoothed) == 0 {
		return nil, nil
	}

	valueAt := func(i int) (float64, float64) {
		return statValue(smoothed[i]), thresholdValue(thresholds[i])
	}
	isFinitePoint := func(i int) bool {
		stat, threshold := valueAt(i)
		return !math.IsNaN(stat) && !math.IsInf(stat, 0) &&
			!math.IsNaN(threshold) && !math.IsInf(threshold, 0)
	}
	sideAt := func(i int) int {
		if !isFinitePoint(i) {
			return 0
		}
		stat, threshold := valueAt(i)
		switch {
		case stat > threshold:
			return 1
		case stat < threshold:
			return -1
		default:
			return 0
		}
	}
	peakAt := func(i int) ThresholdPoint {
		stat, threshold := valueAt(i)
		return ThresholdPoint{
			POS:       float64(smoothed[i].POS),
			Stat:      stat,
			Threshold: threshold,
		}
	}
	intersection := func(left, right int) ThresholdPoint {
		leftStat, leftThreshold := valueAt(left)
		rightStat, rightThreshold := valueAt(right)
		leftDiff := leftStat - leftThreshold
		rightDiff := rightStat - rightThreshold

		fraction := 0.0
		if leftDiff != rightDiff {
			fraction = leftDiff / (leftDiff - rightDiff)
		}
		fraction = min(max(fraction, 0), 1)

		pos := float64(smoothed[left].POS) + fraction*float64(smoothed[right].POS-smoothed[left].POS)
		threshold := leftThreshold + fraction*(rightThreshold-leftThreshold)
		return ThresholdPoint{
			POS:       pos,
			Stat:      threshold,
			Threshold: threshold,
		}
	}

	var peaks []ThresholdPeak
	for start := 0; start < len(smoothed); {
		stop := start + 1
		for stop < len(smoothed) && smoothed[stop].CHROM == smoothed[start].CHROM {
			stop++
		}

		type activePeak struct {
			start          ThresholdPoint
			peakIndex      int
			aboveThreshold bool
		}

		var active *activePeak
		previousIndex := -1
		previousSide := 0

		for i := start; i < stop; i++ {
			currentSide := sideAt(i)
			if currentSide == 0 {
				continue
			}
			if previousSide == 0 {
				previousIndex = i
				previousSide = currentSide
				continue
			}

			if currentSide == previousSide {
				if active != nil {
					stat, peakStat := statValue(smoothed[i]), statValue(smoothed[active.peakIndex])
					if (active.aboveThreshold && stat > peakStat) || (!active.aboveThreshold && stat < peakStat) {
						active.peakIndex = i
					}
				}
				previousIndex = i
				previousSide = currentSide
				continue
			}

			crossing := intersection(previousIndex, i)
			if active != nil {
				peaks = append(peaks, ThresholdPeak{
					Chrom:          smoothed[active.peakIndex].CHROM,
					Start:          active.start,
					Peak:           peakAt(active.peakIndex),
					Stop:           crossing,
					AboveThreshold: active.aboveThreshold,
				})
			}
			active = &activePeak{
				start:          crossing,
				peakIndex:      i,
				aboveThreshold: currentSide > 0,
			}
			previousIndex = i
			previousSide = currentSide
		}

		start = stop
	}

	return peaks, nil
}
