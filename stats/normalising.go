package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

func NormalizeSmoothedStats(cfg utils.AnalysisConfig, bsaType string, smoothed []BSASmoothStats) error {
	hasHighBulk := strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLowBulk := strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	hasBothBulks := hasHighBulk && hasLowBulk
	hasOneBulk := hasHighBulk != hasLowBulk

	zScore := func(values []float64) []float64 {
		z := make([]float64, len(values))
		if len(values) == 0 {
			return z
		}

		sortedValues := append([]float64(nil), values...)
		sort.Float64s(sortedValues)
		median := sortedValues[len(sortedValues)/2]
		if len(sortedValues)%2 == 0 {
			median = (sortedValues[len(sortedValues)/2-1] + sortedValues[len(sortedValues)/2]) / 2
		}

		deviations := make([]float64, len(values))
		for i, value := range values {
			deviations[i] = math.Abs(value - median)
		}
		sort.Float64s(deviations)
		mad := deviations[len(deviations)/2]
		if len(deviations)%2 == 0 {
			mad = (deviations[len(deviations)/2-1] + deviations[len(deviations)/2]) / 2
		}

		if mad > 0 {
			for i, value := range values {
				z[i] = math.Round((0.6745*(value-median)/mad)*1e6) / 1e6
			}
			return z
		}

		var mean float64
		for _, value := range values {
			mean += value
		}
		mean /= float64(len(values))

		var variance float64
		for _, value := range values {
			d := value - mean
			variance += d * d
		}
		sd := math.Sqrt(variance / float64(len(values)))
		if sd == 0 {
			return z
		}
		for i, value := range values {
			z[i] = math.Round(((value-mean)/sd)*1e6) / 1e6
		}
		return z
	}

	highSIValues := make([]float64, len(smoothed))
	lowSIValues := make([]float64, len(smoothed))
	deltaSIValues := make([]float64, len(smoothed))
	gprimeValues := make([]float64, len(smoothed))
	edValues := make([]float64, len(smoothed))
	lodValues := make([]float64, len(smoothed))
	bbLogBFValues := make([]float64, len(smoothed))
	afDevValues := make([]float64, len(smoothed))
	oneBulkGprimeValues := make([]float64, len(smoothed))
	oneBulkLODValues := make([]float64, len(smoothed))
	oneBulkBBLogBFValues := make([]float64, len(smoothed))

	for i, s := range smoothed {
		highSIValues[i] = s.HighSI
		lowSIValues[i] = s.LowSI
		deltaSIValues[i] = s.DeltaSI
		gprimeValues[i] = s.Gprime
		edValues[i] = s.ED
		lodValues[i] = s.LOD
		bbLogBFValues[i] = s.BBLogBF
		afDevValues[i] = s.OneBulkAFDev
		oneBulkGprimeValues[i] = s.OneBulkGprime
		oneBulkLODValues[i] = s.OneBulkLOD
		oneBulkBBLogBFValues[i] = s.OneBulkBBLogBF
	}

	highSIZ := zScore(highSIValues)
	lowSIZ := zScore(lowSIValues)
	deltaSIZ := zScore(deltaSIValues)
	gprimeZ := zScore(gprimeValues)
	edZ := zScore(edValues)
	lodZ := zScore(lodValues)
	bbLogBFZ := zScore(bbLogBFValues)
	afDevZ := zScore(afDevValues)
	oneBulkGprimeZ := zScore(oneBulkGprimeValues)
	oneBulkLODZ := zScore(oneBulkLODValues)
	oneBulkBBLogBFZ := zScore(oneBulkBBLogBFValues)

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return err
	}

	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.smoothed_and_normalised.tsv", bsaType))
	f, err := os.Create(outPath)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	header := []string{"CHROM", "WindowStart", "WindowEnd", "WindowCenter", "NSites", "MeanDepth"}
	if hasHighBulk {
		header = append(header, "SmoothedHighSI", "HighSI_Z")
	}
	if hasLowBulk {
		header = append(header, "SmoothedLowSI", "LowSI_Z")
	}
	if hasBothBulks {
		header = append(header, "SmoothedDeltaSI", "DeltaSI_Z", "Gprime", "Gprime_Z", "SmoothedED4", "ED4_Z", "WindowLOD", "WindowLOD_Z", "WindowBBLogBF", "WindowBBLogBF_Z")
	}
	if hasOneBulk {
		header = append(header, "P0", "SmoothedAFDev", "AFDev_Z", "Gprime1", "Gprime1_Z", "WindowLOD1", "WindowLOD1_Z", "WindowBBLogBF1", "WindowBBLogBF1_Z")
	}
	fmt.Fprintln(w, strings.Join(header, "\t"))

	for i, s := range smoothed {
		row := []string{
			s.CHROM,
			fmt.Sprintf("%d", s.WindowStart),
			fmt.Sprintf("%d", s.WindowEnd),
			fmt.Sprintf("%d", s.WindowCenter),
			fmt.Sprintf("%d", s.NSites),
			fmt.Sprintf("%.6f", s.MeanDepth),
		}
		if hasHighBulk {
			row = append(row, fmt.Sprintf("%.6f", s.HighSI), fmt.Sprintf("%.6f", highSIZ[i]))
		}
		if hasLowBulk {
			row = append(row, fmt.Sprintf("%.6f", s.LowSI), fmt.Sprintf("%.6f", lowSIZ[i]))
		}
		if hasBothBulks {
			row = append(row,
				fmt.Sprintf("%.6f", s.DeltaSI), fmt.Sprintf("%.6f", deltaSIZ[i]),
				fmt.Sprintf("%.6f", s.Gprime), fmt.Sprintf("%.6f", gprimeZ[i]),
				fmt.Sprintf("%.6f", s.ED), fmt.Sprintf("%.6f", edZ[i]),
				fmt.Sprintf("%.6f", s.LOD), fmt.Sprintf("%.6f", lodZ[i]),
				fmt.Sprintf("%.6f", s.BBLogBF), fmt.Sprintf("%.6f", bbLogBFZ[i]),
			)
		}
		if hasOneBulk {
			row = append(row,
				fmt.Sprintf("%.6f", s.OneBulkP0),
				fmt.Sprintf("%.6f", s.OneBulkAFDev), fmt.Sprintf("%.6f", afDevZ[i]),
				fmt.Sprintf("%.6f", s.OneBulkGprime), fmt.Sprintf("%.6f", oneBulkGprimeZ[i]),
				fmt.Sprintf("%.6f", s.OneBulkLOD), fmt.Sprintf("%.6f", oneBulkLODZ[i]),
				fmt.Sprintf("%.6f", s.OneBulkBBLogBF), fmt.Sprintf("%.6f", oneBulkBBLogBFZ[i]),
			)
		}
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}

	color.Green("Smoothed and normalised stats written to %s (%d windows)", outPath, len(smoothed))
	return nil
}
