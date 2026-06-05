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

// Thresholds calculates per-chromosome thresholds for a selection of smoothed
// statistics (both on the smoothed scale and on the robust z-score / normalized
// scale). It writes a TSV file with one row per (chrom, stat, metric).
// The implementation keeps readability high and avoids unnecessary helpers.
func Thresholds(cfg utils.AnalysisConfig, bsaType string, smoothed []BSASmoothStats) error {
	hasHighBulk := strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLowBulk := strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	hasBothBulks := hasHighBulk && hasLowBulk
	hasOneBulk := hasHighBulk != hasLowBulk

	// group by chromosome
	byChrom := make(map[string][]BSASmoothStats)
	chroms := make([]string, 0)
	seen := make(map[string]bool)
	for _, s := range smoothed {
		byChrom[s.CHROM] = append(byChrom[s.CHROM], s)
		if !seen[s.CHROM] {
			chroms = append(chroms, s.CHROM)
			seen[s.CHROM] = true
		}
	}

	if len(chroms) == 0 {
		color.Yellow("No smoothed windows found; skipping threshold calculation")
		return nil
	}

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.thresholds.tsv", bsaType))
	f, err := os.Create(outPath)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()

	header := []string{"CHROM", "STAT", "SMOOTHED_P95", "SMOOTHED_P99", "SMOOTHED_P999", "NORML_Z3", "NORML_P95", "NORML_P99", "NORML_P999"}
	fmt.Fprintln(w, strings.Join(header, "\t"))

	// local helpers kept inside function for readability
	percentile := func(vals []float64, p float64) float64 {
		if len(vals) == 0 {
			return math.NaN()
		}
		s := append([]float64(nil), vals...)
		sort.Float64s(s)
		if p <= 0 {
			return s[0]
		}
		if p >= 1 {
			return s[len(s)-1]
		}
		idx := p*float64(len(s)-1)
		lo := int(math.Floor(idx))
		hi := int(math.Ceil(idx))
		if lo == hi {
			return s[lo]
		}
		fract := idx - float64(lo)
		return s[lo]*(1-fract) + s[hi]*fract
	}

	// robust z-score (median + MAD fallback to sd) — same logic as normalising.go
	robustZ := func(values []float64) []float64 {
		z := make([]float64, len(values))
		if len(values) == 0 {
			return z
		}

		sorted := append([]float64(nil), values...)
		sort.Float64s(sorted)
		median := sorted[len(sorted)/2]
		if len(sorted)%2 == 0 {
			median = (sorted[len(sorted)/2-1] + sorted[len(sorted)/2]) / 2
		}

		dev := make([]float64, len(values))
		for i, v := range values {
			dev[i] = math.Abs(v - median)
		}
		sort.Float64s(dev)
		mad := dev[len(dev)/2]
		if len(dev)%2 == 0 {
			mad = (dev[len(dev)/2-1] + dev[len(dev)/2]) / 2
		}

		if mad > 0 {
			for i, v := range values {
				z[i] = math.Round((0.6745*(v-median)/mad)*1e6) / 1e6
			}
			return z
		}

		// fallback to mean/std
		mean := 0.0
		for _, v := range values {
			mean += v
		}
		mean /= float64(len(values))
		variance := 0.0
		for _, v := range values {
			d := v - mean
			variance += d * d
		}
		sd := math.Sqrt(variance / float64(len(values)))
		if sd == 0 {
			return z
		}
		for i, v := range values {
			z[i] = math.Round(((v-mean)/sd)*1e6) / 1e6
		}
		return z
	}

	for _, chrom := range chroms {
		windows := byChrom[chrom]
		if len(windows) == 0 {
			continue
		}

		// collect arrays for each stat
		delta := make([]float64, 0, len(windows))
		gprime := make([]float64, 0, len(windows))
		ed := make([]float64, 0, len(windows))
		lod := make([]float64, 0, len(windows))
		bb := make([]float64, 0, len(windows))
		high := make([]float64, 0, len(windows))
		low := make([]float64, 0, len(windows))
		afdev := make([]float64, 0, len(windows))
		oneG := make([]float64, 0, len(windows))
		oneLOD := make([]float64, 0, len(windows))
		oneBB := make([]float64, 0, len(windows))

		for _, w := range windows {
			if hasBothBulks {
				delta = append(delta, w.DeltaSI)
				gprime = append(gprime, w.Gprime)
				ed = append(ed, w.ED)
				lod = append(lod, w.LOD)
				bb = append(bb, w.BBLogBF)
			}
			if hasHighBulk {
				high = append(high, w.HighSI)
			}
			if hasLowBulk {
				low = append(low, w.LowSI)
			}
			if hasOneBulk {
				afdev = append(afdev, w.OneBulkAFDev)
				oneG = append(oneG, w.OneBulkGprime)
				oneLOD = append(oneLOD, w.OneBulkLOD)
				oneBB = append(oneBB, w.OneBulkBBLogBF)
			}
		}

		// compute normalized z-scores per-chrom
		deltaZ := robustZ(delta)
		gZ := robustZ(gprime)
		edZ := robustZ(ed)
		lodZ := robustZ(lod)
		bbZ := robustZ(bb)
		highZ := robustZ(high)
		lowZ := robustZ(low)
		afdevZ := robustZ(afdev)
		oneGZ := robustZ(oneG)
		oneLODZ := robustZ(oneLOD)
		oneBBZ := robustZ(oneBB)

		// for each stat write thresholds: smoothed percentiles and normalized z thresholds
		writeRow := func(stat string, smoothedVals []float64, normZ []float64) {
			s95 := percentile(smoothedVals, 0.95)
			s99 := percentile(smoothedVals, 0.99)
			s999 := percentile(smoothedVals, 0.999)
			// normalized thresholds: fixed Z=3 and percentile thresholds on Z
			z3 := 3.0
			z95 := percentile(normZ, 0.95)
			z99 := percentile(normZ, 0.99)
			z999 := percentile(normZ, 0.999)
			fmt.Fprintf(w, "%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", chrom, stat, s95, s99, s999, z3, z95, z99, z999)
		}

		if hasBothBulks {
			writeRow("DeltaSI", delta, deltaZ)
			writeRow("Gprime", gprime, gZ)
			writeRow("ED4", ed, edZ)
			writeRow("LOD", lod, lodZ)
			writeRow("BBLogBF", bb, bbZ)
			// also write negative tail for DeltaSI (two-sided)
			// report negative percentiles as lower tail values
			if len(delta) > 0 {
				small95 := percentile(delta, 0.05)
				small99 := percentile(delta, 0.01)
				small999 := percentile(delta, 0.001)
				fmt.Fprintf(w, "%s\t%s_LO\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", chrom, "DeltaSI", small95, small99, small999, -3.0, -percentile(deltaZ, 0.05), -percentile(deltaZ, 0.01), -percentile(deltaZ, 0.001))
			}
		}
		if hasHighBulk {
			writeRow("HighSI", high, highZ)
			// negative tail for HighSI
			if len(high) > 0 {
				small95 := percentile(high, 0.05)
				small99 := percentile(high, 0.01)
				small999 := percentile(high, 0.001)
				fmt.Fprintf(w, "%s\t%s_LO\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", chrom, "HighSI", small95, small99, small999, -3.0, -percentile(highZ, 0.05), -percentile(highZ, 0.01), -percentile(highZ, 0.001))
			}
		}
		if hasLowBulk {
			writeRow("LowSI", low, lowZ)
			// negative tail for LowSI
			if len(low) > 0 {
				small95 := percentile(low, 0.05)
				small99 := percentile(low, 0.01)
				small999 := percentile(low, 0.001)
				fmt.Fprintf(w, "%s\t%s_LO\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", chrom, "LowSI", small95, small99, small999, -3.0, -percentile(lowZ, 0.05), -percentile(lowZ, 0.01), -percentile(lowZ, 0.001))
			}
		}
		if hasOneBulk {
			writeRow("AFDev", afdev, afdevZ)
			writeRow("Gprime1", oneG, oneGZ)
			writeRow("LOD1", oneLOD, oneLODZ)
			writeRow("BBLogBF1", oneBB, oneBBZ)
		}
	}

	color.Green("Thresholds written to %s", outPath)
	return nil
}
