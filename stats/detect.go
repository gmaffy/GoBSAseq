package stats

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"github.com/gmaffy/GoBSAseq/utils"
)

type QTL struct {
	CHROM string
	START float64
	STOP  float64
	PEAK  float64
}

func detectPeaks(chrom string, x []int64, stats []float64, threshData []float64) ([]QTL, error) {
	if len(x) != len(stats) || len(x) != len(threshData) {
		return nil, fmt.Errorf("detectPeaks: mismatched lengths x=%d stats=%d threshData=%d", len(x), len(stats), len(threshData))
	}
	if len(x) == 0 {
		return nil, nil
	}

	var qtls []QTL
	startIdx := -1
	if stats[0] > threshData[0] {
		startIdx = 0 // already inside a peak region at the start
	}

	for i := 1; i < len(x); i++ {
		switch {
		case stats[i-1] < threshData[i-1] && stats[i] > threshData[i]:
			startIdx = i

		case stats[i-1] > threshData[i-1] && stats[i] < threshData[i] && startIdx != -1:
			stopIdx := i
			peakIdx := startIdx
			for j := startIdx; j <= stopIdx; j++ {
				if stats[j] > stats[peakIdx] {
					peakIdx = j
				}
			}
			qtls = append(qtls, QTL{
				CHROM: chrom,
				START: float64(x[startIdx]),
				STOP:  float64(x[stopIdx]),
				PEAK:  float64(x[peakIdx]),
			})
			startIdx = -1
		}
	}

	if startIdx != -1 {
		stopIdx := len(x) - 1
		peakIdx := startIdx
		for j := startIdx; j <= stopIdx; j++ {
			if stats[j] > stats[peakIdx] {
				peakIdx = j
			}
		}
		qtls = append(qtls, QTL{
			CHROM: chrom,
			START: float64(x[startIdx]),
			STOP:  float64(x[stopIdx]),
			PEAK:  float64(x[peakIdx]),
		})
	}

	return qtls, nil
}

func detectTroughs(chrom string, x []int64, stats []float64, threshData []float64) ([]QTL, bool) {
	if len(x) != len(stats) || len(x) != len(threshData) {
		return nil, false
	}
	if len(x) == 0 {
		return nil, false
	}

	var qtls []QTL
	startIdx := -1
	if stats[0] < threshData[0] {
		startIdx = 0 // already inside a trough region at the start
	}

	for i := 1; i < len(x); i++ {
		switch {
		case stats[i-1] > threshData[i-1] && stats[i] < threshData[i]:
			startIdx = i

		case stats[i-1] < threshData[i-1] && stats[i] > threshData[i] && startIdx != -1:
			stopIdx := i
			peakIdx := startIdx
			for j := startIdx; j <= stopIdx; j++ {
				if stats[j] < stats[peakIdx] {
					peakIdx = j
				}
			}
			qtls = append(qtls, QTL{
				CHROM: chrom,
				START: float64(x[startIdx]),
				STOP:  float64(x[stopIdx]),
				PEAK:  float64(x[peakIdx]),
			})
			startIdx = -1
		}
	}

	if startIdx != -1 {
		stopIdx := len(x) - 1
		peakIdx := startIdx
		for j := startIdx; j <= stopIdx; j++ {
			if stats[j] < stats[peakIdx] {
				peakIdx = j
			}
		}
		qtls = append(qtls, QTL{
			CHROM: chrom,
			START: float64(x[startIdx]),
			STOP:  float64(x[stopIdx]),
			PEAK:  float64(x[peakIdx]),
		})
	}

	return qtls, len(qtls) > 0
}

func DetectIndividualStatQTLs(smoothed []SmoothedStats, thresholds []Thresholds, bsaType string, cfg *utils.AnalysisConfig) ([]QTL, error) {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.individual_stats_qtls.tsv", bsaType))

	// Group smoothed stats by chromosome
	chromToStats := make(map[string][]SmoothedStats)
	for _, s := range smoothed {
		chromToStats[s.CHROM] = append(chromToStats[s.CHROM], s)
	}

	// Build O(1) index: (CHROM, POS) → original slice index for threshold lookup.
	type chromPos struct{ chrom string; pos int64 }
	posToIdx := make(map[chromPos]int, len(smoothed))
	for i, s := range smoothed {
		posToIdx[chromPos{s.CHROM, s.POS}] = i
	}

	// Sorted chromosome list for deterministic output order.
	chroms := make([]string, 0, len(chromToStats))
	for c := range chromToStats {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// Process chromosomes concurrently.
	type chromResult struct {
		qtls []QTL
		err  error
	}
	results := make([]chromResult, len(chroms))
	var wg sync.WaitGroup
	wg.Add(len(chroms))

	for ci, chrom := range chroms {
		go func(ci int, chrom string) {
			defer wg.Done()

			chromStats := chromToStats[chrom]
			sort.Slice(chromStats, func(i, j int) bool { return chromStats[i].POS < chromStats[j].POS })

			// Build per-chromosome thresholds via the O(1) map.
			chromThresholds := make([]Thresholds, len(chromStats))
			for i, s := range chromStats {
				chromThresholds[i] = thresholds[posToIdx[chromPos{s.CHROM, s.POS}]]
			}

			positions := make([]int64, len(chromStats))
			var qtls []QTL

			// Two-bulk stats
			if hasBothBulks {
				for _, sf := range []struct {
					name      string
					getVal    func(SmoothedStats) float64
					getThresh func(Thresholds) float64
				}{
					{"DeltaSI", func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP95 }},
					{"Gstat", func(s SmoothedStats) float64 { return s.SmGstat }, func(t Thresholds) float64 { return t.TwoBulk.GstatP95 }},
					{"ED", func(s SmoothedStats) float64 { return s.SmED }, func(t Thresholds) float64 { return t.TwoBulk.ED4P95 }},
					{"LOD", func(s SmoothedStats) float64 { return s.SmLOD }, func(t Thresholds) float64 { return t.TwoBulk.LODP95 }},
					{"BBLogBF", func(s SmoothedStats) float64 { return s.SmBBLogBF }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP95 }},
				} {
					for i, s := range chromStats {
						positions[i] = s.POS
					}
					thresh := make([]float64, len(chromStats))
					for i, t := range chromThresholds {
						thresh[i] = sf.getThresh(t)
					}
					statVals := make([]float64, len(chromStats))
					for i, s := range chromStats {
						statVals[i] = sf.getVal(s)
					}
					q, err := detectPeaks(chrom, positions, statVals, thresh)
					if err != nil {
						results[ci] = chromResult{err: err}
						return
					}
					qtls = append(qtls, q...)
				}
			}

			// One-bulk stats
			if hasOneBulk {
				for _, sf := range []struct {
					name      string
					getVal    func(SmoothedStats) float64
					getThresh func(Thresholds) float64
				}{
					{"AFDev", func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevP95 }},
					{"OneBulkG", func(s SmoothedStats) float64 { return s.SmOneBulkG }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP95 }},
					{"OneBulkLOD", func(s SmoothedStats) float64 { return s.SmOneBulkLOD }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP95 }},
					{"OneBulkBBLogBF", func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP95 }},
				} {
					for i, s := range chromStats {
						positions[i] = s.POS
					}
					thresh := make([]float64, len(chromStats))
					for i, t := range chromThresholds {
						thresh[i] = sf.getThresh(t)
					}
					statVals := make([]float64, len(chromStats))
					for i, s := range chromStats {
						statVals[i] = sf.getVal(s)
					}
					q, err := detectPeaks(chrom, positions, statVals, thresh)
					if err != nil {
						results[ci] = chromResult{err: err}
						return
					}
					qtls = append(qtls, q...)
				}
			}

			results[ci] = chromResult{qtls: qtls}
		}(ci, chrom)
	}
	wg.Wait()

	// Collect results in chromosome order.
	var allQtls []QTL
	for _, r := range results {
		if r.err != nil {
			return nil, r.err
		}
		allQtls = append(allQtls, r.qtls...)
	}

	// Write results
	if err := writeQTLTSV(outPath, allQtls); err != nil {
		return nil, err
	}

	return allQtls, nil
}

func DetectCompositeZQTLs(smoothed []SmoothedStats, thresholds []Thresholds, bsaType string, cfg *utils.AnalysisConfig) ([]QTL, error) {
	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.composite.qtls.tsv", bsaType))

	var allQtls []QTL

	// Group smoothed stats by chromosome
	chromToStats := make(map[string][]SmoothedStats)
	for _, s := range smoothed {
		chromToStats[s.CHROM] = append(chromToStats[s.CHROM], s)
	}

	// Process each chromosome
	for chrom, stats := range chromToStats {
		// Sort by position
		sort.Slice(stats, func(i, j int) bool { return stats[i].POS < stats[j].POS })

		// Extract positions and CompositeZ values
		positions := make([]int64, len(stats))
		compositeZ := make([]float64, len(stats))
		
		// Get the thresholds for this chromosome
		var chromThresholds []Thresholds
		for _, s := range smoothed {
			if s.CHROM == chrom {
				// Find the index of this variant in the original smoothed slice
				for i, orig := range smoothed {
					if orig.CHROM == s.CHROM && orig.POS == s.POS {
						chromThresholds = append(chromThresholds, thresholds[i])
						break
					}
				}
			}
		}
		
		for i, s := range stats {
			positions[i] = s.POS
			compositeZ[i] = s.CompositeZ
		}

		// Use CompositeZ threshold
		thresh := make([]float64, len(stats))
		for i, t := range chromThresholds {
			thresh[i] = t.Z.CompositeZP95
		}

		// Detect peaks in CompositeZ
		qtls, err := detectPeaks(chrom, positions, compositeZ, thresh)
		if err != nil {
			return nil, err
		}
		allQtls = append(allQtls, qtls...)
	}

	// Write results
	if err := writeQTLTSV(outPath, allQtls); err != nil {
		return nil, err
	}

	return allQtls, nil
}

func FinalQTLs(smoothed []SmoothedStats, thresholds []Thresholds, bsaType string, cfg *utils.AnalysisConfig) ([]QTL, error) {
	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.qtls.tsv", bsaType))

	// Get QTLs from individual stats
	individualQtls, err := DetectIndividualStatQTLs(smoothed, thresholds, bsaType, cfg)
	if err != nil {
		return nil, err
	}

	// Get QTLs from composite Z
	compositeQtls, err := DetectCompositeZQTLs(smoothed, thresholds, bsaType, cfg)
	if err != nil {
		return nil, err
	}

	// Combine and deduplicate QTLs
	// Simple approach: merge overlapping QTLs and track supporting stats
	allQtls := append(individualQtls, compositeQtls...)

	// Sort by chromosome and start position
	sort.Slice(allQtls, func(i, j int) bool {
		if allQtls[i].CHROM != allQtls[j].CHROM {
			return allQtls[i].CHROM < allQtls[j].CHROM
		}
		return allQtls[i].START < allQtls[j].START
	})

	// Merge overlapping QTLs
	var finalQtls []QTL
	if len(allQtls) > 0 {
		current := allQtls[0]
		for i := 1; i < len(allQtls); i++ {
			next := allQtls[i]
			// If overlapping or within merge distance, merge them
			if next.CHROM == current.CHROM && next.START <= current.STOP {
				// Extend the current QTL
				if next.STOP > current.STOP {
					current.STOP = next.STOP
				}
				if next.PEAK > current.PEAK {
					current.PEAK = next.PEAK
				}
			} else {
				finalQtls = append(finalQtls, current)
				current = next
			}
		}
		finalQtls = append(finalQtls, current)
	}

	// Write results
	if err := writeQTLTSV(outPath, finalQtls); err != nil {
		return nil, err
	}

	return finalQtls, nil
}

func writeQTLTSV(filename string, qtls []QTL) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	header := []string{"CHROM", "START", "STOP", "PEAK"}
	fmt.Fprintln(w, strings.Join(header, "\t"))

	for _, q := range qtls {
		row := []string{
			q.CHROM,
			fmt.Sprintf("%.0f", q.START),
			fmt.Sprintf("%.0f", q.STOP),
			fmt.Sprintf("%.0f", q.PEAK),
		}
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}
	return nil
}
