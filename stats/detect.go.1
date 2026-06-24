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

type QTLRecord struct {
	Chrom     string
	Start     int64
	Stop      int64
	PeakPos   int64
	PeakVal   float64
	Stat      string // "Stouffer" or "MaxAbsZ"
	Threshold string // "P99" or "P95"
}

type DirectQTLRecord struct {
	Chrom      string
	Start      int64
	Stop       int64
	PeakPos    int64
	ZPeak      float64
	Statistics string
	Threshold  string
}

// MergedQTL is a unified QTL record that may originate from CompositeZ
// detection, BRM block detection, or both.  When both methods agree on a
// chromosome the intervals are unioned and both peak values are recorded.
type MergedQTL struct {
	Chrom string
	Start int64
	Stop  int64

	// CompositeZ peak (from DetectQTLs); math.NaN() when absent.
	ZPeak float64
	// BRM peak delta-SI / AF-deviation (from RunBRM); math.NaN() when absent.
	BRMPeak float64
	// BRM significance threshold at the peak position; math.NaN() when absent.
	BRMThreshold float64

	// Source summarises which method(s) contributed.
	// Values: "ZScore", "BRM", "ZScore+BRM"
	Source string
}

// geneSpaceRegion satisfies the GeneSpaceInterval interface in genespace.go.
func (m MergedQTL) geneSpaceRegion() (string, int, int) {
	return m.Chrom, int(m.Start), int(m.Stop)
}

// Z-score reference thresholds (identical to legacy values).
const (
	zSig = 3.0 // ~p99 equivalent — significant
	//zSugg = 2.0 // ~p95 equivalent — suggestive
)

//func DetectQTLsWithMCDirect(cfg utils.AnalysisConfig, bsaType string, smoothed []SmoothedStats, thresholds []Thresholds) ([]QTLRecord, error) {
//	color.Cyan("===================================== Detecting QTLs with MC Thresholds =========================================")
//	if len(smoothed) == 0 {
//		fmt.Println("No smoothed stats to detect QTLs from")
//		return nil, nil
//	}
//	if len(thresholds) != len(smoothed) {
//		return nil, fmt.Errorf("smoothed and thresholds slices must have the same length")
//	}
//
//	// ── group by chromosome ──────────────────────────────────────────────────
//
//	byChr := make(map[string][]SmoothedStats)
//	byChrThresholds := make(map[string][]Thresholds)
//	for i, s := range smoothed {
//		byChr[s.CHROM] = append(byChr[s.CHROM], s)
//		byChrThresholds[s.CHROM] = append(byChrThresholds[s.CHROM], thresholds[i])
//	}
//	chroms := make([]string, 0, len(byChr))
//	for c := range byChr {
//		chroms = append(chroms, c)
//	}
//	sort.Strings(chroms)
//
//	// ── scan each chromosome ─────────────────────────────────────────────────
//
//	var qtls []QTLRecord
//
//	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
//
//	for _, chrom := range chroms {
//		stats := byChr[chrom]
//		threshs := byChrThresholds[chrom]
//
//		// Collect every contiguous run where key statistics exceed their MC thresholds.
//		type run struct {
//			start, stop int64
//			peak        float64 // the most extreme statistic value at the peak
//		}
//
//		var runs []run
//		inRun := false
//		var cur run
//
//		for i, s := range stats {
//			t := threshs[i]
//
//			// Check which statistics are significant based on MC thresholds
//			// We use a combined approach: any of the main statistics exceeding threshold
//			isSignificant := false
//			var maxStat float64
//
//			if hasBothBulks {
//				// Check DeltaSI against MC threshold
//				if s.SmDeltaSI >= t.TwoBulk.DeltaSIP99 || s.SmDeltaSI <= t.TwoBulk.DeltaSIMp99 {
//					isSignificant = true
//					maxStat = math.Abs(s.SmDeltaSI)
//				}
//				// Check G-statistic
//				if !isSignificant && s.SmGstat >= t.TwoBulk.GstatP99 {
//					isSignificant = true
//					maxStat = s.SmGstat
//				}
//				// Check LOD
//				if !isSignificant && s.SmLOD >= t.TwoBulk.LODP99 {
//					isSignificant = true
//					maxStat = s.SmLOD
//				}
//				// Check BBLogBF
//				if !isSignificant && s.SmBBLogBF >= t.TwoBulk.BBLogBFP99 {
//					isSignificant = true
//					maxStat = s.SmBBLogBF
//				}
//				// Check ED4
//				if !isSignificant && s.SmED >= t.TwoBulk.ED4P99 {
//					isSignificant = true
//					maxStat = s.SmED
//				}
//			}
//
//			if hasOneBulk {
//				// Check AFDev
//				if s.SmAFDev >= t.OneBulk.AFDevP99 || s.SmAFDev <= t.OneBulk.AFDevMp99 {
//					isSignificant = true
//					maxStat = math.Abs(s.SmAFDev)
//				}
//				// Check one-bulk statistics
//				if !isSignificant && s.SmOneBulkG >= t.OneBulk.OneBulkGstatP99 {
//					isSignificant = true
//					maxStat = s.SmOneBulkG
//				}
//				if !isSignificant && s.SmOneBulkLOD >= t.OneBulk.OneBulkLODP99 {
//					isSignificant = true
//					maxStat = s.SmOneBulkLOD
//				}
//				if !isSignificant && s.SmOneBulkBBLogBF >= t.OneBulk.OneBulkBBLogBFP99 {
//					isSignificant = true
//					maxStat = s.SmOneBulkBBLogBF
//				}
//			}
//
//			if isSignificant {
//				if !inRun {
//					inRun = true
//					cur = run{start: s.POS, stop: s.POS, peak: maxStat}
//				} else {
//					cur.stop = s.POS
//					if maxStat > math.Abs(cur.peak) {
//						cur.peak = maxStat
//					}
//				}
//			} else {
//				if inRun {
//					runs = append(runs, cur)
//					inRun = false
//				}
//			}
//		}
//		if inRun {
//			runs = append(runs, cur)
//		}
//
//		if len(runs) == 0 {
//			continue
//		}
//
//		// Keep only the run with the largest peak on this chromosome.
//		best := runs[0]
//		for _, r := range runs[1:] {
//			if math.Abs(r.peak) > math.Abs(best.peak) {
//				best = r
//			}
//		}
//
//		// Find the peak position and CompositeZ value for this run
//		peakZ := 0.0
//		for _, s := range stats {
//			if s.POS >= best.start && s.POS <= best.stop {
//				if math.Abs(s.CompositeZ) > math.Abs(peakZ) {
//					peakZ = s.CompositeZ
//				}
//			}
//		}
//
//		qtls = append(qtls, QTLRecord{
//			Chrom: chrom,
//			Start: best.start,
//			Stop:  best.stop,
//			Peak:  peakZ,
//		})
//	}
//
//	// ── write QTL TSV ────────────────────────────────────────────────────────
//	err := os.MkdirAll(filepath.Join(cfg.OutputDir, "qtls"), 0775)
//	if err != nil {
//		fmt.Println("DetectQTLsWithMCDirect: mkdir: ", err)
//		return nil, err
//	}
//	outPath := filepath.Join(cfg.OutputDir, "qtls", fmt.Sprintf("GoBSAseq.%s.mc_qtls.tsv", bsaType))
//	f, err := os.Create(outPath)
//	if err != nil {
//		return nil, fmt.Errorf("DetectQTLsWithMCDirect: create %s: %w", outPath, err)
//	}
//	defer f.Close()
//
//	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK")
//	for _, q := range qtls {
//		fmt.Fprintf(f, "%s\t%d\t%d\t%.6f\n", q.Chrom, q.Start, q.Stop, q.Peak)
//	}
//
//	color.Green("MC-based QTLs written to %s (%d chromosomes with signal)", outPath, len(qtls))
//	return qtls, nil
//}

func DetectQTLsWithMCDirect(smoothed []SmoothedStats, thresh []Thresholds, bsaType string, cfg *utils.AnalysisConfig) error {
	color.Cyan("--> Running direct single-statistic Monte Carlo evaluations...")

	hasBothBulks := strings.Contains(strings.ToLower(bsaType), "twobulk") || strings.Contains(strings.ToLower(bsaType), "2p2b")
	hasOneBulk := strings.Contains(strings.ToLower(bsaType), "onebulk") || strings.Contains(strings.ToLower(bsaType), "1p1b")

	var activeRuns []SmoothedStats
	var records []DirectQTLRecord

	// Helper function to process and condense a single continuous track block
	flushDirectRun := func(run []SmoothedStats, statsNames []string, tierLabel string) {
		if len(run) == 0 {
			return
		}
		var peakPos int64
		var maxPeak float64
		maxAbs := -1.0

		for _, s := range run {
			if math.Abs(s.CompositeZ) > maxAbs {
				maxAbs = math.Abs(s.CompositeZ)
				peakPos = s.POS
				maxPeak = s.CompositeZ
			}
		}

		records = append(records, DirectQTLRecord{
			Chrom:      run[0].CHROM,
			Start:      run[0].POS,
			Stop:       run[len(run)-1].POS,
			PeakPos:    peakPos,
			ZPeak:      maxPeak,
			Statistics: strings.Join(statsNames, ","),
			Threshold:  tierLabel,
		})
	}

	for i, s := range smoothed {
		t := thresh[i]
		var crossedP99 []string
		var crossedP95 []string

		if hasBothBulks {
			tb := t.TwoBulk
			// P99 Evaluators (Signed boundary fixes included)
			if s.SmDeltaSI >= tb.DeltaSIP99 || s.SmDeltaSI <= tb.DeltaSIMp99 {
				crossedP99 = append(crossedP99, "DeltaSI")
			}
			if s.SmGstat >= tb.GstatP99 {
				crossedP99 = append(crossedP99, "Gstat")
			}
			if s.SmLOD >= tb.LODP99 {
				crossedP99 = append(crossedP99, "LOD")
			}
			if s.SmBBLogBF >= tb.BBLogBFP99 {
				crossedP99 = append(crossedP99, "BBLogBF")
			}

			// P95 Evaluators
			if s.SmDeltaSI >= tb.DeltaSIP95 || s.SmDeltaSI <= tb.DeltaSIMp95 {
				crossedP95 = append(crossedP95, "DeltaSI")
			}
			if s.SmGstat >= tb.GstatP95 {
				crossedP95 = append(crossedP95, "Gstat")
			}
			if s.SmLOD >= tb.LODP95 {
				crossedP95 = append(crossedP95, "LOD")
			}
			if s.SmBBLogBF >= tb.BBLogBFP95 {
				crossedP95 = append(crossedP95, "BBLogBF")
			}
		}

if hasOneBulk {
			ob := t.OneBulk
			// P99 Evaluators
			if s.SmAFDev >= ob.AFDevP99 || s.SmAFDev <= ob.AFDevMp99 {
				crossedP99 = append(crossedP99, "AFDev")
			}
			if s.SmOneBulkG >= ob.OneBulkGstatP99 {
				crossedP99 = append(crossedP99, "Gstat")
			}
			if s.SmOneBulkLOD >= ob.OneBulkLODP99 {
				crossedP99 = append(crossedP99, "LOD")
			}
			if s.SmOneBulkBBLogBF >= ob.OneBulkBBLogBFP99 {
				crossedP99 = append(crossedP99, "BBLogBF")
			}

			// P95 Evaluators
			if s.SmAFDev >= ob.AFDevP95 || s.SmAFDev <= ob.AFDevMp95 {
				crossedP95 = append(crossedP95, "AFDev")
			}
			if s.SmOneBulkG >= ob.OneBulkGstatP95 {
				crossedP95 = append(crossedP95, "Gstat")
			}
			if s.SmOneBulkLOD >= ob.OneBulkLODP95 {
				crossedP95 = append(crossedP95, "LOD")
			}
			if s.SmOneBulkBBLogBF >= ob.OneBulkBBLogBFP95 {
				crossedP95 = append(crossedP95, "BBLogBF")
			}
		}
			if s.SmOneBulkG >= ob.GstatP99 {
				crossedP99 = append(crossedP99, "Gstat")
			}
			if s.SmOneBulkLOD >= ob.LODP99 {
				crossedP99 = append(crossedP99, "LOD")
			}
			if s.SmOneBulkBBLogBF >= ob.BBLogBFP99 {
				crossedP99 = append(crossedP99, "BBLogBF")
			}

			// P95 Evaluators
			if s.SmAFDev >= ob.AFDevP95 || s.SmAFDev <= ob.AFDevMp95 {
				crossedP95 = append(crossedP95, "AFDev")
			}
			if s.SmOneBulkG >= ob.GstatP95 {
				crossedP95 = append(crossedP95, "Gstat")
			}
			if s.SmOneBulkLOD >= ob.LODP95 {
				crossedP95 = append(crossedP95, "LOD")
			}
			if s.SmOneBulkBBLogBF >= ob.BBLogBFP95 {
				crossedP95 = append(crossedP95, "BBLogBF")
			}
		}

		// Hierarchy determination logic: P99 overrides P95 completely
		if len(crossedP99) > 0 {
			activeRuns = append(activeRuns, s)
			// If at the end of a contiguous chunk or chromosome transitions, dump data
			if i == len(smoothed)-1 || smoothed[i+1].CHROM != s.CHROM {
				flushDirectRun(activeRuns, crossedP99, "P99")
				activeRuns = nil
			}
		} else if len(crossedP95) > 0 {
			activeRuns = append(activeRuns, s)
			if i == len(smoothed)-1 || smoothed[i+1].CHROM != s.CHROM {
				flushDirectRun(activeRuns, crossedP95, "P95")
				activeRuns = nil
			}
		} else {
			// No thresholds crossed, terminate active contiguous window block
			if len(activeRuns) > 0 {
				// Reconstruct context of what the active run contained
				flushDirectRun(activeRuns, []string{"MultiStat"}, "P95")
				activeRuns = nil
			}
		}
	}

	// Save data into individual_stats_qtls.tsv
	outDir := filepath.Join(cfg.OutputDir, "qtls")
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.individual_stats_qtls.tsv", bsaType))

	f, err := os.Create(outPath)
	if err != nil {
		return err
	}
	defer f.Close()

	// Updated column naming convention: ZPEAK, STATISTICS, and THRESHOLD
	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK_POS\tZPEAK\tSTATISTICS\tTHRESHOLD")
	for _, r := range records {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%.6f\t%s\t%s\n",
			r.Chrom, r.Start, r.Stop, r.PeakPos, r.ZPeak, r.Statistics, r.Threshold)
	}

	color.Green("    [OK] Direct single-statistic analysis complete. Written to: %s", outPath)
	return nil
}

// DetectQTLs is an EXPLORATORY / SECONDARY caller based on the CompositeZ track.
//
// It scans for regions where |CompositeZ| exceeds the empirical MC threshold
// (CompositeZP99) and keeps the most extreme run per chromosome.
//
// CompositeZ = Σ Zᵢ / √k  (Stouffer, k=3: ΔSI, G-stat, LOD).
// The empirical threshold is derived from a pipeline-null simulation that models
// the actual two-stage bulk+read sampling process and the genome-wide MAD
// normalisation step.  Even so, CompositeZ combines correlated statistics, so
// its null variance is > 1; residual anti-conservatism is possible.
//
// Use DetectQTLsWithMCDirect for formal, publication-grade QTL calls.
// Use this function for a complementary visual sanity-check or when you want a
// single signed composite score for downstream ranking.

//func DetectQTLs(cfg utils.AnalysisConfig, bsaType string, smoothed []SmoothedStats, thresholds []Thresholds) ([]QTLRecord, error) {
//	color.Cyan("===================================== Detecting QTLs =========================================")
//	if len(smoothed) == 0 {
//		fmt.Println("No smoothed stats to detect QTLs from")
//		return nil, nil
//	}
//
//	// Determine threshold to use
//	useEmpiricalThresholds := len(thresholds) == len(smoothed) && len(thresholds) > 0
//
//	// ── group by chromosome ──────────────────────────────────────────────────
//
//	byChr := make(map[string][]SmoothedStats)
//	byChrThresholds := make(map[string][]Thresholds)
//	for i, s := range smoothed {
//		byChr[s.CHROM] = append(byChr[s.CHROM], s)
//		if useEmpiricalThresholds {
//			byChrThresholds[s.CHROM] = append(byChrThresholds[s.CHROM], thresholds[i])
//		}
//	}
//	chroms := make([]string, 0, len(byChr))
//	for c := range byChr {
//		chroms = append(chroms, c)
//	}
//	sort.Strings(chroms)
//
//	// ── scan each chromosome ─────────────────────────────────────────────────
//
//	var qtls []QTLRecord
//
//	for _, chrom := range chroms {
//		stats := byChr[chrom]
//		threshs := byChrThresholds[chrom]
//
//		// Collect every contiguous run where |CompositeZ| ≥ threshold.
//		type run struct {
//			start, stop int64
//			peak        float64 // signed CompositeZ at the most extreme position
//		}
//
//		var runs []run
//		inRun := false
//		var cur run
//
//		for i, s := range stats {
//			z := s.CompositeZ
//
//			// Determine threshold for this variant
//			threshold := zSig
//			if useEmpiricalThresholds {
//				// Use empirical CompositeZ threshold (P99 = 99.5th percentile for two-tailed)
//				threshold = threshs[i].Z.CompositeZP99
//			}
//
//			if math.Abs(z) >= threshold {
//				if !inRun {
//					inRun = true
//					cur = run{start: s.POS, stop: s.POS, peak: z}
//				} else {
//					cur.stop = s.POS
//					if math.Abs(z) > math.Abs(cur.peak) {
//						cur.peak = z
//					}
//				}
//			} else {
//				if inRun {
//					runs = append(runs, cur)
//					inRun = false
//				}
//			}
//		}
//		if inRun {
//			runs = append(runs, cur)
//		}
//
//		if len(runs) == 0 {
//			continue
//		}
//
//		// Keep only the run with the largest |peak| on this chromosome.
//		best := runs[0]
//		for _, r := range runs[1:] {
//			if math.Abs(r.peak) > math.Abs(best.peak) {
//				best = r
//			}
//		}
//
//		qtls = append(qtls, QTLRecord{
//			Chrom: chrom,
//			Start: best.start,
//			Stop:  best.stop,
//			Peak:  best.peak,
//		})
//	}
//
//	// ── write QTL TSV ────────────────────────────────────────────────────────
//	err := os.MkdirAll(filepath.Join(cfg.OutputDir, "qtls"), 0775)
//	if err != nil {
//		fmt.Println("DetectQTLs: mkdir: ", err)
//		return nil, err
//	}
//	outPath := filepath.Join(cfg.OutputDir, "qtls", fmt.Sprintf("GoBSAseq.%s.qtls.tsv", bsaType))
//	f, err := os.Create(outPath)
//	if err != nil {
//		return nil, fmt.Errorf("DetectQTLs: create %s: %w", outPath, err)
//	}
//	defer f.Close()
//
//	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK")
//	for _, q := range qtls {
//		fmt.Fprintf(f, "%s\t%d\t%d\t%.6f\n", q.Chrom, q.Start, q.Stop, q.Peak)
//	}
//
//	color.Green("QTLs written to %s (%d chromosomes with signal)", outPath, len(qtls))
//	return qtls, nil
//}

func DetectQTLs(smoothed []SmoothedStats, thresh *Thresholds, bsaType string, cfg *utils.AnalysisConfig) ([]QTLRecord, error) {
	color.Cyan("--> Detecting QTLs via hierarchical CompositeZ / MaxAbsZ tracks...")

	// Group variants by chromosome
	chrMap := make(map[string][]SmoothedStats)
	var chroms []string
	for _, s := range smoothed {
		if _, ok := chrMap[s.CHROM]; !ok {
			chroms = append(chroms, s.CHROM)
		}
		chrMap[s.CHROM] = append(chrMap[s.CHROM], s)
	}
	sort.Strings(chroms)

	var allQTLs []QTLRecord

	// Process each chromosome independently using the fallback hierarchy
	for _, chrom := range chroms {
		variants := chrMap[chrom]
		if len(variants) == 0 {
			continue
		}

		var chrQTLs []QTLRecord
		calledStat := ""
		calledThresh := ""

		// Hierarchy Tier 1: Stouffer (CompositeZ) P99
		chrQTLs = findZSegments(variants, thresh.Z.CompositeZP99, true)
		if len(chrQTLs) > 0 {
			calledStat = "Stouffer"
			calledThresh = "P99"
		} else {
			// Hierarchy Tier 2: Stouffer (CompositeZ) P95
			chrQTLs = findZSegments(variants, thresh.Z.CompositeZP95, true)
			if len(chrQTLs) > 0 {
				calledStat = "Stouffer"
				calledThresh = "P95"
			} else {
				// Hierarchy Tier 3: MaxAbsZ P99
				chrQTLs = findZSegments(variants, thresh.Z.ZP99, false)
				if len(chrQTLs) > 0 {
					calledStat = "MaxAbsZ"
					calledThresh = "P99"
				} else {
					// Hierarchy Tier 4: MaxAbsZ P95
					chrQTLs = findZSegments(variants, thresh.Z.ZP95, false)
					if len(chrQTLs) > 0 {
						calledStat = "MaxAbsZ"
						calledThresh = "P95"
					}
				}
			}
		}

		// Fill in metadata for identified peaks on this chromosome
		for i := range chrQTLs {
			chrQTLs[i].Stat = calledStat
			chrQTLs[i].Threshold = calledThresh
			allQTLs = append(allQTLs, chrQTLs[i])
		}
	}

	// Write output to qtls.tsv
	outDir := filepath.Join(cfg.OutputDir, "qtls")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.qtls.tsv", bsaType))
	f, err := os.Create(outPath)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	// Added STAT and THRESHOLD columns
	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK_VAL\tSTAT\tTHRESHOLD")
	for _, q := range allQTLs {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%.6f\t%s\t%s\n",
			q.Chrom, q.Start, q.Stop, q.PeakPos, q.PeakVal, q.Stat, q.Threshold)
	}

	color.Green("    [OK] Hierarchical QTL analysis completed. Written to: %s", outPath)
	return allQTLs, nil
}

func findZSegments(variants []SmoothedStats, cutoff float64, useComposite bool) []QTLRecord {
	var records []QTLRecord
	var currentRun []SmoothedStats

	appendRun := func() {
		if len(currentRun) == 0 {
			return
		}
		start := currentRun[0].POS
		stop := currentRun[len(currentRun)-1].POS

		// Find peak position based on max absolute value score
		var peakPos int64
		var peakVal float64
		maxAbs := -1.0

		for _, v := range currentRun {
			val := v.MaxAbsZ
			if useComposite {
				val = v.CompositeZ
			}
			if math.Abs(val) > maxAbs {
				maxAbs = math.Abs(val)
				peakPos = v.POS
				peakVal = val
			}
		}

		records = append(records, QTLRecord{
			Chrom:   currentRun[0].CHROM,
			Start:   start,
			Stop:    stop,
			PeakPos: peakPos,
			PeakVal: peakVal,
		})
		currentRun = nil
	}

	for _, v := range variants {
		score := v.MaxAbsZ
		if useComposite {
			score = v.CompositeZ
		}

		if math.Abs(score) >= cutoff {
			currentRun = append(currentRun, v)
		} else {
			appendRun()
		}
	}
	appendRun()

	return records
}

// ── Merging ───────────────────────────────────────────────────────────────────

// MergeQTLsAndBRM combines QTLRecord and BRMBlock results into a single
// []MergedQTL, unions overlapping / same-chromosome intervals, and writes a
// merged TSV.
//
// Merging rules (per chromosome):
//   - Intervals from both sources on the same chromosome are unioned
//     (Start = min, Stop = max) and labelled "ZScore+BRM".
//   - Intervals present in only one source keep their source label.
//   - Multiple BRM blocks on the same chromosome are first reduced to their
//     bounding interval (widest span, highest |peak|) before merging with the
//     Z-score QTL.  This mirrors the legacy behaviour where only the most
//     extreme block was forwarded to gene space analysis.
func MergeQTLsAndBRM(cfg utils.AnalysisConfig, bsaType string, qtls []QTLRecord, blocks []BRMBlock) ([]MergedQTL, error) {

	nan := math.NaN()

	// Index QTL records by chromosome.
	qtlByChrom := make(map[string]QTLRecord, len(qtls))
	for _, q := range qtls {
		qtlByChrom[q.Chrom] = q
	}

	// Reduce BRM blocks to one representative per chromosome:
	// widest bounding interval, highest |peak| block retained.
	type brmSummary struct {
		start, stop int64
		peak        float64
		threshold   float64
	}
	brmByChrom := make(map[string]brmSummary)
	for _, b := range blocks {
		prev, exists := brmByChrom[b.Chrom]
		if !exists {
			brmByChrom[b.Chrom] = brmSummary{b.Start, b.Stop, b.Peak, b.Threshold}
			continue
		}
		s := prev
		if b.Start < s.start {
			s.start = b.Start
		}
		if b.Stop > s.stop {
			s.stop = b.Stop
		}
		if math.Abs(b.Peak) > math.Abs(s.peak) {
			s.peak = b.Peak
			s.threshold = b.Threshold
		}
		brmByChrom[b.Chrom] = s
	}

	// Build the union set of chromosomes.
	chromSet := make(map[string]struct{})
	for c := range qtlByChrom {
		chromSet[c] = struct{}{}
	}
	for c := range brmByChrom {
		chromSet[c] = struct{}{}
	}
	chroms := make([]string, 0, len(chromSet))
	for c := range chromSet {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	var merged []MergedQTL
	for _, chrom := range chroms {
		q, hasQ := qtlByChrom[chrom]
		b, hasB := brmByChrom[chrom]

		m := MergedQTL{
			Chrom:        chrom,
			ZPeak:        nan,
			BRMPeak:      nan,
			BRMThreshold: nan,
		}

		switch {
		case hasQ && hasB:
			// Union the two intervals.
			m.Start = minI64(q.Start, b.start)
			m.Stop = maxI64(q.Stop, b.stop)
			m.ZPeak = q.PeakVal
			m.BRMPeak = b.peak
			m.BRMThreshold = b.threshold
			m.Source = "ZScore+BRM"

		case hasQ:
			m.Start = q.Start
			m.Stop = q.Stop
			m.ZPeak = q.Peak
			m.Source = "ZScore"

		case hasB:
			m.Start = b.start
			m.Stop = b.stop
			m.BRMPeak = b.peak
			m.BRMThreshold = b.threshold
			m.Source = "BRM"
		}

		merged = append(merged, m)
	}

	// ── write merged TSV ─────────────────────────────────────────────────────

	outDir := filepath.Join(cfg.OutputDir, "qtls")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, fmt.Errorf("MergeQTLsAndBRM: mkdir: %w", err)
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.merged_qtls.tsv", bsaType))

	f, err := os.Create(outPath)
	if err != nil {
		return nil, fmt.Errorf("MergeQTLsAndBRM: create %s: %w", outPath, err)
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tSOURCE\tZ_PEAK\tBRM_PEAK\tBRM_THRESHOLD")
	for _, m := range merged {
		fmt.Fprintf(f, "%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
			m.Chrom, m.Start, m.Stop, m.Source,
			nanOrFloat(m.ZPeak),
			nanOrFloat(m.BRMPeak),
			nanOrFloat(m.BRMThreshold),
		)
	}

	color.Green("Merged QTL+BRM intervals written to %s (%d intervals)", outPath, len(merged))
	return merged, nil
}

// ── helpers ───────────────────────────────────────────────────────────────────

func minI64(a, b int64) int64 {
	if a < b {
		return a
	}
	return b
}

func maxI64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}

// nanOrFloat formats v as "NA" when it is NaN, otherwise as %.6f.
func nanOrFloat(v float64) string {
	if math.IsNaN(v) {
		return "NA"
	}
	return fmt.Sprintf("%.6f", v)
}
