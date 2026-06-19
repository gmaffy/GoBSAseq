package stats

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// ── QTL detection ─────────────────────────────────────────────────────────────

// QTLRecord holds the single dominant QTL interval detected per chromosome.
type QTLRecord struct {
	Chrom string
	Start int64
	Stop  int64
	Peak  float64 // CompositeZ value at the peak position
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

// DetectQTLsWithMCDirect uses the per-variant Monte Carlo thresholds directly for QTL detection.
// This provides an alternative to the CompositeZ ≥ 3.0 rule, using the fully sound Monte Carlo
// thresholds that account for depth and bulk size.
// Returns QTL regions where statistics exceed their MC-derived thresholds.
func DetectQTLsWithMCDirect(cfg utils.AnalysisConfig, bsaType string, smoothed []SmoothedStats, thresholds []Thresholds) ([]QTLRecord, error) {
	color.Cyan("===================================== Detecting QTLs with MC Thresholds =========================================")
	if len(smoothed) == 0 {
		fmt.Println("No smoothed stats to detect QTLs from")
		return nil, nil
	}
	if len(thresholds) != len(smoothed) {
		return nil, fmt.Errorf("smoothed and thresholds slices must have the same length")
	}

	// ── group by chromosome ──────────────────────────────────────────────────

	byChr := make(map[string][]SmoothedStats)
	byChrThresholds := make(map[string][]Thresholds)
	for i, s := range smoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
		byChrThresholds[s.CHROM] = append(byChrThresholds[s.CHROM], thresholds[i])
	}
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// ── scan each chromosome ─────────────────────────────────────────────────

	var qtls []QTLRecord

	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	for _, chrom := range chroms {
		stats := byChr[chrom]
		threshs := byChrThresholds[chrom]

		// Collect every contiguous run where key statistics exceed their MC thresholds.
		type run struct {
			start, stop int64
			peak        float64 // the most extreme statistic value at the peak
		}

		var runs []run
		inRun := false
		var cur run

		for i, s := range stats {
			t := threshs[i]
			
			// Check which statistics are significant based on MC thresholds
			// We use a combined approach: any of the main statistics exceeding threshold
			isSignificant := false
			var maxStat float64
			
			if hasBothBulks {
				// Check DeltaSI against MC threshold
				if math.Abs(s.SmDeltaSI) >= t.TwoBulk.DeltaSIP99 || math.Abs(s.SmDeltaSI) <= t.TwoBulk.DeltaSIMp99 {
					isSignificant = true
					maxStat = math.Abs(s.SmDeltaSI)
				}
				// Check G-statistic
				if !isSignificant && s.SmGstat >= t.TwoBulk.GstatP99 {
					isSignificant = true
					maxStat = s.SmGstat
				}
				// Check LOD
				if !isSignificant && s.SmLOD >= t.TwoBulk.LODP99 {
					isSignificant = true
					maxStat = s.SmLOD
				}
				// Check BBLogBF
				if !isSignificant && s.SmBBLogBF >= t.TwoBulk.BBLogBFP99 {
					isSignificant = true
					maxStat = s.SmBBLogBF
				}
				// Check ED4
				if !isSignificant && s.SmED >= t.TwoBulk.ED4P99 {
					isSignificant = true
					maxStat = s.SmED
				}
			}
			
			if hasOneBulk {
				// Check AFDev
				if math.Abs(s.SmAFDev) >= t.OneBulk.AFDevP99 || math.Abs(s.SmAFDev) <= t.OneBulk.AFDevMp99 {
					isSignificant = true
					maxStat = math.Abs(s.SmAFDev)
				}
				// Check one-bulk statistics
				if !isSignificant && s.SmOneBulkG >= t.OneBulk.OneBulkGstatP99 {
					isSignificant = true
					maxStat = s.SmOneBulkG
				}
				if !isSignificant && s.SmOneBulkLOD >= t.OneBulk.OneBulkLODP99 {
					isSignificant = true
					maxStat = s.SmOneBulkLOD
				}
				if !isSignificant && s.SmOneBulkBBLogBF >= t.OneBulk.OneBulkBBLogBFP99 {
					isSignificant = true
					maxStat = s.SmOneBulkBBLogBF
				}
			}

			if isSignificant {
				if !inRun {
					inRun = true
					cur = run{start: s.POS, stop: s.POS, peak: maxStat}
				} else {
					cur.stop = s.POS
					if maxStat > math.Abs(cur.peak) {
						cur.peak = maxStat
					}
				}
			} else {
				if inRun {
					runs = append(runs, cur)
					inRun = false
				}
			}
		}
		if inRun {
			runs = append(runs, cur)
		}

		if len(runs) == 0 {
			continue
		}

		// Keep only the run with the largest peak on this chromosome.
		best := runs[0]
		for _, r := range runs[1:] {
			if math.Abs(r.peak) > math.Abs(best.peak) {
				best = r
			}
		}

		// Find the peak position and CompositeZ value for this run
		peakZ := 0.0
		for _, s := range stats {
			if s.POS >= best.start && s.POS <= best.stop {
				if math.Abs(s.CompositeZ) > math.Abs(peakZ) {
					peakZ = s.CompositeZ
				}
			}
		}

		qtls = append(qtls, QTLRecord{
			Chrom: chrom,
			Start: best.start,
			Stop:  best.stop,
			Peak:  peakZ,
		})
	}

	// ── write QTL TSV ────────────────────────────────────────────────────────
	err := os.MkdirAll(filepath.Join(cfg.OutputDir, "qtls"), 0775)
	if err != nil {
		fmt.Println("DetectQTLsWithMCDirect: mkdir: ", err)
		return nil, err
	}
	outPath := filepath.Join(cfg.OutputDir, "qtls", fmt.Sprintf("GoBSAseq.%s.mc_qtls.tsv", bsaType))
	f, err := os.Create(outPath)
	if err != nil {
		return nil, fmt.Errorf("DetectQTLsWithMCDirect: create %s: %w", outPath, err)
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK")
	for _, q := range qtls {
		fmt.Fprintf(f, "%s\t%d\t%d\t%.6f\n", q.Chrom, q.Start, q.Stop, q.Peak)
	}

	color.Green("MC-based QTLs written to %s (%d chromosomes with signal)", outPath, len(qtls))
	return qtls, nil
}

// DetectQTLs scans the smoothed CompositeZ track for regions that exceed a threshold
// and keeps the most extreme run per chromosome. Results are written to a TSV.
// If thresholds are provided, uses empirical CompositeZ thresholds; otherwise uses zSig.
func DetectQTLs(cfg utils.AnalysisConfig, bsaType string, smoothed []SmoothedStats, thresholds []Thresholds) ([]QTLRecord, error) {
	color.Cyan("===================================== Detecting QTLs =========================================")
	if len(smoothed) == 0 {
		fmt.Println("No smoothed stats to detect QTLs from")
		return nil, nil
	}

	// Determine threshold to use
	useEmpiricalThresholds := len(thresholds) == len(smoothed) && len(thresholds) > 0
	
	// ── group by chromosome ──────────────────────────────────────────────────

	byChr := make(map[string][]SmoothedStats)
	byChrThresholds := make(map[string][]Thresholds)
	for i, s := range smoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
		if useEmpiricalThresholds {
			byChrThresholds[s.CHROM] = append(byChrThresholds[s.CHROM], thresholds[i])
		}
	}
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// ── scan each chromosome ─────────────────────────────────────────────────

	var qtls []QTLRecord

	for _, chrom := range chroms {
		stats := byChr[chrom]
		threshs := byChrThresholds[chrom]

		// Collect every contiguous run where |CompositeZ| ≥ threshold.
		type run struct {
			start, stop int64
			peak        float64 // signed CompositeZ at the most extreme position
		}

		var runs []run
		inRun := false
		var cur run

		for i, s := range stats {
			z := s.CompositeZ
			
			// Determine threshold for this variant
			threshold := zSig
			if useEmpiricalThresholds {
				// Use empirical CompositeZ threshold (P99 = 99.5th percentile for two-tailed)
				threshold = threshs[i].Z.CompositeZP99
			}
			
			if math.Abs(z) >= threshold {
				if !inRun {
					inRun = true
					cur = run{start: s.POS, stop: s.POS, peak: z}
				} else {
					cur.stop = s.POS
					if math.Abs(z) > math.Abs(cur.peak) {
						cur.peak = z
					}
				}
			} else {
				if inRun {
					runs = append(runs, cur)
					inRun = false
				}
			}
		}
		if inRun {
			runs = append(runs, cur)
		}

		if len(runs) == 0 {
			continue
		}

		// Keep only the run with the largest |peak| on this chromosome.
		best := runs[0]
		for _, r := range runs[1:] {
			if math.Abs(r.peak) > math.Abs(best.peak) {
				best = r
			}
		}

		qtls = append(qtls, QTLRecord{
			Chrom: chrom,
			Start: best.start,
			Stop:  best.stop,
			Peak:  best.peak,
		})
	}

	// ── write QTL TSV ────────────────────────────────────────────────────────
	err := os.MkdirAll(filepath.Join(cfg.OutputDir, "qtls"), 0775)
	if err != nil {
		fmt.Println("DetectQTLs: mkdir: ", err)
		return nil, err
	}
	outPath := filepath.Join(cfg.OutputDir, "qtls", fmt.Sprintf("GoBSAseq.%s.qtls.tsv", bsaType))
	f, err := os.Create(outPath)
	if err != nil {
		return nil, fmt.Errorf("DetectQTLs: create %s: %w", outPath, err)
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tSTART\tSTOP\tPEAK")
	for _, q := range qtls {
		fmt.Fprintf(f, "%s\t%d\t%d\t%.6f\n", q.Chrom, q.Start, q.Stop, q.Peak)
	}

	color.Green("QTLs written to %s (%d chromosomes with signal)", outPath, len(qtls))
	return qtls, nil
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
			m.ZPeak = q.Peak
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
