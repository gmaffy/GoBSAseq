package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// Thresholds are the genome-wide significance thresholds of one analysis,
// from the genome-wide null simulation (hmmnull.go). Each is a quantile of
// the most extreme value a statistic reaches anywhere in a simulated no-QTL
// genome, as scanned (smoothed tracks after smoothing, the HMM LOD as is), so
// P99 means a 1% chance that a genome without QTLs crosses it anywhere.
//
// One-sided statistics (G′, ED⁴, G1, |ΔSI|, HMM LOD, CompositeZ) use the
// level whole: P99/P95 are the 99%/95% quantiles of the genome maxima.
// Two-tailed ones (ΔSI, AFDev, SI) split it between the tails: P99/P95 are
// the 99.5%/97.5% quantiles of the maxima, and Mp99/Mp95 the matching lower
// quantiles of the minima.
type Thresholds struct {
	TwoBulk  TwoBulkThresholds
	OneBulk  OneBulkThresholds
	Z        ZThresholds
	HMM      HMMThresholds
	Rho      float64 // read overdispersion simulated in the null (estimateRho)
	NullReps int     // simulated no-QTL genomes
	// BRMU is the BRM threshold on the standardised statistic |z| (brmModel):
	// the (1 − BRMAlpha) quantile of the genome-wide maximum |z| of the null
	// genomes, so BRMAlpha (--brm-alpha) is a genome-wide error rate. 0 when
	// not computed (BRMAlpha outside (0, 1)).
	BRMU, BRMAlpha float64
}

// HMMThresholds are the HMM LOD thresholds; zero when it is not computed.
type HMMThresholds struct {
	LODP99, LODP95 float64
}

// TwoBulkThresholds are the smoothed two-bulk tracks' thresholds. The SI
// fields are also used for the bulk's SI in one-bulk modes. With unknown
// polarity (2b, --ref-parent none) the ΔSI track is |ΔSI|: one-sided, and the
// lower-tail fields are NaN.
type TwoBulkThresholds struct {
	HighSIP99, HighSIP95     float64
	HighSIMp99, HighSIMp95   float64
	LowSIP99, LowSIP95       float64
	LowSIMp99, LowSIMp95     float64
	DeltaSIP99, DeltaSIP95   float64
	DeltaSIMp99, DeltaSIMp95 float64

	GstatP99, GstatP95 float64
	ED4P99, ED4P95     float64
}

// OneBulkThresholds are the smoothed one-bulk tracks' thresholds.
type OneBulkThresholds struct {
	AFDevP99, AFDevP95               float64
	AFDevMp99, AFDevMp95             float64 // lower tail
	OneBulkGstatP99, OneBulkGstatP95 float64
}

// ZThresholds are the CompositeZ thresholds; NumStats is the number of
// Z-scores CompositeZ combines.
type ZThresholds struct {
	CompositeZP99, CompositeZP95 float64
	NumStats                     int
}

// countZStats returns k, the number of Z-scores consolidate combines into
// CompositeZ (smoothLayout.compositeCount): 2 in every mode, |ZΔSI| + ZG′ or
// |ZAFDev| + ZG1. The HMM LOD is a separate track.
func countZStats(bsaType string) int {
	hb, lb, _, _ := BulkFlags(bsaType)
	return smoothLayout{highBulk: hb, lowBulk: lb, signedDSI: true}.compositeCount()
}

// CalculateThresholds runs the genome-wide null simulation over raw's markers
// and depths and returns the analysis's thresholds, also written to
// stats/GoBSAseq.<bsaType>.thresholds.tsv.
func CalculateThresholds(cfg utils.AnalysisConfig, bsaType string, raw []BSAstats) (Thresholds, error) {
	if _, _, _, err := cfg.BulkSizes(bsaType); err != nil {
		return Thresholds{}, err
	}
	if !(cfg.BrmAlpha > 0 && cfg.BrmAlpha < 1) {
		return Thresholds{}, fmt.Errorf("--brm-alpha must be in (0, 1), got %v", cfg.BrmAlpha)
	}
	color.Cyan("\n============================ Calculating genome-wide thresholds (%s) ============================\n", bsaType)
	th, err := genomeNullThresholds(cfg, bsaType, raw)
	if err != nil {
		return Thresholds{}, fmt.Errorf("genome-wide null simulation: %w", err)
	}
	outPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.thresholds.tsv", bsaType))
	if err := writeThresholdsTSV(outPath, th, bsaType, cfg.SignedDeltaSI(bsaType), HMMAvailable(cfg, bsaType)); err != nil {
		return th, fmt.Errorf("failed to write thresholds TSV: %w", err)
	}
	return th, nil
}

// writeThresholdsTSV writes one row per scanned statistic and tail.
func writeThresholdsTSV(path string, t Thresholds, bsaType string, signedDSI, hmm bool) (err error) {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	w := bufio.NewWriter(f)
	fmt.Fprintf(w, "# Genome-wide thresholds: quantiles of the extreme value over %d simulated no-QTL genomes (generalized extreme-value fit); read overdispersion rho = %.4f\n", t.NullReps, t.Rho)
	fmt.Fprintln(w, "Statistic\tTail\tP99\tP95")
	row := func(stat, tail string, p99, p95 float64) {
		fmt.Fprintf(w, "%s\t%s\t%.6g\t%.6g\n", stat, tail, p99, p95)
	}
	hb, lb, both, one := BulkFlags(bsaType)
	tb, ob := t.TwoBulk, t.OneBulk
	if hb {
		row("Sm_HighSI", "upper", tb.HighSIP99, tb.HighSIP95)
		row("Sm_HighSI", "lower", tb.HighSIMp99, tb.HighSIMp95)
	}
	if lb {
		row("Sm_LowSI", "upper", tb.LowSIP99, tb.LowSIP95)
		row("Sm_LowSI", "lower", tb.LowSIMp99, tb.LowSIMp95)
	}
	if both {
		if signedDSI {
			row("Sm_DeltaSI", "upper", tb.DeltaSIP99, tb.DeltaSIP95)
			row("Sm_DeltaSI", "lower", tb.DeltaSIMp99, tb.DeltaSIMp95)
		} else {
			row("Sm_AbsDeltaSI", "upper", tb.DeltaSIP99, tb.DeltaSIP95)
		}
		row("Sm_Gstat", "upper", tb.GstatP99, tb.GstatP95)
		row("Sm_ED4", "upper", tb.ED4P99, tb.ED4P95)
	}
	if one {
		row("Sm_AFDev", "upper", ob.AFDevP99, ob.AFDevP95)
		row("Sm_AFDev", "lower", ob.AFDevMp99, ob.AFDevMp95)
		row("Sm_OneBulkGstat", "upper", ob.OneBulkGstatP99, ob.OneBulkGstatP95)
	}
	if hmm {
		row("HMMLOD", "upper", t.HMM.LODP99, t.HMM.LODP95)
	}
	row("CompositeZ", "upper", t.Z.CompositeZP99, t.Z.CompositeZP95)
	if t.BRMU > 0 {
		fmt.Fprintf(w, "# BRM: blocks where |z| = |signal|/σ ≥ %.6g (genome-wide alpha = %g)\n", t.BRMU, t.BRMAlpha)
	}
	if err := w.Flush(); err != nil {
		return err
	}
	color.Green("✔ Thresholds written to %s\n", path)
	return nil
}

// r6 rounds a float64 to 6 decimal places (NaN and ±Inf unchanged).
func r6(v float64) float64 {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return v
	}
	return math.Round(v*1e6) / 1e6
}
