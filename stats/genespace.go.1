package stats

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/gmaffy/genome-whisperer/annotation"
	"github.com/gmaffy/genome-whisperer/genespace"
)

// GeneSpaceInterval is a genomic region to analyse — satisfied by both
// QTLRecord and BRMBlock so callers can pass either type.
type GeneSpaceInterval interface {
	geneSpaceRegion() (chrom string, start, stop int)
}

func (q QTLRecord) geneSpaceRegion() (string, int, int) {
	return q.Chrom, int(q.Start), int(q.Stop)
}

func (b BRMBlock) geneSpaceRegion() (string, int, int) {
	return b.Chrom, int(b.Start), int(b.Stop)
}

// RunGeneSpace performs SnpEff annotation of the filtered VCF produced for
// bsaType, then runs genespace.GeneSpace over every supplied interval.
//
// intervals accepts []QTLRecord, []BRMBlock, or any mix via the
// GeneSpaceInterval interface.  Pass the merged slice from MergeQTLsAndBRM
// (detect.go) for the most complete coverage.
//
// The function is a no-op (with a yellow warning) when any of the four
// required config fields (SnpEffDB, GeneDesc, Prg, Gff) is absent.
func RunGeneSpace(cfg utils.AnalysisConfig, bsaType string, filteredVcfPath string, intervals []GeneSpaceInterval) error {
	color.Cyan("\n============================= Gene space analysis (%s) ==============================\n", bsaType)

	// ── prerequisite check ───────────────────────────────────────────────────
	if cfg.SnpEffDB == "" || cfg.GeneDesc == "" || cfg.Prg == "" || cfg.Gff == "" {
		color.Yellow("Skipping gene space analysis: one or more required parameters are missing " +
			"(SnpEffDB, GeneDesc, Prg, Gff).")
		return nil
	}
	if len(intervals) == 0 {
		color.Yellow("Skipping gene space analysis: no intervals to analyse.")
		return nil
	}

	outDir := filepath.Join(cfg.OutputDir, "gene_space")
	err := os.MkdirAll(outDir, 0775)
	if err != nil {
		return err
	}
	//outDir = filepath.Join(outDir, bsaType,)
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return fmt.Errorf("RunGeneSpace: create output dir: %w", err)
	}

	// ── Step 1: SnpEff annotation ─────────────────────────────────────────────
	// Only annotate if the VCF does not already carry EFF annotations.  The
	// legacy twoBulk.go inspected the live vcfgo reader header; here we check
	// via cfg.Rdr when it is set, and skip the guard when it is nil so that
	// callers that do not preserve the reader can still use this function.
	if cfg.Rdr != nil {
		if _, hasEFF := cfg.Rdr.Header.Infos["EFF"]; hasEFF {
			color.Yellow("EFF column already present in VCF; skipping SnpEff annotation step.")
			color.Yellow("Gene space analysis requires a freshly annotated TSV — skipping.")
			return nil
		}
	}

	color.Green("Step 1: Annotating variants with SnpEff ...")
	err, annotatedTsvFiles := annotation.CreateSuperVcf(
		[]string{filteredVcfPath},
		cfg.SnpEffDB,
		true,
		cfg.GeneDesc,
		cfg.Prg,
	)
	if err != nil {
		return fmt.Errorf("RunGeneSpace: SnpEff annotation failed: %w", err)
	}
	if len(annotatedTsvFiles) == 0 {
		color.Yellow("Skipping gene space analysis: SnpEff produced no annotated TSV files.")
		return nil
	}
	color.Green("SnpEff annotation complete.  Annotated TSV: %v", annotatedTsvFiles)

	annotatedTsv := annotatedTsvFiles[0]

	// ── Step 2: GeneSpace per interval ───────────────────────────────────────
	color.Cyan("Step 2: Running gene space analysis for %d interval(s) ...", len(intervals))

	for _, iv := range intervals {
		chrom, start, stop := iv.geneSpaceRegion()
		color.Blue("  Interval: %s:%d-%d", chrom, start, stop)

		_, err := genespace.GeneSpace(
			cfg.Gff,
			annotatedTsv,
			chrom,
			start,
			stop,
			[]string{}, // included effect types  — empty = all
			[]string{}, // excluded effect types  — empty = none
			cfg.GeneDesc,
			cfg.Prg,
			outDir,
		)
		if err != nil {
			// Log and continue rather than aborting the whole run.
			color.Red("  gene space analysis failed for %s:%d-%d: %v", chrom, start, stop, err)
		}
	}

	color.Green("Gene space analysis complete.")
	return nil
}

// GeneSpaceFromQTLs is a convenience wrapper that converts []QTLRecord to the
// GeneSpaceInterval slice expected by RunGeneSpace.
func GeneSpaceFromQTLs(
	cfg utils.AnalysisConfig,
	bsaType string,
	filteredVcfPath string,
	qtls []QTLRecord,
) error {
	ivs := make([]GeneSpaceInterval, len(qtls))
	for i, q := range qtls {
		ivs[i] = q
	}
	return RunGeneSpace(cfg, bsaType, filteredVcfPath, ivs)
}

// GeneSpaceFromBRM is a convenience wrapper that converts []BRMBlock to the
// GeneSpaceInterval slice expected by RunGeneSpace.
func GeneSpaceFromBRM(
	cfg utils.AnalysisConfig,
	bsaType string,
	filteredVcfPath string,
	blocks []BRMBlock,
) error {
	ivs := make([]GeneSpaceInterval, len(blocks))
	for i, b := range blocks {
		ivs[i] = b
	}
	return RunGeneSpace(cfg, bsaType, filteredVcfPath, ivs)
}

// GeneSpaceFromMerged is a convenience wrapper that accepts the MergedQTL
// slice produced by MergeQTLsAndBRM (detect.go) and runs gene space analysis
// over all merged intervals in one pass.
func GeneSpaceFromMerged(cfg utils.AnalysisConfig, bsaType string, filteredVcfPath string, merged []MergedQTL) error {
	ivs := make([]GeneSpaceInterval, len(merged))
	for i, m := range merged {
		ivs[i] = m
	}
	return RunGeneSpace(cfg, bsaType, filteredVcfPath, ivs)
}
