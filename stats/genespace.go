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

type GeneSpaceInterval interface {
	geneSpaceRegion() (chrom string, start, stop int)
}

func (b BRMBlock) geneSpaceRegion() (string, int, int) {
	return b.Chrom, int(b.Start), int(b.Stop)
}

// gtCol returns the VariantsToTable genotype column name for a VCF sample.
func gtCol(sample string) string {
	if sample == "" || sample == "None" {
		return ""
	}
	return sample + ".GT"
}

func nonEmptyGTCols(names ...string) []string {
	var cols []string
	for _, name := range names {
		if col := gtCol(name); col != "" {
			cols = append(cols, col)
		}
	}
	return cols
}

// geneSpaceGTLines maps the current BSA mode to resistant and susceptible
// genotype columns in the annotated SuperVCF TSV. High-parent/high-bulk lines
// are treated as resistant; low-parent/low-bulk as susceptible.
func geneSpaceGTLines(cfg utils.AnalysisConfig, bsaType string) (res, sus []string) {
	switch bsaType {
	case "2p2b":
		return nonEmptyGTCols(cfg.HighParentName, cfg.HighBulkName),
			nonEmptyGTCols(cfg.LowParentName, cfg.LowBulkName)
	case "2b":
		return nonEmptyGTCols(cfg.HighBulkName), nonEmptyGTCols(cfg.LowBulkName)
	case "2phb":
		return nonEmptyGTCols(cfg.HighParentName, cfg.HighBulkName),
			nonEmptyGTCols(cfg.LowParentName)
	case "2plb":
		return nonEmptyGTCols(cfg.HighParentName),
			nonEmptyGTCols(cfg.LowParentName, cfg.LowBulkName)
	case "hp2b":
		return nonEmptyGTCols(cfg.HighParentName, cfg.HighBulkName),
			nonEmptyGTCols(cfg.LowBulkName)
	case "lp2b":
		return nonEmptyGTCols(cfg.HighBulkName),
			nonEmptyGTCols(cfg.LowParentName, cfg.LowBulkName)
	case "hphb":
		return nonEmptyGTCols(cfg.HighParentName, cfg.HighBulkName), nil
	case "hplb":
		return nonEmptyGTCols(cfg.HighParentName), nonEmptyGTCols(cfg.LowBulkName)
	case "lphb":
		return nonEmptyGTCols(cfg.HighBulkName), nonEmptyGTCols(cfg.LowParentName)
	case "lplb":
		return nonEmptyGTCols(cfg.LowParentName), nonEmptyGTCols(cfg.LowBulkName)
	default:
		return nil, nil
	}
}

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
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return fmt.Errorf("RunGeneSpace: create output dir: %w", err)
	}

	if cfg.Rdr != nil {
		if _, hasEFF := cfg.Rdr.Header.Infos["EFF"]; hasEFF {
			color.Yellow("EFF column already present in VCF; skipping SnpEff annotation step.")
			color.Yellow("Gene space analysis requires a freshly annotated TSV — skipping.")
			return nil
		}
	}

	color.Green("Step 1: Annotating variants with SnpEff ...")
	annotErr, annotatedTsvFiles := annotation.CreateSuperVcf([]string{filteredVcfPath}, cfg.SnpEffDB, true, cfg.GeneDesc, cfg.Prg)
	
	if annotErr != nil {
		return fmt.Errorf("RunGeneSpace: SnpEff annotation failed: %w", annotErr)
	}
	if len(annotatedTsvFiles) == 0 {
		color.Yellow("Skipping gene space analysis: SnpEff produced no annotated TSV files.")
		return nil
	}
	color.Green("SnpEff annotation complete.  Annotated TSV: %v", annotatedTsvFiles)

	annotatedTsv := annotatedTsvFiles[0]

	resLines, susLines := geneSpaceGTLines(cfg, bsaType)
	if len(resLines) == 0 || len(susLines) == 0 {
		color.Yellow("Gene space SNP counts require both resistant and susceptible sample columns; got res=%v sus=%v",
			resLines, susLines)
	} else {
		color.Blue("  Resistant GT columns: %v", resLines)
		color.Blue("  Susceptible GT columns: %v", susLines)
	}

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
			resLines,
			susLines,
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

func GeneSpaceFromBRM(cfg utils.AnalysisConfig, bsaType string, filteredVcfPath string, blocks []BRMBlock) error {
	ivs := make([]GeneSpaceInterval, len(blocks))
	for i, b := range blocks {
		ivs[i] = b
	}
	return RunGeneSpace(cfg, bsaType, filteredVcfPath, ivs)
}

func GeneSpaceFromMerged(cfg utils.AnalysisConfig, bsaType string, filteredVcfPath string, merged []MergedQTL) error {
	ivs := make([]GeneSpaceInterval, len(merged))
	for i, m := range merged {
		ivs[i] = m
	}
	return RunGeneSpace(cfg, bsaType, filteredVcfPath, ivs)
}
