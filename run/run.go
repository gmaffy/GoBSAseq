package run

import (
	"compress/gzip"
	"fmt"
	"io"
	"maps"
	"os"
	"path/filepath"
	"slices"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/plots"
	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
)

func openVCF(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	cleanup := func() { f.Close() }

	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, err
		}
		cleanup = func() { gz.Close(); f.Close() }
		return gz, cleanup, nil
	}
	return f, cleanup, nil
}

func promptForSampleChoice(sampleMap map[int]string, label string) (int, error) {
	keys := slices.Sorted(maps.Keys(sampleMap))
	for _, i := range keys {
		fmt.Printf("%v : %v\n", i, sampleMap[i])
	}

	color.Blue("\nEnter %s number:", label)
	var choice int
	if _, err := fmt.Scan(&choice); err != nil {
		return 0, fmt.Errorf("%s number should be numerical and part of the list above: %w", label, err)
	}
	if !slices.Contains(keys, choice) {
		return 0, fmt.Errorf("invalid input")
	}
	return choice, nil
}

func sampleIndexByName(sampleMap map[int]string, name string) (int, bool) {
	for k, v := range sampleMap {
		if v == name {
			return k, true
		}
	}
	return 0, false
}

func missingGeneSpaceParams(cfg *utils.AnalysisConfig) []string {
	var missing []string
	if cfg.SnpEffDB == "" {
		missing = append(missing, "SnpEffDB")
	}
	if cfg.GeneDesc == "" {
		missing = append(missing, "GeneDesc")
	}
	if cfg.Prg == "" {
		missing = append(missing, "Prg")
	}
	if cfg.Gff == "" {
		missing = append(missing, "Gff")
	}
	return missing
}

func bsaseqType(cfg *utils.AnalysisConfig) (string, []int, error) {
	hp := cfg.HighParentName != "" && cfg.HighParentName != "None"
	lp := cfg.LowParentName != "" && cfg.LowParentName != "None"
	hb := cfg.HighBulkName != "" && cfg.HighBulkName != "None"
	lb := cfg.LowBulkName != "" && cfg.LowBulkName != "None"

	switch {
	case hp && lp && hb && lb:
		fmt.Println("================================ Running 2 parent 2 bulk analysis ==========================================")
		return "2p2b", []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}, nil

	case hp && lp && hb && !lb:
		fmt.Println("=================================== Running 2 parent High Bulk ============================================")
		return "2phb", []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx}, nil
	case hp && lp && !hb && lb:
		fmt.Println("=================================== Running 2 parent Low bulk =============================================")
		return "2plb", []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.LowBulkIdx}, nil
	case hp && !lp && hb && lb:
		fmt.Println("=================================== Running High parent 2 bulks ===========================================")
		return "hp2b", []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}, nil
	case hp && !lp && hb && !lb:
		fmt.Println("=================================== Running high parent high bulk ==========================================")
		return "hphb", []int{cfg.HighParentIdx, cfg.HighBulkIdx}, nil
	case hp && !lp && !hb && lb:
		fmt.Println("=================================== Running High parent low bulk ===========================================")
		return "hplb", []int{cfg.HighParentIdx, cfg.LowBulkIdx}, nil
	case !hp && lp && hb && lb:
		fmt.Println("=================================== Running Low parent 2 bulks =============================================")
		return "lp2b", []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}, nil
	case !hp && lp && hb && !lb:
		fmt.Println("Running Low parent high bulk")
		return "lphb", []int{cfg.LowParentIdx, cfg.HighBulkIdx}, nil
	case !hp && lp && !hb && lb:
		fmt.Println("Running Low parent low bulk")
		return "lplb", []int{cfg.LowParentIdx, cfg.LowBulkIdx}, nil
	case !hp && !lp && hb && lb:
		fmt.Println("Running bulks only")
		return "2b", []int{cfg.HighBulkIdx, cfg.LowBulkIdx}, nil
	default:
		return "bad", []int{}, fmt.Errorf("invalid combination — at least one bulk is required")
	}

}

func getRunType(cfg *utils.AnalysisConfig) (string,error) {

	if cfg.VCF != "" {
		if cfg.HighParBam == "" && cfg.LowParBam == "" && cfg.HighBulkBam == "" && cfg.LowBulkBam == ""  && cfg.HighParFwdReads == "" && cfg.HighParRevReads == "" && cfg.LowParFwdReads == "" && cfg.LowParRevReads == "" && cfg.HighBulkFwdReads == "" && cfg.HighBulkRevReads == "" && cfg.LowBulkFwdReads == "" && cfg.LowBulkRevReads == ""{
			return "vcf", nil
		}
		return "", fmt.Errorf("if VCF flag is passed then no other flag should be passed")
	}
	if cfg.HighBulkBam != "" || cfg.LowBulkBam != "" {
		if cfg.VCF == "" && cfg.HighParFwdReads == "" && cfg.HighParRevReads == "" && cfg.LowParFwdReads == "" && cfg.LowParRevReads == "" && cfg.HighBulkFwdReads == "" && cfg.HighBulkRevReads == "" && cfg.LowBulkFwdReads == "" && cfg.LowBulkRevReads == "" {
			return "bams", nil
		}
		return "", fmt.Errorf("if bulk bams flags are passed then neither VCFs or reads should be passed")
	}
	if cfg.HighParFwdReads != "" || cfg.HighParRevReads != "" || cfg.LowParFwdReads != "" || cfg.LowParRevReads != "" || cfg.HighBulkFwdReads != "" && cfg.HighBulkRevReads != "" && cfg.LowBulkFwdReads != "" && cfg.LowBulkRevReads != "" {
		if cfg.VCF == "" && cfg.HighParBam == "" && cfg.LowParBam == "" && cfg.HighBulkBam == "" && cfg.LowBulkBam == "" {
			return "reads", nil
		}
		return "", fmt.Errorf("if reads flags are passed then neither VCFs or bams should be passed")
	}

	return "", fmt.Errorf("invalid run type")
}

func bsaseq(cfg *utils.AnalysisConfig, hfcfg *utils.HardFilterConfig, btype string, idxs []int) error {
	//--------------------------------------- Filter -----------------------------------------------------------------//
	fmt.Printf("Filtering %s with bsaseq Type %v\n", cfg.VCF, btype)
	passedVariants, original, passed, err := filter.HardFilterVcf(*cfg, *hfcfg, btype, idxs)

	if err != nil {
		return err
	}
	color.Cyan("Original variants: %v\nFiltered Variants: %v", original, passed)

	// ----------------------------------------- Stats ---------------------------------------------------------------//
	rawStats, err := stats.RawStats(*cfg, btype, idxs, passedVariants)
	if err != nil {
		return err
	}

	// ---------------------------------------- smoothing ------------------------------------------------------------//
	smoothedStats, err := stats.SmoothAndNormalise(*cfg, btype, rawStats)
	if err != nil {
		return err
	}

	fmt.Println(len(smoothedStats))

	// ----------------------------------------- Threshold calculation -----------------------------------------------//
	thresholds, err := stats.CalculateThresholds(*cfg, btype, smoothedStats)
	if err != nil {
		return err
	}
	fmt.Println(len(thresholds))

	// ------------------------------------------ BRM blocks ------------------------------------------------------- //
	brmBlocks, err := stats.RunBRM(*cfg, btype, smoothedStats)
	if err != nil {
		return err
	}
	fmt.Printf("Detected %d BRM blocks\n", len(brmBlocks))

	// ------------------------------------------------ Plots ------------------------------------------------------- //
	if err := plots.GeneratePlots(*cfg, btype, smoothedStats, thresholds, brmBlocks); err != nil {
		fmt.Println("Error generating plots:", err)
		return err
	}

	// ---------------------------------------------- Detect QTLs --------------------------------------------------- //
	qtls, err := stats.DetectQTLs(*cfg, btype, smoothedStats)
	if err != nil {
		return err
	}
	fmt.Printf("Detected %d QTLs\n", len(qtls))

	// -------------------------------------------- Merge QTLs + BRM ------------------------------------------------ //
	merged, err := stats.MergeQTLsAndBRM(*cfg, btype, qtls, brmBlocks)
	if err != nil {
		return err
	}
	fmt.Printf("Merged intervals: %d\n", len(merged))

	// --------------------------------------------------- Gene Space ---------------------------------------------- //
	hardFilteredVcfPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.hardfiltered.vcf.gz", btype))
	if err := stats.GeneSpaceFromMerged(*cfg, btype, hardFilteredVcfPath, merged); err != nil {
		color.Red("Gene space analysis error: %v", err)
		// Non-fatal — log and continue.
	}

	return nil
}

func Run(cfg *utils.AnalysisConfig, hf utils.HardFilterConfig) error {
	// ============================================ Run Type ====================================================== //
	runType, err := getRunType(cfg)
	if err != nil {
		return err
	}
	fmt.Printf("Run type: %s\n", runType)

	if runType == "reads" {
		fmt.Printf("Running BSAseq with read-based analysis\n")
		return nil
	}
	if runType == "bams" {
		fmt.Printf("Running BSAseq with bam-based analysis\n")
		return nil
	}
	if runType == "vcf" {
		fmt.Printf("Running BSAseq with VCF-based analysis\n")
		return nil
	}


	fmt.Printf("HighParent: %s\n", cfg.HighParentName)
	fmt.Printf("LowParent: %s\n", cfg.LowParentName)
	fmt.Printf("HighBulk: %s\n", cfg.HighBulkName)
	fmt.Printf("LowBulk: %s\n", cfg.LowBulkName)



	// =========================================== Gene Space Check ================================================= //
	if missing := missingGeneSpaceParams(cfg); len(missing) > 0 {
		color.Yellow("Gene space analysis parameters missing: %s", strings.Join(missing, ", "))
		color.Blue("Continue without gene space analysis? [y/N]: ")
		var answer string
		if _, err := fmt.Scan(&answer); err != nil {
			return err
		}
		switch strings.ToLower(strings.TrimSpace(answer)) {
		case "y", "yes":
			color.Yellow("Continuing without gene space analysis.")
		default:
			return fmt.Errorf("missing gene space parameters: %s", strings.Join(missing, ", "))
		}
	}

	// ========================================== Open VCF ========================================================== //
	fmt.Println("Running BSAseq using a VCF file ...")

	bold := color.New(color.Bold).SprintFunc()
	color.Cyan("=============================== Checking parameters =====================================================\n")
	if !strings.HasSuffix(cfg.VCF, ".vcf.gz") && !strings.HasSuffix(cfg.VCF, ".vcf") {
		return fmt.Errorf("VCF file must be a .vcf or .vcf.gz file")
	}

	f, cleanup, err := openVCF(cfg.VCF)
	if err != nil {
		return fmt.Errorf("failed to open VCF %q: %w", cfg.VCF, err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(f, false)
	if err != nil {
		return fmt.Errorf("failed to create VCF reader: %w", err)
	}
	cfg.Rdr = rdr

	// ================================================ Sample Map ================================================== //
	sampleNames := rdr.Header.SampleNames
	sampleMap := map[int]string{0: "None"}
	for i, name := range sampleNames {
		sampleMap[i+1] = name
	}

	// ================================================= Sample selection =========================================== //
	color.Cyan("\n========================================== SAMPLE SELECTION =================================================\n\n")
	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	fmt.Println(sampleNames)

	var highParentChoice int
	var lowParentChoice int
	var highBulkChoice int
	var lowBulkChoice int

	//fmt.Printf("------------------------------------- SAMPLE CHOICES ----------------------------------------\n\n")

	// ------------------------------------------- High Parent ------------------------------------------------- //

	if cfg.HighParentName == "" {
		choice, err := promptForSampleChoice(sampleMap, "HIGH PARENT")
		if err != nil {
			color.Red("%v", err)
			return err
		}
		highParentChoice = choice

		cfg.HighParentName = sampleMap[highParentChoice]
		fmt.Printf("\n-----------------------------------------------------------\nHIGH Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighParentName))
		if highParentChoice != 0 {
			delete(sampleMap, highParentChoice)
		}
		cfg.HighParentIdx = highParentChoice - 1

	} else {
		fmt.Printf("HIGH parent is: %s \n\n", cfg.HighParentName)
		if !slices.Contains(sampleNames, cfg.HighParentName) {
			color.Yellow(" HIGH PARENT %s is not part of the VCF sample list\n", cfg.HighParentName)
			color.Cyan("Choose the number corresponding to the appropriate HIGH parent")
			choice, err := promptForSampleChoice(sampleMap, "HIGH PARENT")
			if err != nil {
				color.Red("%v", err)
				return err
			}
			highParentChoice = choice
			cfg.HighParentName = sampleMap[highParentChoice]
			if highParentChoice != 0 {
				delete(sampleMap, highParentChoice)
			}
			cfg.HighParentIdx = highParentChoice - 1
			fmt.Printf("\n-----------------------------------------------------------\nHIGH Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighParentName))
		} else if k, ok := sampleIndexByName(sampleMap, cfg.HighParentName); ok {
			cfg.HighParentIdx = k - 1
			delete(sampleMap, k)
		}
	}

	// ---------------------------------------------- Low Parent ------------------------------------------------- //

	if cfg.LowParentName == "" {
		choice, err := promptForSampleChoice(sampleMap, "LOW PARENT")
		if err != nil {
			color.Red("%v", err)
			return err
		}
		lowParentChoice = choice
		cfg.LowParentName = sampleMap[lowParentChoice]
		if lowParentChoice != 0 {
			delete(sampleMap, lowParentChoice)
		}
		cfg.LowParentIdx = lowParentChoice - 1
		fmt.Printf("\n-----------------------------------------------------------\nLOW Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowParentName))
	} else {
		fmt.Printf("LOW parent is: %s \n\n", cfg.LowParentName)

		if !slices.Contains(sampleNames, cfg.LowParentName) {
			color.Yellow("LOW PARENT %s is not part of the VCF sample list\n", cfg.LowParentName)
			color.Cyan("Choose the number corresponding to the appropriate LOW parent")
			choice, err := promptForSampleChoice(sampleMap, "LOW PARENT")
			if err != nil {
				color.Red("%v", err)
				return err
			}
			lowParentChoice = choice
			cfg.LowParentName = sampleMap[lowParentChoice]
			if lowParentChoice != 0 {
				delete(sampleMap, lowParentChoice)
			}
			cfg.LowParentIdx = lowParentChoice - 1
			fmt.Printf("\n-----------------------------------------------------------\nLOW Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowParentName))
		} else if k, ok := sampleIndexByName(sampleMap, cfg.LowParentName); ok {
			cfg.LowParentIdx = k - 1
			delete(sampleMap, k)
		}

	}

	// --------------------------------------------------- High Bulk ------------------------------------------------ //

	if cfg.HighBulkName == "" {
		color.Cyan("Choose the number corresponding to the appropriate HIGH BULK")
		choice, err := promptForSampleChoice(sampleMap, "HIGH BULK")
		if err != nil {
			color.Red("%v", err)
			return err
		}
		highBulkChoice = choice

		cfg.HighBulkName = sampleMap[highBulkChoice]
		fmt.Printf("\n-----------------------------------------------------------\nHIGH BULK is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighBulkName))
		if highBulkChoice != 0 {
			delete(sampleMap, highBulkChoice)
		}
		cfg.HighBulkIdx = highBulkChoice - 1
	} else {
		fmt.Printf("HIGH bulk is: %s \n\n", cfg.HighBulkName)
		if !slices.Contains(sampleNames, cfg.HighBulkName) {
			color.Yellow(" HIGH BULK %s is not part of the VCF sample list\n", cfg.HighBulkName)
			color.Cyan("Choose the number corresponding to the appropriate HIGH BULK")
			choice, err := promptForSampleChoice(sampleMap, "HIGH BULK")
			if err != nil {
				color.Red("%v", err)
				return err
			}
			highBulkChoice = choice
			cfg.HighBulkName = sampleMap[highBulkChoice]
			if highBulkChoice != 0 {
				delete(sampleMap, highBulkChoice)
			}
			cfg.HighBulkIdx = highBulkChoice - 1
		} else if k, ok := sampleIndexByName(sampleMap, cfg.HighBulkName); ok {
			cfg.HighBulkIdx = k - 1
			delete(sampleMap, k)
		}
	}

	// ---------------------------------------------- Low Bulk ------------------------------------------------ //

	if cfg.LowBulkName == "" {
		color.Cyan("Choose the number corresponding to the appropriate LOW BULK")
		choice, err := promptForSampleChoice(sampleMap, "LOW BULK")
		if err != nil {
			color.Red("%v", err)
			return err
		}
		lowBulkChoice = choice
		cfg.LowBulkName = sampleMap[lowBulkChoice]
		fmt.Printf("\n-----------------------------------------------------------\nLOW BULK is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowBulkName))
		if lowBulkChoice != 0 {
			delete(sampleMap, lowBulkChoice)
		}
		cfg.LowBulkIdx = lowBulkChoice - 1
	} else {
		fmt.Printf("LOW bulk is: %s \n\n", cfg.LowBulkName)
		if !slices.Contains(sampleNames, cfg.LowBulkName) {
			color.Yellow(" LOW BULK %s is not part of the VCF sample list\n", cfg.LowBulkName)
			color.Cyan("Choose the number corresponding to the appropriate LOW BULK")
			choice, err := promptForSampleChoice(sampleMap, "LOW BULK")
			if err != nil {
				color.Red("%v", err)
				return err
			}
			lowBulkChoice = choice
			cfg.LowBulkName = sampleMap[lowBulkChoice]
			if lowBulkChoice != 0 {
				delete(sampleMap, lowBulkChoice)
			}
			cfg.LowBulkIdx = lowBulkChoice - 1
		} else if k, ok := sampleIndexByName(sampleMap, cfg.LowBulkName); ok {
			cfg.LowBulkIdx = k - 1
			delete(sampleMap, k)
		}
	}

	// ================================================== Run ================================================== //

	bType, idxs, err := bsaseqType(cfg)
	if err != nil {
		return err
	}

	err = bsaseq(cfg, &hf, bType, idxs)
	if err != nil {
		return err
	}

	return nil
}
