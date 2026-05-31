package run

import (
	"compress/gzip"
	"fmt"
	"io"
	"maps"
	"os"
	"slices"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/crefactor/onebulk"
	"github.com/gmaffy/GoBSAseq/crefactor/twobulk"
	"github.com/gmaffy/GoBSAseq/crefactor/utils"
)

func openVCF(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}

	cleanup := func() { f.Close() }

	// Check suffix
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, err
		}

		cleanup = func() {
			gz.Close()
			f.Close()
		}

		return gz, cleanup, nil
	}

	// Plain text VCF
	return f, cleanup, nil
}

func missingGeneSpaceParams(cfg utils.AnalysisConfig) []string {
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

func Run(cfg utils.AnalysisConfig, hfCfg utils.HardFilterConfig) error { //, vcf string, highParentDepth int, lowParentDepth int, oneParentDepth int, highBulkDepth int, lowBulkDepth int, oneBulkDepth int, highBulkSize int, lowBulkSize int, oneBulkSize int, windowSize int, population string, recurrent bool, rep int, alpha float64, minQTL int64, mergeDist int64, outputDir string) error {
	bold := color.New(color.Bold).SprintFunc()
	color.Cyan("=============================== Checking parameters =====================================================\n")

	fmt.Printf("VCF: %s\n", cfg.VCF)
	fmt.Printf("Min High Parent Depth: %d\n", cfg.HighParentDepth)
	fmt.Printf("Min Low Parent Depth: %d\n", cfg.LowParentDepth)
	fmt.Printf("One Min Parent Depth: %d\n", cfg.OneParentDepth)
	fmt.Printf("High Bulk Depth: %d\n", cfg.HighBulkDepth)
	fmt.Printf("Low Bulk Depth: %d\n", cfg.LowBulkDepth)
	fmt.Printf("One Bulk Depth: %d\n", cfg.OneBulkDepth)
	fmt.Printf("High Bulk Size: %d\n", cfg.HighBulkSize)
	fmt.Printf("Low Bulk Size: %d\n", cfg.LowBulkSize)
	fmt.Printf("One Bulk Size: %d\n", cfg.OneBulkSize)
	fmt.Printf("Window Size: %d\n", cfg.WindowSize)
	fmt.Printf("Step Size: %d\n", cfg.StepSize)
	fmt.Printf("Population: %s\n", cfg.Population)
	//fmt.Printf("Recurrent: %v\n")
	fmt.Printf("Simulations: %d\n", cfg.Rep)
	fmt.Printf("Alphas: %v\n", cfg.Alphas)
	fmt.Printf("Min QTL Length: %d\n", cfg.MinQTLWidth)
	fmt.Printf("Merge Distance: %d\n", cfg.MergeDistance)
	fmt.Printf("Output Dir: %s\n", cfg.OutputDir)

	if !strings.HasSuffix(cfg.VCF, ".vcf.gz") && !strings.HasSuffix(cfg.VCF, ".vcf") {
		color.Red("VCF file must be a .vcf or .vcf.gz file")
		return fmt.Errorf("VCF file must be a .vcf or .vcf.gz file")
	}

	f, cleanup, err := openVCF(cfg.VCF)
	if err != nil {
		panic(err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(f, false)
	if err != nil {
		panic(err)
	}

	cfg.Rdr = rdr

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
			color.Red("Analysis cancelled. Re-run with the missing gene space parameters.")
			return fmt.Errorf("missing gene space parameters: %s", strings.Join(missing, ", "))
		}
	}

	var highParentChoice int
	var lowParentChoice int
	var highBulkChoice int
	var lowBulkChoice int
	//var oneParentChoice int

	color.Cyan("\n========================================== SAMPLE SELECTION =================================================\n\n")
	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	sampleNames := rdr.Header.SampleNames
	sampleNamesDic := make(map[int]string)
	sampleNamesDic[0] = "None"
	for i, _ := range sampleNames {
		sampleNamesDic[i+1] = sampleNames[i]
	}
	sampleIndex := func(name string) (int, bool) {
		for i, sampleName := range sampleNames {
			if sampleName == name {
				return i, true
			}
		}
		return -1, false
	}
	resolveNamedSample := func(label, name string) (int, error) {
		idx, ok := sampleIndex(name)
		if !ok {
			return -1, fmt.Errorf("%s %q is not part of the VCF sample list", label, name)
		}
		delete(sampleNamesDic, idx+1)
		return idx, nil
	}
	parentsProvided := cfg.HighParentName != "" || cfg.LowParentName != "" || cfg.OneParentName != ""
	bulksOnlyFromCLI := !parentsProvided && cfg.HighBulkName != "" && cfg.LowBulkName != "" && cfg.OneBulkName == ""

	fmt.Printf("------------------------------------- PARENT CHOICES ----------------------------------------\n\n")
	if bulksOnlyFromCLI {
		color.Yellow("No parent samples supplied; running in bulks-only mode from CLI bulk names.")
	} else if cfg.HighParentName == "" && cfg.LowParentName == "" && cfg.OneParentName == "" {
		fmt.Printf("No parent samples specified ...\n\n")
		color.Cyan("Enter number corresponding to the sample ...\n\n")

		keys := slices.Sorted(maps.Keys(sampleNamesDic))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleNamesDic[i])
		}

		color.Blue("\nEnter HIGH PARENT number:")
		_, highParErr := fmt.Scan(&highParentChoice)
		if highParErr != nil {
			color.Red("HIGH PARENT number should be numerical and part of the list above: %s\n", highParErr)
			return highParErr
		}

		if !slices.Contains(keys, highParentChoice) {
			color.Red("The selected number is not in the list.")
			return fmt.Errorf("invalid input")
		}

		cfg.HighParentName = sampleNamesDic[highParentChoice]
		fmt.Printf("\n-----------------------------------------------------------\nHIGH Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighParentName))
		if highParentChoice != 0 {
			delete(sampleNamesDic, highParentChoice)
		}
		cfg.HighParentIdx = highParentChoice - 1

		keys = slices.Sorted(maps.Keys(sampleNamesDic))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleNamesDic[i])
		}

		color.Blue("\nEnter LOW PARENT number:")
		_, lowParErr := fmt.Scan(&lowParentChoice)
		if lowParErr != nil {
			fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", lowParErr)
			return lowParErr
		}

		if lowParentChoice == highParentChoice && lowParentChoice != 0 {
			fmt.Println("LOW PARENT should not be the same as HIGH PARENT")
			return fmt.Errorf("invalid input")
		}

		if !slices.Contains(keys, lowParentChoice) {
			color.Red("The selected number is not in the list.")
			return fmt.Errorf("invalid input")
		}

		cfg.LowParentName = sampleNamesDic[lowParentChoice]

		if lowParentChoice != 0 {
			delete(sampleNamesDic, lowParentChoice)
		}
		cfg.LowParentIdx = lowParentChoice - 1
		fmt.Printf("\n-----------------------------------------------------------\nLOW Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowParentName))

		// ============================= PARENTS ARGUMENTS PASSED ======================================================= //
	} else if cfg.HighParentName != "" && cfg.LowParentName != "" {
		fmt.Printf("HIGH parent is: %s \n\n", cfg.HighParentName)
		fmt.Printf("LOW parent is: %s \n\n", cfg.LowParentName)

		if !slices.Contains(sampleNames, cfg.HighParentName) {
			color.Yellow(" HIGH PARENT %s is not part of the VCF sample list\n", cfg.HighParentName)
			color.Cyan("Choose the number corresponding to the appropriate HIGH parent")
			keys := slices.Sorted(maps.Keys(sampleNamesDic))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleNamesDic[i])
			}
			_, highParErr := fmt.Scan(&highParentChoice)
			if highParErr != nil {
				fmt.Printf("HIGH PARENT number should be numerical and part of the list above: %s\n", highParErr)
				return highParErr
			}

			if !slices.Contains(keys, highParentChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}

			cfg.HighParentName = sampleNamesDic[highParentChoice]
			if highParentChoice != 0 {
				delete(sampleNamesDic, highParentChoice)
			}
			cfg.HighParentIdx = highParentChoice - 1
			fmt.Printf("\n-----------------------------------------------------------\nHIGH Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighParentName))

		} else {
			idx, err := resolveNamedSample("HIGH PARENT", cfg.HighParentName)
			if err != nil {
				return err
			}
			cfg.HighParentIdx = idx

		}

		if !slices.Contains(sampleNames, cfg.LowParentName) {
			fmt.Printf(" LOW PARENT %s is not part of the VCF sample list\n", cfg.LowParentName)
			color.Cyan("Choose the number corresponding to the appropriate LOW PARENT")
			keys := slices.Sorted(maps.Keys(sampleNamesDic))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleNamesDic[i])
			}
			fmt.Println("Enter LOW PARENT number:")
			_, lowParErr := fmt.Scan(&lowParentChoice)
			if lowParErr != nil {
				fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", lowParErr)
				return lowParErr
			}

			if lowParentChoice == highParentChoice && lowParentChoice != 0 {
				fmt.Println("LOW PARENT should not be the same as HIGH PARENT")
				return fmt.Errorf("invalid input")
			}

			if !slices.Contains(keys, lowParentChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}

			cfg.LowParentName = sampleNamesDic[lowParentChoice]
			fmt.Printf("LOW parent is: %s \n\n", cfg.LowParentName)
			if lowParentChoice != 0 {
				delete(sampleNamesDic, lowParentChoice)
			}
			cfg.LowParentIdx = lowParentChoice - 1
			fmt.Printf("\n-----------------------------------------------------------\nLOW Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowParentName))
		} else {
			idx, err := resolveNamedSample("LOW PARENT", cfg.LowParentName)
			if err != nil {
				return err
			}
			cfg.LowParentIdx = idx

		}

	}

	fmt.Printf("\n====================================================================================================================================\n\n")
	fmt.Printf("------------------------------------- BULK CHOICES ----------------------------------------\n\n")

	if cfg.HighBulkName == "" && cfg.LowBulkName == "" && cfg.OneBulkName == "" {
		color.Cyan("Choose the number corresponding to the appropriate HIGH BULK")
		keys := slices.Sorted(maps.Keys(sampleNamesDic))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleNamesDic[i])
		}
		color.Blue("\nEnter HIGH BULK number:")
		_, highBulkErr := fmt.Scan(&highBulkChoice)
		if highBulkErr != nil {
			fmt.Printf("HIGH BULK number should be numerical and part of the list above: %s\n", highBulkErr)
			return fmt.Errorf("invalid input")
		}

		if highBulkChoice == highParentChoice || highBulkChoice == lowParentChoice {
			fmt.Println("Your HIGH bulk cannot be the same as any of the parents")
			return fmt.Errorf("invalid input")
		}

		cfg.HighBulkName = sampleNamesDic[highBulkChoice]
		fmt.Printf("\n-----------------------------------------------------------\nHIGH BULK is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighBulkName))
		if highBulkChoice != 0 {
			delete(sampleNamesDic, highBulkChoice)
		}
		cfg.HighBulkIdx = highBulkChoice - 1

		keys = slices.Sorted(maps.Keys(sampleNamesDic))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleNamesDic[i])
		}

		color.Blue("Enter LOW BULK number:")
		_, lowBulkErr := fmt.Scan(&lowBulkChoice)
		if lowBulkErr != nil {
			fmt.Printf("LOW BULK number should be numerical and part of the list above: %s\n", lowBulkErr)
			return fmt.Errorf("invalid input")
		}

		// i dont think we will ever get here with the choice deletes
		if lowBulkChoice == highBulkChoice || lowBulkChoice == highParentChoice || lowBulkChoice == lowParentChoice {
			fmt.Println("Your LOW bulk cannot be the same as any of the parents OR the HIGH bulk")
			return fmt.Errorf("invalid input")
		}
		cfg.LowBulkName = sampleNamesDic[lowBulkChoice]
		fmt.Printf("\n-----------------------------------------------------------\nLOW BULK is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowBulkName))
		if lowBulkChoice != 0 {
			delete(sampleNamesDic, lowBulkChoice)
		}
		cfg.LowBulkIdx = lowBulkChoice - 1

	} else {
		// CLI path. A single --bulks value is treated as a low/affected bulk
		// for compatibility; use "--bulks HighBulk," to force high-bulk mode.
		if cfg.OneBulkName != "" && cfg.HighBulkName == "" && cfg.LowBulkName == "" {
			cfg.LowBulkName = cfg.OneBulkName
			if cfg.OneBulkDepth > 0 {
				cfg.LowBulkDepth = cfg.OneBulkDepth
			}
			if cfg.OneBulkSize > 0 {
				cfg.LowBulkSize = cfg.OneBulkSize
			}
			color.Yellow("Single --bulks value %q treated as one-bulk LOW/affected sample.", cfg.OneBulkName)
		}

		if cfg.HighBulkName != "" {
			idx, err := resolveNamedSample("HIGH BULK", cfg.HighBulkName)
			if err != nil {
				return err
			}
			cfg.HighBulkIdx = idx
		}
		if cfg.LowBulkName != "" {
			idx, err := resolveNamedSample("LOW BULK", cfg.LowBulkName)
			if err != nil {
				return err
			}
			cfg.LowBulkIdx = idx
		}
		if cfg.OneBulkName != "" {
			cfg.OneBulkIdx = cfg.LowBulkIdx
			if cfg.HighBulkName != "" && cfg.LowBulkName == "" {
				cfg.OneBulkIdx = cfg.HighBulkIdx
			}
		}
	}

	runDir, err := utils.PrepareRunDir(cfg.OutputDir)
	if err != nil {
		return err
	}
	cfg.OutputDir = runDir
	color.Green("Results directory: %s", runDir)

	mode := resolveAnalysisMode(cfg, highParentChoice, lowParentChoice, highBulkChoice, lowBulkChoice)
	switch mode {
	case modeBulksOnly:
		fmt.Println("Running bulks only")
		fmt.Printf("High Bulk: %s, Index: %v\n", cfg.HighBulkName, cfg.HighBulkIdx)
		fmt.Printf("Low Bulk: %s, Index: %v\n", cfg.LowBulkName, cfg.LowBulkIdx)
		twobulk.RunTwoBulksOnly(cfg, hfCfg)
	case modeOneBulkHigh:
		fmt.Println("Working with one bulk BSAseq (HIGH bulk)...")
		if err := onebulk.RunTwoParentsHighBulk(cfg, hfCfg); err != nil {
			color.Red("Error running one bulk analysis: %v", err)
			return err
		}
	case modeOneBulkLow:
		fmt.Println("Working with one bulk BSAseq (LOW bulk)")
		if err := onebulk.RunTwoParentsLowBulk(cfg, hfCfg); err != nil {
			color.Red("Error running one bulk analysis: %v", err)
			return err
		}
	case modeTwoBulkTwoParents:
		fmt.Println("Working with two bulks")
		color.Green("=================================== Running Two Bulk Analysis =============================================\n\n")
		fmt.Printf("High Parent: %s, Index: %v\n", cfg.HighParentName, cfg.HighParentIdx)
		fmt.Printf("Low Parent: %s, Index: %v\n", cfg.LowParentName, cfg.LowParentIdx)
		fmt.Printf("High Bulk: %s, Index: %v\n", cfg.HighBulkName, cfg.HighBulkIdx)
		fmt.Printf("Low Bulk: %s, Index: %v\n", cfg.LowBulkName, cfg.LowBulkIdx)
		twobulk.RunTwoBulkTwoParents(cfg, hfCfg)
	default:
		return fmt.Errorf("unsupported sample layout: need two bulks, or two parents plus one bulk")
	}
	return nil
}

type analysisMode int

const (
	modeBulksOnly analysisMode = iota
	modeOneBulkHigh
	modeOneBulkLow
	modeTwoBulkTwoParents
	modeUnsupported
)

func resolveAnalysisMode(cfg utils.AnalysisConfig, highPar, lowPar, highBulk, lowBulk int) analysisMode {
	hasHighPar := cfg.HighParentName != "" || highPar > 0
	hasLowPar := cfg.LowParentName != "" || lowPar > 0
	hasHighBulk := cfg.HighBulkName != "" || highBulk > 0
	hasLowBulk := cfg.LowBulkName != "" || lowBulk > 0
	hasOneBulk := cfg.OneBulkName != ""

	if hasHighBulk && hasLowBulk && !hasHighPar && !hasLowPar {
		return modeBulksOnly
	}
	if hasHighPar && hasLowPar && hasOneBulk && !hasHighBulk && !hasLowBulk {
		return modeOneBulkLow
	}
	if hasHighPar && hasLowPar && hasHighBulk && !hasLowBulk {
		return modeOneBulkHigh
	}
	if hasHighPar && hasLowPar && hasLowBulk && !hasHighBulk {
		return modeOneBulkLow
	}
	if hasHighPar && hasLowPar && hasHighBulk && hasLowBulk {
		return modeTwoBulkTwoParents
	}
	// Interactive fallbacks from numeric choices
	if lowBulk != 0 && highBulk != 0 && highPar == 0 && lowPar == 0 {
		return modeBulksOnly
	}
	if lowBulk == 0 && highBulk != 0 && lowPar != 0 && highPar != 0 {
		return modeOneBulkHigh
	}
	if highBulk == 0 && highPar != 0 && lowPar != 0 {
		return modeOneBulkLow
	}
	return modeUnsupported
}
