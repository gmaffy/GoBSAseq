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
	"github.com/gmaffy/GoBSAseq/filter"
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
		fmt.Println("=================================== Running high parent high bulk ===========================================")
		return "hphb", []int{cfg.HighParentIdx, cfg.HighBulkIdx}, nil
	case hp && !lp && !hb && lb:
		fmt.Println("Running High parent low bulk")
		return "hplb", []int{cfg.HighParentIdx, cfg.LowBulkIdx}, nil
	case !hp && lp && hb && lb:
		fmt.Println("Running Low parent 2 bulks")
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

func bsaseq(cfg *utils.AnalysisConfig, hfcfg *utils.HardFilterConfig, btype string, idxs []int) error {
	//--------------------------------------- Filter -----------------------------------------------------------------//
	fmt.Printf("Filtering %s with bsaseq Type %v\n", cfg.VCF, btype)
	passedVariants, original, passed, err := filter.HardFilterVcf(*cfg, *hfcfg, btype, idxs)

	if err != nil {
		return err
	}
	color.Cyan("Original variants: %v\nFiltered Variants: %v", original, passed)

	// ----------------------------------------- Stats ---------------------------------------------------------------//
	if _, err := stats.RawStats(*cfg, btype, idxs, passedVariants); err != nil {
		return err
	}

	// ---------------------------------------- smoothing -----------------------------------------------------------//

	return nil
}

func Run(cfg *utils.AnalysisConfig, hf utils.HardFilterConfig) error {

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
		keys := slices.Sorted(maps.Keys(sampleMap))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleMap[i])
		}

		color.Blue("\nEnter HIGH PARENT number:")
		_, err := fmt.Scan(&highParentChoice)
		if err != nil {
			color.Red("HIGH PARENT number should be numerical and part of the list above: %s\n", err)
			return err
		}
		if !slices.Contains(keys, highParentChoice) {
			color.Red("The selected number is not in the list.")
			return fmt.Errorf("invalid input")
		}

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
			keys := slices.Sorted(maps.Keys(sampleMap))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleMap[i])
			}
			_, err := fmt.Scan(&highParentChoice)
			if err != nil {
				fmt.Printf("HIGH PARENT number should be numerical and part of the list above: %s\n", err)
				return err
			}
			if !slices.Contains(keys, highParentChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}
			cfg.HighParentName = sampleMap[highParentChoice]
			if highParentChoice != 0 {
				delete(sampleMap, highParentChoice)
			}
			cfg.HighParentIdx = highParentChoice - 1
			fmt.Printf("\n-----------------------------------------------------------\nHIGH Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighParentName))
		} else {
			for k, v := range sampleMap {
				if v == cfg.HighParentName {
					cfg.HighParentIdx = k - 1
					delete(sampleMap, k)
					break
				}
			}
		}
	}

	// ---------------------------------------------- Low Parent ------------------------------------------------- //

	if cfg.LowParentName == "" {
		keys := slices.Sorted(maps.Keys(sampleMap))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleMap[i])
		}

		color.Blue("\nEnter LOW PARENT number:")
		_, err := fmt.Scan(&lowParentChoice)
		if err != nil {
			fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", err)
			return err
		}
		keys = slices.Sorted(maps.Keys(sampleMap))
		if !slices.Contains(keys, lowParentChoice) {
			color.Red("The selected number is not in the list.")
			return fmt.Errorf("invalid input")
		}
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
			keys := slices.Sorted(maps.Keys(sampleMap))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleMap[i])
			}
			_, err := fmt.Scan(&lowParentChoice)
			if err != nil {
				fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", err)
				return err
			}
			if !slices.Contains(keys, lowParentChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}
			cfg.LowParentName = sampleMap[lowParentChoice]
			if lowParentChoice != 0 {
				delete(sampleMap, lowParentChoice)
			}
			cfg.LowParentIdx = lowParentChoice - 1
			fmt.Printf("\n-----------------------------------------------------------\nLOW Parent is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowParentName))
		} else {
			for k, v := range sampleMap {
				if v == cfg.LowParentName {
					cfg.LowParentIdx = k - 1
					delete(sampleMap, k)
					break
				}
			}
		}

	}

	// --------------------------------------------------- High Bulk ------------------------------------------------ //

	if cfg.HighBulkName == "" {
		color.Cyan("Choose the number corresponding to the appropriate HIGH BULK")
		keys := slices.Sorted(maps.Keys(sampleMap))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleMap[i])
		}
		color.Blue("\nEnter HIGH BULK number:")
		_, err := fmt.Scan(&highBulkChoice)
		if err != nil {
			fmt.Printf("HIGH BULK number should be numerical and part of the list above: %s\n", err)
			return fmt.Errorf("invalid input")
		}

		if !slices.Contains(keys, highBulkChoice) {
			color.Red("The selected number is not in the list.")
			return fmt.Errorf("invalid input")
		}

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
			keys := slices.Sorted(maps.Keys(sampleMap))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleMap[i])
			}
			_, err := fmt.Scan(&highBulkChoice)
			if err != nil {
				return fmt.Errorf("HIGH BULK number should be numerical and part of the list above: %w", err)
			}
			if !slices.Contains(keys, highBulkChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}
			cfg.HighBulkName = sampleMap[highBulkChoice]
			if highBulkChoice != 0 {
				delete(sampleMap, highBulkChoice)
			}
			cfg.HighBulkIdx = highBulkChoice - 1
		} else {
			for k, v := range sampleMap {
				if v == cfg.HighBulkName {
					cfg.HighBulkIdx = k - 1
					delete(sampleMap, k)
					break
				}
			}
		}
	}

	// ---------------------------------------------- Low Bulk ------------------------------------------------ //

	if cfg.LowBulkName == "" {
		color.Cyan("Choose the number corresponding to the appropriate LOW BULK")
		keys := slices.Sorted(maps.Keys(sampleMap))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleMap[i])
		}
		color.Blue("Enter LOW BULK number:")
		_, err := fmt.Scan(&lowBulkChoice)
		if err != nil {
			fmt.Printf("LOW BULK number should be numerical and part of the list above: %s\n", err)
			return fmt.Errorf("invalid input")
		}
		if lowBulkChoice == highBulkChoice || lowBulkChoice == highParentChoice || lowBulkChoice == lowParentChoice {
			fmt.Println("Your LOW bulk cannot be the same as any of the parents OR the HIGH bulk")
			return fmt.Errorf("invalid input")
		}
		if !slices.Contains(keys, lowBulkChoice) {
			color.Red("The selected number is not in the list.")
			return fmt.Errorf("invalid input")
		}
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
			keys := slices.Sorted(maps.Keys(sampleMap))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleMap[i])
			}
			_, err := fmt.Scan(&lowBulkChoice)
			if err != nil {
				return fmt.Errorf("LOW BULK number should be numerical and part of the list above: %w", err)
			}
			if lowBulkChoice == highBulkChoice || lowBulkChoice == highParentChoice || lowBulkChoice == lowParentChoice {
				fmt.Println("Your LOW bulk cannot be the same as any of the parents OR the HIGH bulk")
				return fmt.Errorf("invalid input")
			}
			if !slices.Contains(keys, lowBulkChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}
			cfg.LowBulkName = sampleMap[lowBulkChoice]
			if lowBulkChoice != 0 {
				delete(sampleMap, lowBulkChoice)
			}
			cfg.LowBulkIdx = lowBulkChoice - 1
		} else {
			for k, v := range sampleMap {
				if v == cfg.LowBulkName {
					cfg.LowBulkIdx = k - 1
					delete(sampleMap, k)
					break
				}
			}
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
