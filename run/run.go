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

// ─── Analysis category ────────────────────────────────────────────────────────

type analysisCategory int

const (
	TwoParentsTwoBulks analysisCategory = iota
	TwoParentsHighBulk
	TwoParentsLowBulk
	HighParentTwoBulks
	HighParentHighBulk
	HighParentLowBulk
	LowParentTwoBulks
	LowParentHighBulk
	LowParentLowBulk
	BulksOnly
)

func (a analysisCategory) String() string {
	return [...]string{
		"Two Parents, Two Bulks",
		"Two Parents, High Bulk",
		"Two Parents, Low Bulk",
		"High Parent, Two Bulks",
		"High Parent, High Bulk",
		"High Parent, Low Bulk",
		"Low Parent, Two Bulks",
		"Low Parent, High Bulk",
		"Low Parent, Low Bulk",
		"Bulks Only",
	}[a]
}

func categorise(cfg *utils.AnalysisConfig) (analysisCategory, error) {
	hp := cfg.HighParentName != "" && cfg.HighParentName != "None"
	lp := cfg.LowParentName != "" && cfg.LowParentName != "None"
	hb := cfg.HighBulkName != "" && cfg.HighBulkName != "None"
	lb := cfg.LowBulkName != "" && cfg.LowBulkName != "None"

	if !hb && !lb {
		return 0, fmt.Errorf("at least one bulk must be selected")
	}

	switch {
	case hp && lp && hb && lb:
		return TwoParentsTwoBulks, nil
	case hp && lp && hb && !lb:
		return TwoParentsHighBulk, nil
	case hp && lp && !hb && lb:
		return TwoParentsLowBulk, nil
	case hp && !lp && hb && lb:
		return HighParentTwoBulks, nil
	case hp && !lp && hb && !lb:
		return HighParentHighBulk, nil
	case hp && !lp && !hb && lb:
		return HighParentLowBulk, nil
	case !hp && lp && hb && lb:
		return LowParentTwoBulks, nil
	case !hp && lp && hb && !lb:
		return LowParentHighBulk, nil
	case !hp && lp && !hb && lb:
		return LowParentLowBulk, nil
	case !hp && !lp && hb && lb:
		return BulksOnly, nil
	default:
		return 0, fmt.Errorf("invalid combination — at least one bulk is required")
	}
}

// ─── Run ──────────────────────────────────────────────────────────────────────

func Run(cfg *utils.AnalysisConfig, hf utils.HardFilterConfig) error {

	// ── Gene space check ────────────────────────────────────────────────────
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

	// ── VCF setup ───────────────────────────────────────────────────────────
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

	// ── Sample map ──────────────────────────────────────────────────────────
	sampleNames := rdr.Header.SampleNames
	sampleMap := map[int]string{0: "None"}
	for i, name := range sampleNames {
		sampleMap[i+1] = name
	}

	color.Cyan("\n========================================== SAMPLE SELECTION =================================================\n\n")
	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	fmt.Println(sampleNames)

	var highParentChoice int
	var lowParentChoice int
	var highBulkChoice int
	var lowBulkChoice int

	// ── PARENTS ─────────────────────────────────────────────────────────────
	fmt.Printf("------------------------------------- PARENT CHOICES ----------------------------------------\n\n")

	if cfg.HighParentName == "" && cfg.LowParentName == "" {
		fmt.Printf("No parent samples specified ...\n\n")
		color.Cyan("Enter number corresponding to the sample ...\n\n")

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

		keys = slices.Sorted(maps.Keys(sampleMap))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleMap[i])
		}

		color.Blue("\nEnter LOW PARENT number:")
		_, err = fmt.Scan(&lowParentChoice)
		if err != nil {
			fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", err)
			return err
		}
		if lowParentChoice == highParentChoice && lowParentChoice != 0 {
			fmt.Println("LOW PARENT should not be the same as HIGH PARENT")
			return fmt.Errorf("invalid input")
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

	} else if cfg.HighParentName != "" && cfg.LowParentName != "" {
		fmt.Printf("HIGH parent is: %s \n\n", cfg.HighParentName)
		fmt.Printf("LOW parent is: %s \n\n", cfg.LowParentName)

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

		if !slices.Contains(sampleNames, cfg.LowParentName) {
			fmt.Printf(" LOW PARENT %s is not part of the VCF sample list\n", cfg.LowParentName)
			color.Cyan("Choose the number corresponding to the appropriate LOW PARENT")
			keys := slices.Sorted(maps.Keys(sampleMap))
			for _, i := range keys {
				fmt.Printf("%v : %v\n", i, sampleMap[i])
			}
			fmt.Println("Enter LOW PARENT number:")
			_, err := fmt.Scan(&lowParentChoice)
			if err != nil {
				fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", err)
				return err
			}
			if lowParentChoice == highParentChoice && lowParentChoice != 0 {
				fmt.Println("LOW PARENT should not be the same as HIGH PARENT")
				return fmt.Errorf("invalid input")
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

	// ── BULKS ───────────────────────────────────────────────────────────────
	fmt.Printf("\n====================================================================================================================================\n\n")
	fmt.Printf("------------------------------------- BULK CHOICES ----------------------------------------\n\n")

	if cfg.HighBulkName == "" && cfg.LowBulkName == "" {
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
		if highBulkChoice == highParentChoice || highBulkChoice == lowParentChoice {
			fmt.Println("Your HIGH bulk cannot be the same as any of the parents")
			return fmt.Errorf("invalid input")
		}
		cfg.HighBulkName = sampleMap[highBulkChoice]
		fmt.Printf("\n-----------------------------------------------------------\nHIGH BULK is: %s\n-----------------------------------------------------------\n\n", bold(cfg.HighBulkName))
		if highBulkChoice != 0 {
			delete(sampleMap, highBulkChoice)
		}
		cfg.HighBulkIdx = highBulkChoice - 1

		keys = slices.Sorted(maps.Keys(sampleMap))
		for _, i := range keys {
			fmt.Printf("%v : %v\n", i, sampleMap[i])
		}

		color.Blue("Enter LOW BULK number:")
		_, err = fmt.Scan(&lowBulkChoice)
		if err != nil {
			fmt.Printf("LOW BULK number should be numerical and part of the list above: %s\n", err)
			return fmt.Errorf("invalid input")
		}
		if lowBulkChoice == highBulkChoice || lowBulkChoice == highParentChoice || lowBulkChoice == lowParentChoice {
			fmt.Println("Your LOW bulk cannot be the same as any of the parents OR the HIGH bulk")
			return fmt.Errorf("invalid input")
		}
		cfg.LowBulkName = sampleMap[lowBulkChoice]
		fmt.Printf("\n-----------------------------------------------------------\nLOW BULK is: %s\n-----------------------------------------------------------\n\n", bold(cfg.LowBulkName))
		if lowBulkChoice != 0 {
			delete(sampleMap, lowBulkChoice)
		}
		cfg.LowBulkIdx = lowBulkChoice - 1
	}

	// If bulks were provided via flags, resolve their indices
	if cfg.HighBulkName != "" && cfg.LowBulkName != "" {
		if slices.Contains(sampleNames, cfg.HighBulkName) {
			for k, v := range sampleMap {
				if v == cfg.HighBulkName {
					cfg.HighBulkIdx = k - 1
					delete(sampleMap, k)
					break
				}
			}
		} else {
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
		}

		if slices.Contains(sampleNames, cfg.LowBulkName) {
			for k, v := range sampleMap {
				if v == cfg.LowBulkName {
					cfg.LowBulkIdx = k - 1
					delete(sampleMap, k)
					break
				}
			}
		} else {
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
			if !slices.Contains(keys, lowBulkChoice) {
				color.Red("The selected number is not in the list.")
				return fmt.Errorf("invalid input")
			}
			cfg.LowBulkName = sampleMap[lowBulkChoice]
			if lowBulkChoice != 0 {
				delete(sampleMap, lowBulkChoice)
			}
			cfg.LowBulkIdx = lowBulkChoice - 1
		}
	}

	// ── Categorise & summary ────────────────────────────────────────────────
	cat, err := categorise(cfg)
	if err != nil {
		color.Red("%s", err)
		return err
	}

	// ── Output directory ────────────────────────────────────────────────────
	if cfg.OutputDir != "." {
		if err := os.MkdirAll(cfg.OutputDir, 0755); err != nil {
			return fmt.Errorf("failed to create output directory: %w", err)
		}
	}

	// ── Dispatch ────────────────────────────────────────────────────────────
	switch cat {
	case TwoParentsTwoBulks:
		//fmt.Println("Working with two bulks")
		color.Green("=================================== Running Two Parents 2 Bulk Analysis =============================================\n\n")
		fmt.Printf("High Parent: %s, Index: %v\n", cfg.HighParentName, cfg.HighParentIdx)
		fmt.Printf("Low Parent: %s, Index: %v\n", cfg.LowParentName, cfg.LowParentIdx)
		fmt.Printf("High Bulk: %s, Index: %v\n", cfg.HighBulkName, cfg.HighBulkIdx)
		fmt.Printf("Low Bulk: %s, Index: %v\n", cfg.LowBulkName, cfg.LowBulkIdx)
		// twobulk.RunTwoBulkTwoParents(cfg, hf)
	case TwoParentsHighBulk:
		fmt.Println("Working with one bulk BSAseq (HIGH bulk)...")
		// onebulk.RunTwoParentsHighBulk(cfg, hf)
	case TwoParentsLowBulk:
		fmt.Println("Working with one bulk BSAseq (LOW bulk)")
		// onebulk.RunTwoParentsLowBulk(cfg, hf)
	case HighParentTwoBulks:
		// onebulk.RunHighParentTwoBulks(cfg, hf)
	case HighParentHighBulk:
		// onebulk.RunHighParentHighBulk(cfg, hf)
	case HighParentLowBulk:
		// onebulk.RunHighParentLowBulk(cfg, hf)
	case LowParentTwoBulks:
		// onebulk.RunLowParentTwoBulks(cfg, hf)
	case LowParentHighBulk:
		// onebulk.RunLowParentHighBulk(cfg, hf)
	case LowParentLowBulk:
		// onebulk.RunLowParentLowBulk(cfg, hf)
	case BulksOnly:
		fmt.Println("Running bulks only")
		fmt.Printf("High Bulk: %s, Index: %v\n", cfg.HighBulkName, cfg.HighBulkIdx)
		fmt.Printf("Low Bulk: %s, Index: %v\n", cfg.LowBulkName, cfg.LowBulkIdx)
		// twobulk.RunBulksOnly(cfg, hf)
	}

	return nil
}
