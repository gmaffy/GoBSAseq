package run

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/twobulk"
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

func Run(vcf string, highParentDepth int, lowParentDepth int, oneParentDepth int, highBulkDepth int, lowBulkDepth int, oneBulkDepth int, highBulkSize int, lowBulkSize int, oneBulkSize int, windowSize int, population string, recurrent bool, rep int, alpha float64, minQTL int64, mergeDist int64, outputDir string) error {
	fmt.Printf("VCF: %s\n", vcf)
	fmt.Printf("High Parent Depth: %d\n", highParentDepth)
	fmt.Printf("Low Parent Depth: %d\n", lowParentDepth)
	fmt.Printf("One Parent Depth: %d\n", oneParentDepth)
	fmt.Printf("High Bulk Depth: %d\n", highBulkDepth)
	fmt.Printf("Low Bulk Depth: %d\n", lowBulkDepth)
	fmt.Printf("One Bulk Depth: %d\n", oneBulkDepth)
	fmt.Printf("High Bulk Size: %d\n", highBulkSize)
	fmt.Printf("Low Bulk Size: %d\n", lowBulkSize)
	fmt.Printf("One Bulk Size: %d\n", oneBulkSize)
	fmt.Printf("Window Size: %d\n", windowSize)
	fmt.Printf("Population: %s\n", population)
	fmt.Printf("Recurrent: %v\n", recurrent)
	fmt.Printf("Simulations: %d\n", rep)
	fmt.Printf("Alpha: %v\n", alpha)
	fmt.Printf("Min QTL Length: %d\n", minQTL)
	fmt.Printf("Merge Distance: %d\n", mergeDist)
	fmt.Printf("Output Dir/Prefix: %s\n", outputDir)

	color.Cyan("\n========================================== SAMPLE SELECTION =================================================\n\n")

	f, cleanup, err := openVCF(vcf)
	if err != nil {
		panic(err)
	}
	defer cleanup()

	rdr, err := vcfgo.NewReader(f, false)
	if err != nil {
		panic(err)
	}

	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	sampleNames := rdr.Header.SampleNames
	sampleNamesDic := make(map[int]string)
	sampleNamesDic[0] = "None"
	fmt.Printf("0 : None\n")
	for i, _ := range sampleNames {
		sampleNamesDic[i+1] = sampleNames[i]
		fmt.Printf("%v : %v\n", i+1, sampleNames[i])
	}
	fmt.Printf("\n=================================================================================================\n\n")
	fmt.Printf("Enter number corresponding to the sample ...\n\n")
	var highParentChoice int
	var lowParentChoice int
	var highBulkChoice int
	var lowBulkChoice int

	fmt.Println("Enter HIGH PARENT number:")
	_, highParErr := fmt.Scan(&highParentChoice)
	if highParErr != nil {
		fmt.Printf("HIGH PARENT number should be numerical and part of the list above: %s\n", highParErr)
		return highParErr
	}

	highParent := sampleNamesDic[highParentChoice]
	fmt.Printf("HIGH Parent is: %s \n\n", highParent)

	fmt.Println("Enter LOW PARENT number:")
	_, lowParErr := fmt.Scan(&lowParentChoice)
	if lowParErr != nil {
		fmt.Printf("LOW PARENT number should be numerical and part of the list above: %s\n", lowParErr)
		return lowParErr
	}

	if lowParentChoice == highParentChoice && lowParentChoice != 0 {
		fmt.Println("LOW PARENT should not be the same as HIGH PARENT")
		return fmt.Errorf("Invalid input")
	}

	lowParent := sampleNamesDic[lowParentChoice]
	fmt.Printf("LOW parent is: %s \n\n", lowParent)

	fmt.Printf("------------------------------------- BULK CHOICES ----------------------------------------\n\n")
	fmt.Println("Enter HIGH BULK number:")
	_, highBulkErr := fmt.Scan(&highBulkChoice)
	if highBulkErr != nil {
		fmt.Printf("HIGH BULK number should be numerical and part of the list above: %s\n", highBulkErr)
		return fmt.Errorf("invalid input")
	}

	if highBulkChoice == highParentChoice || highBulkChoice == lowParentChoice {
		fmt.Println("Your HIGH bulk cannot be the same as any of the parents")
		return fmt.Errorf("invalid input")
	}

	highBulk := sampleNamesDic[highBulkChoice]
	fmt.Printf("HIGH bulk is: %s \n\n", highBulk)

	fmt.Println("Enter LOW BULK number:")
	_, lowBulkErr := fmt.Scan(&lowBulkChoice)
	if lowBulkErr != nil {
		fmt.Printf("LOW BULK number should be numerical and part of the list above: %s\n", lowBulkErr)
		return fmt.Errorf("invalid input")
	}

	if lowBulkChoice == highBulkChoice || lowBulkChoice == highParentChoice || lowBulkChoice == lowParentChoice {
		fmt.Println("Your LOW bulk cannot be the same as any of the parents OR the HIGH bulk")
		return fmt.Errorf("invalid input")
	}
	lowBulk := sampleNamesDic[lowBulkChoice]
	fmt.Printf("LOW bulk is: %s \n\n", lowBulk)

	// Ensure output directory exists
	if outputDir != "" {
		dir := filepath.Dir(outputDir)
		if dir != "." {
			if err := os.MkdirAll(dir, 0755); err != nil {
				return fmt.Errorf("failed to create output directory: %w", err)
			}
		}
	}

	config := twobulk.AnalysisConfig{
		Population:       twobulk.PopulationType(strings.ToUpper(population)),
		BCAltIsRecurrent: recurrent,
		WindowSize:       windowSize,
		NSimulations:     rep,
		Alpha:            alpha,
		MinQTLWidth:      minQTL,
		MergeDistance:    mergeDist,
		OutputFile:       outputDir,
	}

	if lowBulkChoice != 0 && highBulkChoice != 0 && highParentChoice == 0 && lowParentChoice == 0 {
		fmt.Println("Running bulks only")
	} else if lowBulkChoice == 0 && highBulkChoice != 0 && lowParentChoice != 0 && highParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (HIGH bulk)...")
	} else if highBulkChoice == 0 && highParentChoice != 0 && lowParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (LOW bulk)")
	} else {
		fmt.Println("Working with two bulks")
		twobulk.RunTwoBulkTwoParentsWithConfig(rdr, highParentChoice-1, highParentDepth, lowParentChoice-1, lowParentDepth, highBulkChoice-1, highBulkDepth, lowBulkChoice-1, lowBulkDepth, config)
		//if err != nil {
		//	return err
		//}

	}
	return nil
}
