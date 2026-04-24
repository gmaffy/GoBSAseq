package utils

import (
	"bufio"
	"fmt"
	"os/exec"
	"sort"

	"github.com/gmaffy/GoBSAseq/twobulk"
)

func GetSampleNames(vcf string) (map[int]string, []int, error) {
	sampleParametersDic := map[int]string{0: "None"}
	sampleID := 0
	var ids []int
	ids = append(ids, 0)
	cmd := exec.Command("bcftools", "query", "-l", vcf)

	stdout, err := cmd.StdoutPipe()
	if err != nil {
		return nil, nil, err
	}

	if err := cmd.Start(); err != nil {
		return nil, nil, err
	}

	scanner := bufio.NewScanner(stdout)
	for scanner.Scan() {
		line := scanner.Text()

		sampleID++
		ids = append(ids, sampleID)
		sampleParametersDic[sampleID] = line
	}

	if err := cmd.Wait(); err != nil {
		return nil, nil, err
	}

	return sampleParametersDic, ids, nil

}

func ParseVCF(vcf string) error {
	//fmt.Printf("VCF: %s\n", vcf)
	sampleParametersDic, ids, err := GetSampleNames(vcf)

	if err != nil {
		return err
	}
	fmt.Printf("================================ SAMPLE SELECTION ==========================================\n\n")

	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	sort.Ints(ids)
	for id := range ids {
		fmt.Printf("%v : %v\n", id, sampleParametersDic[id])
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

	highParent := sampleParametersDic[highParentChoice]
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

	lowParent := sampleParametersDic[lowParentChoice]
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

	highBulk := sampleParametersDic[highBulkChoice]
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
	lowBulk := sampleParametersDic[lowBulkChoice]
	fmt.Printf("LOW bulk is: %s \n\n", lowBulk)

	if lowBulkChoice != 0 && highBulkChoice != 0 && highParentChoice == 0 && lowParentChoice == 0 {
		fmt.Println("Running bulks only")
	} else if lowBulkChoice == 0 && highBulkChoice != 0 && lowParentChoice != 0 && highParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (HIGH bulk)...")
	} else if highBulkChoice == 0 && highParentChoice != 0 && lowParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (LOW bulk)")
	} else {
		fmt.Println("Working with two bulks")
		twobulk.TwoParentsTwoBulkRun(vcf, highParent, 5, lowParent, 5, highBulk, 40, lowBulk, 40)

	}

	return nil

}

func Run(vcf string, highParentDepth int, lowParentDepth int, oneParentDepth int, highBulkDepth int, lowBulkDepth int, oneBulkDepth int, highBulkSize int, lowBulkSize int, oneBulkSize int, windowSize int) error {
	fmt.Printf("VCF: %s\n", vcf)
	sampleParametersDic, ids, err := GetSampleNames(vcf)

	if err != nil {
		return err
	}
	fmt.Printf("================================ SAMPLE SELECTION ==========================================\n\n")

	fmt.Printf("Here are the samples found in your VCF file ...\n\n")
	sort.Ints(ids)
	for id := range ids {
		fmt.Printf("%v : %v\n", id, sampleParametersDic[id])
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

	highParent := sampleParametersDic[highParentChoice]
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

	lowParent := sampleParametersDic[lowParentChoice]
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

	highBulk := sampleParametersDic[highBulkChoice]
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
	lowBulk := sampleParametersDic[lowBulkChoice]
	fmt.Printf("LOW bulk is: %s \n\n", lowBulk)

	if lowBulkChoice != 0 && highBulkChoice != 0 && highParentChoice == 0 && lowParentChoice == 0 {
		fmt.Println("Running bulks only")
	} else if lowBulkChoice == 0 && highBulkChoice != 0 && lowParentChoice != 0 && highParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (HIGH bulk)...")
	} else if highBulkChoice == 0 && highParentChoice != 0 && lowParentChoice != 0 {
		fmt.Println("Working with one bulk BSAseq (LOW bulk)")
	} else {
		fmt.Println("Working with two bulks")
		err = twobulk.TwoParentsTwoBulkRun(vcf, highParent, highParentDepth, lowParent, lowParentDepth, highBulk, highBulkDepth, lowBulk, lowBulkDepth)
		if err != nil {
			return err
		}

	}

}