package twobulk

import (
	"bytes"
	"fmt"
	"os/exec"
	"strings"
	"time"
)

type VCF struct {
	CHROM string
	POS   string
	REF   string
	ALT   string
}

func filterVCFTwoBulkTwoParents(vcfFile string, highPar string, highParDP int, lowPar string, lowParDP int, highBulk string, highBulkDP int, lowBulk string, lowBulkDP int, outVCF string) error {

	fmt.Printf("Filtering VCF %s\n", vcfFile)
	fmt.Printf("Replacing | with / ....")

	unphasedVCF := strings.TrimSuffix(vcfFile, ".vcf.gz") + ".unphased.vcf.gz"
	cmdStr := fmt.Sprintf("bcftools view %s \\\n| sed 's/|/\\//g' \\\n| bcftools view -Oz -o %s\nbcftools index %s\n``", vcfFile, unphasedVCF, unphasedVCF)

	cmd := exec.Command("bash", "-c", cmdStr)
	err := cmd.Run()
	if err != nil {
		fmt.Println("CMD error:", err)
		return err
	}

	start := time.Now()

	expr := fmt.Sprintf(`
        GT[%s]!="./." && GT[%s]!="./." && GT[%s]!="./." && GT[%s]!="./." &&

        ((GT[%s]="0/0" && GT[%s]="1/1") || (GT[%s]="1/1" && GT[%s]="0/0")) &&

        DP[%s] >= %d &&
        DP[%s] >= %d &&
        DP[%s] >= %d &&
        DP[%s] >= %d &&

        AD[%s][1]>=0  && AD[%s][2]=="." &&
        AD[%s][1]>=0  && AD[%s][2]=="." &&
        AD[%s][1]>=0  && AD[%s][2]=="." &&
        AD[%s][1]>=0  && AD[%s][2]=="."
    `,
		highPar, lowPar, highBulk, lowBulk,
		highPar, lowPar, highPar, lowPar,
		highPar, highParDP,
		lowPar, lowParDP,
		highBulk, highBulkDP,
		lowBulk, lowBulkDP,
		highPar, highPar,
		lowPar, lowPar,
		highBulk, highBulk,
		lowBulk, lowBulk,
	)

	expr = strings.Join(strings.Fields(expr), " ")

	args := []string{
		"view",
		"-m2", "-M2", "-v", "snps",
		"-i", expr,
		"-Oz",
		"-o", outVCF,
		unphasedVCF,
	}

	cmd = exec.Command("bcftools", args...)

	var stderr bytes.Buffer
	cmd.Stderr = &stderr

	if err := cmd.Run(); err != nil {
		//fmt.Fprintf(log, "[ERROR] bcftools view failed:\n%s\n", stderr.String())
		return fmt.Errorf("bcftools view error: %w", err)
	}

	// index output VCF
	idxCmd := exec.Command("bcftools", "index", "-f", outVCF)
	idxCmd.Stderr = &stderr

	if err := idxCmd.Run(); err != nil {
		//fmt.Fprintf(log, "[ERROR] bcftools index failed:\n%s\n", stderr.String())
		return fmt.Errorf("bcftools index error: %w", err)
	}

	elapsed := time.Since(start)
	//fmt.Fprintf(log, "[INFO] bcftools filtering completed in %s\n", elapsed)
	fmt.Printf("Elapsed time: %s\n", elapsed)

	return nil
}

func TwoParentsTwoBulkRun(vcfFile string, highPar string, highParDP int) error {
	fmt.Printf("Filtering VCF %s\n", vcfFile)
	return nil
}
