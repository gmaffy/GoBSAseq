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

func filterVCFTwoBulkTwoParents(vcfFile string, highPar int, highParDP int, lowPar int, lowParDP int, highBulk int, highBulkDP int, lowBulk int, lowBulkDP int, outVCF string) error {
	expr := fmt.Sprintf(`
		GT[%d]!="./." && GT[%d]!=".|." &&
		GT[%d]!="./." && GT[%d]!=".|." &&
		GT[%d]!="./." && GT[%d]!=".|." &&
		GT[%d]!="./." && GT[%d]!=".|." &&

		(((GT[%d]="0/0" || GT[%d]="0|0") && (GT[%d]="1/1" || GT[%d]="1|1")) ||
		((GT[%d]="1/1" || GT[%d]="1|1") && (GT[%d]="0/0" || GT[%d]="0|0"))) &&

		FMT/DP[%d] >= %d &&
		FMT/DP[%d] >= %d &&
		FMT/DP[%d] >= %d &&
		FMT/DP[%d] >= %d &&

		AD[%d:1]>=0 && AD[%d:2]=="." &&
		AD[%d:1]>=0 && AD[%d:2]=="." &&
		AD[%d:1]>=0 && AD[%d:2]=="." &&
		AD[%d:1]>=0 && AD[%d:2]=="."
	`,
		highPar, highPar,
		lowPar, lowPar,
		highBulk, highBulk,
		lowBulk, lowBulk,

		highPar, highPar, lowPar, lowPar,
		highPar, highPar, lowPar, lowPar,

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

	args := []string{"view", "-m2", "-M2", "-v", "snps", "-i", expr, "-Oz", "-o", outVCF, vcfFile}

	start := time.Now()
	cmd := exec.Command("bcftools", args...)

	var stderr bytes.Buffer
	cmd.Stderr = &stderr

	if err := cmd.Run(); err != nil {
		return fmt.Errorf("bcftools view error: %w: %s", err, strings.TrimSpace(stderr.String()))
	}

	// index output VCF
	stderr.Reset()
	idxCmd := exec.Command("bcftools", "index", "-f", outVCF)
	idxCmd.Stderr = &stderr

	if err := idxCmd.Run(); err != nil {
		return fmt.Errorf("bcftools index error: %w: %s", err, strings.TrimSpace(stderr.String()))
	}

	elapsed := time.Since(start)
	//fmt.Fprintf(log, "[INFO] bcftools filtering completed in %s\n", elapsed)
	fmt.Printf("Elapsed time: %s\n", elapsed)

	return nil
}

func TwoParentsTwoBulkRun(vcfFile string, highPar int, highParDP int, lowPar int, lowParDP int, highBulk int, highBulkDP int, lowBulk int, lowBulkDP int) error {
	fmt.Printf("Filtering VCF %s\n", vcfFile)

	filteredVCF := strings.TrimSuffix(vcfFile, ".vcf.gz") + ".filtered.vcf.gz"
	err := filterVCFTwoBulkTwoParents(vcfFile, highPar, highParDP, lowPar, lowParDP, highBulk, highBulkDP, lowBulk, lowBulkDP, filteredVCF)
	if err != nil {
		fmt.Println("Error filtering VCF:", err)
		return err
	}
	fmt.Printf("Filtered VCF saved to %s\n", filteredVCF)
	return nil
}
