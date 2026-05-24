package run

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/gmaffy/GoBSAseq/utils"
)

func TestRunTwoParentTwoBulkSyntheticVCF(t *testing.T) {
	dir := t.TempDir()
	vcf := filepath.Join(dir, "input.vcf")
	out := filepath.Join(dir, "out")
	data := `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HP	LP	HB	LB
chr1	100	.	A	G	60	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:AD:DP	1/1:0,20:20	0/0:20,0:20	0/1:8,32:40	0/1:28,12:40
chr1	200	.	A	G	60	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:AD:DP	1/1:0,20:20	0/0:20,0:20	0/1:10,30:40	0/1:30,10:40
chr1	300	.	A	G	60	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:AD:DP	1/1:0,20:20	0/0:20,0:20	0/1:12,28:40	0/1:26,14:40
`
	if err := os.WriteFile(vcf, []byte(data), 0o644); err != nil {
		t.Fatal(err)
	}

	cfg := utils.AnalysisConfig{
		VCF:             vcf,
		WindowSize:      200,
		StepSize:        100,
		Rep:             2,
		Alphas:          []float64{0.05, 0.01},
		MinQTLWidth:     1,
		MergeDistance:   100,
		OutputDir:       out,
		HighParentName:  "HP",
		HighParentDepth: 5,
		LowParentName:   "LP",
		LowParentDepth:  5,
		HighBulkName:    "HB",
		HighBulkDepth:   10,
		LowBulkName:     "LB",
		LowBulkDepth:    10,
	}
	hf := utils.HardFilterConfig{
		SNP_QD_Min:             2,
		SNP_QUAL_Min:           30,
		SNP_SOR_Max:            3,
		SNP_FS_Max:             60,
		SNP_MQ_Min:             40,
		SNP_MQRankSum_Min:      -12.5,
		SNP_ReadPosRankSum_Min: -8,
	}
	if err := Run(cfg, hf); err != nil {
		t.Fatal(err)
	}
	for _, name := range []string{
		"variant_statistics.tsv",
		"smoothed_statistics.tsv",
		"thresholds.tsv",
		"qtls_by_simulated_thresholds.tsv",
		"final_qtls_by_robust_z.tsv",
		filepath.Join("plots", "delta_snp_index.html"),
		filepath.Join("plots", "delta_snp_index_robust_z.html"),
	} {
		if _, err := os.Stat(filepath.Join(out, name)); err != nil {
			t.Fatalf("expected output %s: %v", name, err)
		}
	}
}
