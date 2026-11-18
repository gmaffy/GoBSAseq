package filter_test

import (
	"fmt"
	"strings"
	"testing"

	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
)

// GATK hard-filtering recommendations ("(How to) Filter variants either with
// VQSR or by hard-filtering"; "Hard-filtering germline short variants"):
//
//	SNPs:   QD < 2.0, QUAL < 30.0, SOR > 3.0, FS > 60.0, MQ < 40.0,
//	        MQRankSum < -12.5, ReadPosRankSum < -8.0   → filtered
//	Indels: QD < 2.0, QUAL < 30.0, FS > 200.0, ReadPosRankSum < -20.0
//
// A record fails only when an expression is true, so a value exactly at a
// cutoff passes, and a record lacking an annotation passes that filter.
// GoBSAseq's defaults (cmd/root.go) must reproduce this in strict mode.
func gatkDefaults() utils.HardFilterConfig {
	return utils.HardFilterConfig{
		SNP_QD_Min: 2.0, SNP_QUAL_Min: 30.0, SNP_SOR_Max: 3.0, SNP_FS_Max: 60.0,
		SNP_MQ_Min: 40.0, SNP_MQRankSum_Min: -12.5, SNP_ReadPosRankSum_Min: -8.0,
		INDEL_QD_Min: 2.0, INDEL_QUAL_Min: 30.0, INDEL_FS_Max: 200.0,
		INDEL_ReadPosRankSum_Min: -20.0, INDEL_SOR_Max: 10.0,
	}
}

const gatkHeader = "##fileformat=VCFv4.2\n" +
	"##INFO=<ID=QD,Number=1,Type=Float,Description=\"QD\">\n" +
	"##INFO=<ID=FS,Number=1,Type=Float,Description=\"FS\">\n" +
	"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"SOR\">\n" +
	"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"MQ\">\n" +
	"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"MQRankSum\">\n" +
	"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"ReadPosRankSum\">\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"

// passes builds one record and runs the strict GATK-default hard filter.
func passes(t *testing.T, ref, alt string, qual float64, info string) bool {
	t.Helper()
	rec := fmt.Sprintf("chr1\t100\t.\t%s\t%s\t%g\tPASS\t%s\tGT\t0/1\n", ref, alt, qual, info)
	r, err := vcfgo.NewReader(strings.NewReader(gatkHeader+rec), false)
	if err != nil {
		t.Fatal(err)
	}
	v := r.Read()
	if v == nil {
		t.Fatalf("could not parse %q", rec)
	}
	return filter.PassesHardFilter(v, gatkDefaults())
}

func TestHardFilterMatchesGATKRecommendations(t *testing.T) {
	const good = "QD=20;FS=1;SOR=0.7;MQ=60;MQRankSum=0;ReadPosRankSum=0"
	set := func(key, val string) string { // replace one annotation in good
		var out []string
		for _, kv := range strings.Split(good, ";") {
			if strings.HasPrefix(kv, key+"=") {
				kv = key + "=" + val
			}
			out = append(out, kv)
		}
		return strings.Join(out, ";")
	}
	drop := func(key string) string {
		var out []string
		for _, kv := range strings.Split(good, ";") {
			if !strings.HasPrefix(kv, key+"=") {
				out = append(out, kv)
			}
		}
		return strings.Join(out, ";")
	}

	type check struct {
		name       string
		ref, alt   string
		qual       float64
		info       string
		wantPasses bool
	}
	snp := func(name string, qual float64, info string, want bool) check { return check{name, "A", "T", qual, info, want} }
	indel := func(name string, qual float64, info string, want bool) check {
		return check{name, "A", "AT", qual, info, want}
	}
	checks := []check{
		snp("SNP clean", 50, good, true),
		snp("SNP QUAL 30 (at cutoff)", 30, good, true),
		snp("SNP QUAL 29.9", 29.9, good, false),
		snp("SNP QD 2.0", 50, set("QD", "2.0"), true),
		snp("SNP QD 1.99", 50, set("QD", "1.99"), false),
		snp("SNP SOR 3.0", 50, set("SOR", "3.0"), true),
		snp("SNP SOR 3.01", 50, set("SOR", "3.01"), false),
		snp("SNP FS 60", 50, set("FS", "60.0"), true),
		snp("SNP FS 60.1", 50, set("FS", "60.1"), false),
		snp("SNP MQ 40", 50, set("MQ", "40.0"), true),
		snp("SNP MQ 39.9", 50, set("MQ", "39.9"), false),
		snp("SNP MQRankSum -12.5", 50, set("MQRankSum", "-12.5"), true),
		snp("SNP MQRankSum -12.6", 50, set("MQRankSum", "-12.6"), false),
		snp("SNP ReadPosRankSum -8", 50, set("ReadPosRankSum", "-8.0"), true),
		snp("SNP ReadPosRankSum -8.1", 50, set("ReadPosRankSum", "-8.1"), false),
		// Rank sums exist only at het sites; a missing annotation must pass.
		snp("SNP no MQRankSum/ReadPosRankSum", 50, drop("MQRankSum")+";ReadPosRankSum=0", true),
		snp("SNP no annotations", 50, ".", true),

		indel("indel clean", 50, good, true),
		indel("indel QUAL 29.9", 29.9, good, false),
		indel("indel QD 1.99", 50, set("QD", "1.99"), false),
		indel("indel FS 200", 50, set("FS", "200.0"), true),
		indel("indel FS 200.1", 50, set("FS", "200.1"), false),
		indel("indel ReadPosRankSum -20", 50, set("ReadPosRankSum", "-20.0"), true),
		indel("indel ReadPosRankSum -20.1", 50, set("ReadPosRankSum", "-20.1"), false),
		// GATK applies no MQ / MQRankSum filter to indels.
		indel("indel MQ 10 (not filtered for indels)", 50, set("MQ", "10"), true),
		indel("indel MQRankSum -30 (not filtered for indels)", 50, set("MQRankSum", "-30"), true),
		indel("indel SNP-level FS 100 (allowed for indels)", 50, set("FS", "100"), true),
	}
	for _, c := range checks {
		if got := passes(t, c.ref, c.alt, c.qual, c.info); got != c.wantPasses {
			t.Errorf("%s: passes=%v, want %v", c.name, got, c.wantPasses)
		}
	}
}
