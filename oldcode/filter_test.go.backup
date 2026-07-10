package filter_test

import (
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
)

const vcfHeader = "##fileformat=VCFv4.2\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n" +
	"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thigh\tlow\n"

func readVCF(t *testing.T, records string) (*vcfgo.Reader, []*vcfgo.Variant) {
	t.Helper()
	r, err := vcfgo.NewReader(strings.NewReader(vcfHeader+records), false)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	var variants []*vcfgo.Variant
	for v := r.Read(); v != nil; v = r.Read() {
		variants = append(variants, v)
	}
	if err := r.Error(); err != nil {
		t.Fatalf("read VCF: %v", err)
	}
	return r, variants
}

func TestSelectedAltFlowsIntoRawStats(t *testing.T) {
	_, variants := readVCF(t, "chr1\t100\t.\tA\tC,G\t60\tPASS\t.\tGT:AD:DP\t0/2:10,0,30:40\t0/2:10,0,30:40\n")
	if len(variants) != 1 {
		t.Fatalf("got %d variants, want 1", len(variants))
	}

	cfg := utils.AnalysisConfig{HighBulkIdx: 0, LowBulkIdx: 1, HighBulkDepth: 1, LowBulkDepth: 1}
	selected := filter.BsaSeqTargetAlt(variants[0], cfg, "2b")
	if selected != 2 {
		t.Fatalf("selected ALT = %d, want 2", selected)
	}

	cfg.OutputDir = t.TempDir()
	raw, err := stats.RawStats(cfg, "2b", []int{0, 1}, []filter.FilteredVariant{{Variant: variants[0], TargetAlt: selected}})
	if err != nil {
		t.Fatalf("RawStats: %v", err)
	}
	if len(raw) != 1 || raw[0].ALT != "G" || raw[0].HighBulkAD != "10,30" {
		t.Fatalf("raw stats used wrong ALT: %+v", raw)
	}
}

func TestHardFilterPreservesInputOrder(t *testing.T) {
	records := "chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\t0/1:10,10:20\n" +
		"chr1\t200\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\t0/1:10,10:20\n" +
		"chr1\t300\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\t0/1:10,10:20\n"
	r, err := vcfgo.NewReader(strings.NewReader(vcfHeader+records), false)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	cfg := utils.AnalysisConfig{
		Rdr:           r,
		OutputDir:     t.TempDir(),
		HighBulkIdx:   0,
		LowBulkIdx:    1,
		HighBulkDepth: 1,
		LowBulkDepth:  1,
	}
	hf := utils.HardFilterConfig{LightFilter: true}
	passed, _, _, err := filter.HardFilterVcf(cfg, hf, "2b", []int{0, 1})
	if err != nil {
		t.Fatalf("HardFilterVcf: %v", err)
	}
	for i, want := range []int{100, 200, 300} {
		if got := passed[i].Variant.Pos; got != uint64(want) {
			t.Fatalf("passed record %d has position %d, want %d", i, got, want)
		}
	}

	f, err := os.Open(filepath.Join(cfg.OutputDir, "stats", "GoBSAseq.2b.hardfiltered.vcf.gz"))
	if err != nil {
		t.Fatalf("open output: %v", err)
	}
	defer f.Close()
	gz, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("open gzip: %v", err)
	}
	defer gz.Close()
	raw, err := io.ReadAll(gz)
	if err != nil {
		t.Fatalf("read output: %v", err)
	}
	if strings.Index(string(raw), "\t100\t") > strings.Index(string(raw), "\t200\t") {
		t.Fatal("hard-filtered VCF is not coordinate-sorted")
	}
}
