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

func TestHardFilterRejectsCallerFailedRecord(t *testing.T) {
	_, variants := readVCF(t, "chr1\t100\t.\tA\tT\t60\tLowQual\t.\tGT:AD:DP\t0/1:10,10:20\t0/1:10,10:20\n")
	if filter.PassesHardFilter(variants[0], utils.HardFilterConfig{LightFilter: true}) {
		t.Fatal("explicitly failed VCF record passed the hard filter")
	}
}

// A stray read on a third allele must not discard an otherwise-clean biallelic
// site, but substantial third-allele support should.
func TestBsaSeqSecondAltTolerance(t *testing.T) {
	cfg := utils.AnalysisConfig{HighBulkIdx: 0, LowBulkIdx: 1, HighBulkDepth: 1, LowBulkDepth: 1}

	// REF=A ALT=C,G. Target G (29 reads), stray C (1 read) = 2.5% of depth → keep.
	_, keep := readVCF(t, "chr1\t100\t.\tA\tC,G\t60\tPASS\t.\tGT:AD:DP\t0/2:10,1,29:40\t0/2:10,1,29:40\n")
	if got := filter.BsaSeqTargetAlt(keep[0], cfg, "2b"); got != 2 {
		t.Fatalf("clean biallelic site with a single stray read was dropped: targetAlt=%d, want 2", got)
	}

	// Balanced three-allele pileup (C and G each 25%) → not analysable biallelically.
	_, drop := readVCF(t, "chr1\t100\t.\tA\tC,G\t60\tPASS\t.\tGT:AD:DP\t0/2:10,10:40\t0/2:10,10,20:40\n")
	if got := filter.BsaSeqTargetAlt(drop[0], cfg, "2b"); got != 0 {
		t.Fatalf("genuinely multi-allelic site was kept: targetAlt=%d, want 0", got)
	}
}

// The max-depth cap must reject coverage outliers (collapsed repeats / CNVs).
func TestBsaSeqMaxDepthCap(t *testing.T) {
	_, variants := readVCF(t, "chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP\t0/1:500,500:1000\t0/1:500,500:1000\n")

	base := utils.AnalysisConfig{HighBulkIdx: 0, LowBulkIdx: 1, HighBulkDepth: 1, LowBulkDepth: 1}
	if got := filter.BsaSeqTargetAlt(variants[0], base, "2b"); got == 0 {
		t.Fatal("high-depth biallelic site rejected with no cap configured")
	}

	capped := base
	capped.HighBulkMaxDepth = 200
	capped.LowBulkMaxDepth = 200
	if got := filter.BsaSeqTargetAlt(variants[0], capped, "2b"); got != 0 {
		t.Fatalf("coverage outlier passed the max-depth cap: targetAlt=%d, want 0", got)
	}
}

func TestSplitMultiallelicReindexesADandGT(t *testing.T) {
	// REF=A, ALT=C,G. Sample is 1/2 with AD=4,10,20 (ref,C,G) and DP=34.
	_, variants := readVCF(t, "chr1\t100\t.\tA\tC,G\t60\tPASS\t.\tGT:AD:DP\t1/2:4,10,20:34\t0/1:4,10,20:34\n")
	recs := filter.SplitMultiallelic(variants[0])
	if len(recs) != 2 {
		t.Fatalf("expected 2 biallelic records, got %d", len(recs))
	}

	// Record 0 → ALT C (allele 1). AD should be ref,C = 4,10. Sample0 GT 1/2 → 1/0.
	if got := recs[0].Alt(); len(got) != 1 || got[0] != "C" {
		t.Fatalf("record 0 ALT = %v, want [C]", got)
	}
	if got := recs[0].Samples[0].Fields["AD"]; got != "4,10" {
		t.Fatalf("record 0 sample0 AD = %q, want 4,10", got)
	}
	if got := recs[0].Samples[0].Fields["GT"]; got != "1/0" {
		t.Fatalf("record 0 sample0 GT = %q, want 1/0", got)
	}

	// Record 1 → ALT G (allele 2). AD should be ref,G = 4,20. Sample0 GT 1/2 → 0/1.
	if got := recs[1].Alt(); len(got) != 1 || got[0] != "G" {
		t.Fatalf("record 1 ALT = %v, want [G]", got)
	}
	if got := recs[1].Samples[0].Fields["AD"]; got != "4,20" {
		t.Fatalf("record 1 sample0 AD = %q, want 4,20", got)
	}
	if got := recs[1].Samples[0].Fields["GT"]; got != "0/1" {
		t.Fatalf("record 1 sample0 GT = %q, want 0/1", got)
	}
	if got := recs[1].Samples[0].GT; len(got) != 2 || got[0] != 0 || got[1] != 1 {
		t.Fatalf("record 1 sample0 parsed GT = %v, want [0 1]", got)
	}

	// The two split records must not share sample objects (independent Fields).
	if recs[0].Samples[0] == recs[1].Samples[0] {
		t.Fatal("split records share the same sample genotype pointer")
	}
}

func TestSplitMultiallelicReindexesGenotypeLikelihoods(t *testing.T) {
	// PL is Number=G. For ALT G (allele 2) the diploid genotypes (0,0),(0,2),(2,2)
	// map to original PL indices 0, 3, 5.
	header := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ths\n"
	r, err := vcfgo.NewReader(strings.NewReader(header+
		"chr1\t100\t.\tA\tC,G\t60\tPASS\t.\tGT:PL\t1/2:100,50,40,30,20,10\n"), false)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	v := r.Read()
	if v == nil {
		t.Fatal("no variant parsed")
	}
	recs := filter.SplitMultiallelic(v)
	if len(recs) != 2 {
		t.Fatalf("expected 2 records, got %d", len(recs))
	}
	// ALT C (allele 1): genotypes (0,0),(0,1),(1,1) → PL indices 0,1,2 → 100,50,40.
	if got := recs[0].Samples[0].Fields["PL"]; got != "100,50,40" {
		t.Fatalf("record 0 PL = %q, want 100,50,40", got)
	}
	// ALT G (allele 2): PL indices 0,3,5 → 100,30,10.
	if got := recs[1].Samples[0].Fields["PL"]; got != "100,30,10" {
		t.Fatalf("record 1 PL = %q, want 100,30,10", got)
	}
}

func TestHardFilterSplitsMultiallelicEndToEnd(t *testing.T) {
	records := "chr1\t100\t.\tA\tC,G\t60\tPASS\t.\tGT:AD:DP\t0/2:20,1,40:61\t0/2:20,1,40:61\n"
	r, err := vcfgo.NewReader(strings.NewReader(vcfHeader+records), false)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	cfg := utils.AnalysisConfig{
		Rdr:               r,
		OutputDir:         t.TempDir(),
		HighBulkIdx:       0,
		LowBulkIdx:        1,
		HighBulkDepth:     1,
		LowBulkDepth:      1,
		SplitMultiallelic: true,
	}
	passed, original, nPass, err := filter.HardFilterVcf(cfg, utils.HardFilterConfig{LightFilter: true}, "2b", []int{0, 1})
	if err != nil {
		t.Fatalf("HardFilterVcf: %v", err)
	}
	// One multi-allelic input record becomes two biallelic records entering the
	// filter; the G allele (40 reads) passes, the C allele (1 read, ~1.6%) as a
	// clean second-ALT-tolerant biallelic site also passes.
	if original != 2 {
		t.Fatalf("records entering filter = %d, want 2 (decomposed)", original)
	}
	if nPass != len(passed) {
		t.Fatalf("passed count %d != slice length %d", nPass, len(passed))
	}
	for _, fv := range passed {
		if len(fv.Variant.Alt()) != 1 {
			t.Fatalf("passed record is not biallelic: ALT=%v", fv.Variant.Alt())
		}
		if fv.TargetAlt != 1 {
			t.Fatalf("decomposed record TargetAlt=%d, want 1", fv.TargetAlt)
		}
	}
}

func TestSplitMultiallelicBiallelicPassthrough(t *testing.T) {
	_, variants := readVCF(t, "chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\t0/1:10,10:20\n")
	recs := filter.SplitMultiallelic(variants[0])
	if len(recs) != 1 || recs[0] != variants[0] {
		t.Fatalf("biallelic record should pass through unchanged, got %d records", len(recs))
	}
}

// The GQ floor must apply to parents but never to bulks (a bulk GQ floor would
// be AF-correlated and bias the SNP-index).
func TestMinGQAppliesToParentsNotBulks(t *testing.T) {
	header := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thp\tlp\thb\tlb\n"
	readOne := func(rec string) *vcfgo.Variant {
		r, err := vcfgo.NewReader(strings.NewReader(header+rec), false)
		if err != nil {
			t.Fatalf("NewReader: %v", err)
		}
		v := r.Read()
		if v == nil {
			t.Fatal("no variant parsed")
		}
		return v
	}

	// hp (idx0) homozygous ref, lp (idx1) homozygous alt → informative parents.
	// Bulks carry a LOW GQ (5); parents carry a HIGH GQ (60).
	base := utils.AnalysisConfig{
		HighParentIdx: 0, LowParentIdx: 1, HighBulkIdx: 2, LowBulkIdx: 3,
		HighParentDepth: 1, LowParentDepth: 1, HighBulkDepth: 1, LowBulkDepth: 1,
	}

	// Low-GQ BULKS, high-GQ parents: with a GQ floor of 20 the site must still
	// pass, proving the floor is not applied to bulks.
	lowGQBulks := readOne("chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP:GQ\t" +
		"0/0:30,0:30:60\t1/1:0,30:30:60\t0/1:20,20:40:5\t0/1:20,20:40:5\n")
	cfg := base
	cfg.MinGQ = 20
	if filter.BsaSeqTargetAlt(lowGQBulks, cfg, "2p2b") == 0 {
		t.Fatal("GQ floor wrongly rejected a site for low BULK GQ")
	}

	// Low-GQ PARENT (hp GQ=5): the same floor must reject the site.
	lowGQParent := readOne("chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:AD:DP:GQ\t" +
		"0/0:30,0:30:5\t1/1:0,30:30:60\t0/1:20,20:40:60\t0/1:20,20:40:60\n")
	if filter.BsaSeqTargetAlt(lowGQParent, cfg, "2p2b") != 0 {
		t.Fatal("GQ floor failed to reject a site for low PARENT GQ")
	}

	// With the floor disabled (MinGQ=0), the low-parent-GQ site passes again.
	cfg0 := base
	if filter.BsaSeqTargetAlt(lowGQParent, cfg0, "2p2b") == 0 {
		t.Fatal("site rejected with GQ floor disabled")
	}
}

func TestDetectFilterProfile(t *testing.T) {
	gatkHeader := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=QD,Number=1,Type=Float,Description=\"QualByDepth\">\n" +
		"##INFO=<ID=FS,Number=1,Type=Float,Description=\"FisherStrand\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thigh\tlow\n"
	r, err := vcfgo.NewReader(strings.NewReader(gatkHeader), false)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	p := filter.DetectFilterProfile(r.Header)
	if !p.HasQD || !p.HasFS || !p.HasGQ {
		t.Fatalf("profile missed declared annotations: %+v", p)
	}
	if !p.HasGATKAnnotations() {
		t.Fatal("HasGATKAnnotations should be true when QD/FS present")
	}

	// A DeepVariant-style header (no GATK annotations, GQ present) must not
	// report GATK annotations.
	dvHeader := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thigh\tlow\n"
	r2, err := vcfgo.NewReader(strings.NewReader(dvHeader), false)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	p2 := filter.DetectFilterProfile(r2.Header)
	if p2.HasGATKAnnotations() {
		t.Fatalf("DeepVariant-style header wrongly reported GATK annotations: %+v", p2)
	}
	if !p2.HasGQ {
		t.Fatal("GQ FORMAT not detected")
	}
}
