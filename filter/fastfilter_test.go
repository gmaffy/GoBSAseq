package filter_test

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/tabix"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
)

const hdr2 = "##fileformat=VCFv4.2\n" +
	"##contig=<ID=chr1>\n##contig=<ID=chr2>\n" +
	"##INFO=<ID=QD,Number=1,Type=Float,Description=\"QualByDepth\">\n" +
	"##INFO=<ID=FS,Number=1,Type=Float,Description=\"FisherStrand\">\n" +
	"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"StrandOddsRatio\">\n" +
	"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMSMappingQuality\">\n" +
	"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"MQRankSum\">\n" +
	"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"ReadPosRankSum\">\n" +
	"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count, per ALT\">\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n" +
	"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n" +
	"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
	"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">\n"

func header(samples ...string) string {
	return hdr2 + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + strings.Join(samples, "\t") + "\n"
}

// gatk returns the GATK default hard-filter thresholds (strict mode).
func gatk() utils.HardFilterConfig {
	return utils.HardFilterConfig{
		SNP_QD_Min: 2, SNP_QUAL_Min: 30, SNP_SOR_Max: 3, SNP_FS_Max: 60, SNP_MQ_Min: 40,
		SNP_MQRankSum_Min: -12.5, SNP_ReadPosRankSum_Min: -8,
		INDEL_QD_Min: 2, INDEL_QUAL_Min: 30, INDEL_FS_Max: 200, INDEL_SOR_Max: 10, INDEL_ReadPosRankSum_Min: -20,
	}
}

func bulksCfg() utils.AnalysisConfig {
	return utils.AnalysisConfig{HighBulkIdx: 0, LowBulkIdx: 1, HighBulkDepth: 1, LowBulkDepth: 1, SplitMultiallelic: true}
}

func parentsCfg() utils.AnalysisConfig {
	return utils.AnalysisConfig{
		HighParentIdx: 0, LowParentIdx: 1, HighBulkIdx: 2, LowBulkIdx: 3,
		HighParentDepth: 1, LowParentDepth: 1, HighBulkDepth: 1, LowBulkDepth: 1, SplitMultiallelic: true,
	}
}

type runOut struct {
	sites    []filter.Site
	entering int
	pass     []string // data lines of the hard-filtered VCF
	rej      []string // data lines of the rejected VCF
	dir      string
}

func readGz(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gz, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip %s: %v", path, err)
	}
	raw, err := io.ReadAll(gz)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	var out []string
	for _, l := range strings.Split(strings.TrimRight(string(raw), "\n"), "\n") {
		if l != "" && !strings.HasPrefix(l, "#") {
			out = append(out, l)
		}
	}
	return out
}

func run(t *testing.T, vcf string, cfg utils.AnalysisConfig, hf utils.HardFilterConfig, mode string, keep []int) runOut {
	t.Helper()
	dir := t.TempDir()
	cfg.VCF = filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(cfg.VCF, []byte(vcf), 0o644); err != nil {
		t.Fatal(err)
	}
	cfg.OutputDir = dir
	sites, entering, nPass, err := filter.HardFilterVcf(cfg, hf, mode, keep)
	if err != nil {
		t.Fatalf("HardFilterVcf: %v", err)
	}
	if nPass != len(sites) {
		t.Fatalf("pass count %d != len(sites) %d", nPass, len(sites))
	}
	return runOut{
		sites: sites, entering: entering, dir: dir,
		pass: readGz(t, filepath.Join(dir, "stats", "GoBSAseq."+mode+".hardfiltered.vcf.gz")),
		rej:  readGz(t, filepath.Join(dir, "stats", "GoBSAseq."+mode+".lowqual.vcf.gz")),
	}
}

func rec(chrom string, pos int, ref, alt, qual, filt, info, format string, samples ...string) string {
	return fmt.Sprintf("%s\t%d\t.\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", chrom, pos, ref, alt, qual, filt, info, format, strings.Join(samples, "\t"))
}

// ── ported legacy behaviours ──────────────────────────────────────────────────

func TestSelectedAltFlowsIntoRawStats(t *testing.T) {
	vcf := header("high", "low") + rec("chr1", 100, "A", "C,G", "60", "PASS", ".", "GT:AD:DP", "0/2:10,0,30:40", "0/2:10,0,30:40")
	for _, split := range []bool{true, false} {
		cfg := bulksCfg()
		cfg.SplitMultiallelic = split
		out := run(t, vcf, cfg, gatk(), "2b", []int{0, 1})
		if len(out.sites) != 1 || out.sites[0].Alt != "G" {
			t.Fatalf("split=%v: want one site with ALT G, got %+v", split, out.sites)
		}
		cfg.OutputDir = t.TempDir()
		raw, err := stats.RawStats(cfg, "2b", []int{0, 1}, out.sites)
		if err != nil {
			t.Fatalf("RawStats: %v", err)
		}
		if len(raw) != 1 || raw[0].ALT != "G" || raw[0].HighBulkAD != "10,30" {
			t.Fatalf("split=%v: raw stats used wrong ALT: %+v", split, raw)
		}
	}
}

func TestOutputPreservesInputOrderAcrossBatches(t *testing.T) {
	var b strings.Builder
	b.WriteString(header("high", "low"))
	const n = 30000 // several worker batches
	for i := 1; i <= n; i++ {
		b.WriteString(rec("chr1", i*10, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/1:10,10:20", "0/1:10,10:20"))
	}
	out := run(t, b.String(), bulksCfg(), gatk(), "2b", []int{0, 1})
	if len(out.sites) != n || len(out.pass) != n {
		t.Fatalf("got %d sites / %d lines, want %d", len(out.sites), len(out.pass), n)
	}
	for i, s := range out.sites {
		if s.Pos != int64((i+1)*10) {
			t.Fatalf("site %d at %d, want %d", i, s.Pos, (i+1)*10)
		}
		if !strings.HasPrefix(out.pass[i], fmt.Sprintf("chr1\t%d\t", (i+1)*10)) {
			t.Fatalf("output line %d out of order: %.30q", i, out.pass[i])
		}
	}
}

func TestCallerFailedRecordRejected(t *testing.T) {
	vcf := header("high", "low") + rec("chr1", 100, "A", "T", "60", "LowQual", ".", "GT:AD:DP", "0/1:10,10:20", "0/1:10,10:20")
	out := run(t, vcf, bulksCfg(), gatk(), "2b", []int{0, 1})
	if len(out.sites) != 0 || len(out.rej) != 1 {
		t.Fatalf("FILTER=LowQual record passed: sites=%d rej=%d", len(out.sites), len(out.rej))
	}
}

func TestSecondAltTolerance(t *testing.T) {
	// Target G 29 reads, stray C 1 read (2.5%) → kept as REF/G.
	keep := header("high", "low") + rec("chr1", 100, "A", "C,G", "60", "PASS", ".", "GT:AD:DP", "0/2:10,1,29:40", "0/2:10,1,29:40")
	out := run(t, keep, bulksCfg(), gatk(), "2b", []int{0, 1})
	if len(out.sites) != 1 || out.sites[0].Alt != "G" {
		t.Fatalf("clean site with a stray third-allele read dropped: %+v", out.sites)
	}
	// Balanced three-allele pileup → neither biallelic view is representative.
	drop := header("high", "low") + rec("chr1", 100, "A", "C,G", "60", "PASS", ".", "GT:AD:DP", "0/2:10,10,20:40", "0/2:10,10,20:40")
	if out := run(t, drop, bulksCfg(), gatk(), "2b", []int{0, 1}); len(out.sites) != 0 {
		t.Fatalf("multi-allelic pileup kept: %+v", out.sites)
	}
}

func TestMaxDepthCap(t *testing.T) {
	vcf := header("high", "low") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/1:500,500:1000", "0/1:500,500:1000")
	if out := run(t, vcf, bulksCfg(), gatk(), "2b", []int{0, 1}); len(out.sites) != 1 {
		t.Fatal("high-depth site rejected with no cap")
	}
	cfg := bulksCfg()
	cfg.HighBulkMaxDepth, cfg.LowBulkMaxDepth = 200, 200
	if out := run(t, vcf, cfg, gatk(), "2b", []int{0, 1}); len(out.sites) != 0 {
		t.Fatal("coverage outlier passed the max-depth cap")
	}
}

func TestSplitRecordReindexesADGTPLAndInfo(t *testing.T) {
	// Sample "high" is 1/2. Both ALTs carry enough reads to be analysed only if
	// the other-ALT tolerance is relaxed, so check the rendering through the G
	// record of a site where C is a stray allele.
	vcf := header("high", "low") + rec("chr1", 100, "A", "C,G", "60", "PASS", "AC=1,3;QD=20", "GT:AD:DP:PL",
		"1/2:20,1,20:41:100,50,40,30,20,10", "0/2:20,0,20:40:100,50,40,30,20,10")
	out := run(t, vcf, bulksCfg(), gatk(), "2b", []int{0, 1})
	if len(out.pass) != 1 {
		t.Fatalf("want 1 passing split record, got %d: %v", len(out.pass), out.pass)
	}
	f := strings.Split(out.pass[0], "\t")
	if f[4] != "G" {
		t.Fatalf("ALT = %q, want G", f[4])
	}
	if f[7] != "QD=20" {
		t.Fatalf("per-ALT INFO not dropped: INFO=%q", f[7])
	}
	// GT 1/2 → 0/1 for the G record (other ALT → 0); AD ref,G; PL (0,0),(0,2),(2,2) = idx 0,3,5.
	if f[9] != "0/1:20,20:41:100,30,10" {
		t.Fatalf("sample high = %q, want 0/1:20,20:41:100,30,10", f[9])
	}
}

func TestMinGQAppliesToParentsNotBulks(t *testing.T) {
	cfg := parentsCfg()
	cfg.MinGQ = 20
	lowGQBulks := header("hp", "lp", "hb", "lb") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP:GQ",
		"0/0:30,0:30:60", "1/1:0,30:30:60", "0/1:20,20:40:5", "0/1:20,20:40:5")
	if out := run(t, lowGQBulks, cfg, gatk(), "2p2b", []int{0, 1, 2, 3}); len(out.sites) != 1 {
		t.Fatal("GQ floor wrongly rejected a site for low BULK GQ")
	}
	lowGQParent := header("hp", "lp", "hb", "lb") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP:GQ",
		"0/0:30,0:30:5", "1/1:0,30:30:60", "0/1:20,20:40:60", "0/1:20,20:40:60")
	if out := run(t, lowGQParent, cfg, gatk(), "2p2b", []int{0, 1, 2, 3}); len(out.sites) != 0 {
		t.Fatal("GQ floor failed to reject a site for low PARENT GQ")
	}
	if out := run(t, lowGQParent, parentsCfg(), gatk(), "2p2b", []int{0, 1, 2, 3}); len(out.sites) != 1 {
		t.Fatal("site rejected with GQ floor disabled")
	}
}

func TestDetectFilterProfile(t *testing.T) {
	h, err := filter.ParseHeader(bufio.NewReader(strings.NewReader(header("a", "b"))))
	if err != nil {
		t.Fatal(err)
	}
	p := filter.DetectFilterProfile(h)
	if !p.HasQD || !p.HasFS || !p.HasGQ || !p.HasGATKAnnotations() {
		t.Fatalf("profile missed declared annotations: %+v", p)
	}
	dv := "##fileformat=VCFv4.2\n##DeepVariant_version=1.6\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, with, commas\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n"
	h2, err := filter.ParseHeader(bufio.NewReader(strings.NewReader(dv)))
	if err != nil {
		t.Fatal(err)
	}
	p2 := filter.DetectFilterProfile(h2)
	if p2.HasGATKAnnotations() || !p2.HasGQ || p2.Caller != "deepvariant" {
		t.Fatalf("DeepVariant header misdetected: %+v", p2)
	}
}

// ── GATK hard-filter semantics ────────────────────────────────────────────────

func TestGATKThresholdBoundaries(t *testing.T) {
	cases := []struct {
		name, ref, alt, info string
		light, want          bool
	}{
		{"QD at threshold passes", "A", "T", "QD=2.0", false, true},
		{"QD below threshold fails", "A", "T", "QD=1.99", false, false},
		{"FS at threshold passes", "A", "T", "FS=60.0", false, true},
		{"FS above threshold fails", "A", "T", "FS=60.01", false, false},
		{"SOR above threshold fails", "A", "T", "SOR=3.1", false, false},
		{"MQ below threshold fails", "A", "T", "MQ=39.9", false, false},
		{"MQRankSum fails in strict mode", "A", "T", "MQRankSum=-13", false, false},
		{"MQRankSum ignored in light mode", "A", "T", "MQRankSum=-13", true, true},
		{"FS kept in light mode", "A", "T", "FS=61", true, false},
		{"missing annotations pass", "A", "T", ".", false, true},
		{"'.' and NaN count as missing", "A", "T", "QD=.;FS=nan", false, true},
		{"indel uses indel FS threshold", "A", "AT", "FS=100", false, true},
		{"indel fails above indel FS", "A", "AT", "FS=201", false, false},
		{"indel ignores MQ", "A", "AT", "MQ=10", false, true},
		{"indel ReadPosRankSum threshold", "A", "AT", "ReadPosRankSum=-21", false, false},
	}
	for _, c := range cases {
		vcf := header("high", "low") + rec("chr1", 100, c.ref, c.alt, "60", ".", c.info, "GT:AD:DP", "0/1:10,10:20", "0/1:10,10:20")
		hf := gatk()
		hf.LightFilter = c.light
		if got := len(run(t, vcf, bulksCfg(), hf, "2b", []int{0, 1}).sites) == 1; got != c.want {
			t.Errorf("%s: passed=%v, want %v", c.name, got, c.want)
		}
	}
	// QUAL below 30 fails.
	vcf := header("high", "low") + rec("chr1", 100, "A", "T", "29.9", "PASS", ".", "GT:AD:DP", "0/1:10,10:20", "0/1:10,10:20")
	if len(run(t, vcf, bulksCfg(), gatk(), "2b", []int{0, 1}).sites) != 0 {
		t.Error("QUAL 29.9 passed")
	}
}

// ── BSA-seq fixes ─────────────────────────────────────────────────────────────

// A parent homozygous for a different ALT must not be read as a confident REF
// parent after splitting (the old split remapped its GT to 0/0).
func TestParentCarryingOtherAltIsNotREF(t *testing.T) {
	vcf := header("hp", "lp", "hb", "lb") + rec("chr1", 100, "A", "C,G", "60", "PASS", ".", "GT:AD:DP",
		"2/2:0,0,30:30", "1/1:0,30,0:30", "1/2:0,20,20:40", "1/2:0,20,20:40")
	out := run(t, vcf, parentsCfg(), gatk(), "2p2b", []int{0, 1, 2, 3})
	if len(out.sites) != 0 {
		t.Fatalf("site kept although neither parent carries REF: %+v", out.sites)
	}
}

// Without allele depths a parent's GT is used, but only a REF/target-ALT
// homozygote counts; a GT naming another ALT is not confident.
func TestParentGTFallbackRejectsOtherAlt(t *testing.T) {
	vcf := header("hp", "lp", "hb", "lb") + rec("chr1", 100, "A", "C,G", "60", "PASS", ".", "GT:AD:DP",
		"2/2:.:30", "1/1:.:30", "0/1:20,20,0:40", "0/1:20,20,0:40")
	out := run(t, vcf, parentsCfg(), gatk(), "2p2b", []int{0, 1, 2, 3})
	if len(out.sites) != 0 {
		t.Fatalf("GT 2/2 parent treated as REF: %+v", out.sites)
	}
}

// Bulk minimum depth uses REF+ALT allele depth (what the SNP-index uses), not
// FORMAT/DP, which also counts uninformative reads.
func TestBulkDepthUsesAlleleDepth(t *testing.T) {
	vcf := header("high", "low") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/1:10,10:60", "0/1:10,10:60")
	cfg := bulksCfg()
	cfg.HighBulkDepth, cfg.LowBulkDepth = 30, 30
	if out := run(t, vcf, cfg, gatk(), "2b", []int{0, 1}); len(out.sites) != 0 {
		t.Fatal("bulk with 20 informative reads passed a 30-read minimum because DP=60")
	}
}

func TestPolarityFromParents(t *testing.T) {
	// High parent carries ALT → HighIsAlt; HighSI must be the ALT fraction.
	vcf := header("hp", "lp", "hb", "lb") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP",
		"1/1:0,30:30", "0/0:30,0:30", "0/1:10,30:40", "0/1:30,10:40")
	out := run(t, vcf, parentsCfg(), gatk(), "2p2b", []int{0, 1, 2, 3})
	if len(out.sites) != 1 || !out.sites[0].HighIsAlt {
		t.Fatalf("want HighIsAlt site, got %+v", out.sites)
	}
	cfg := parentsCfg()
	cfg.OutputDir = t.TempDir()
	raw, err := stats.RawStats(cfg, "2p2b", []int{0, 1, 2, 3}, out.sites)
	if err != nil {
		t.Fatal(err)
	}
	if raw[0].HighSI != 0.75 || raw[0].LowSI != 0.25 {
		t.Fatalf("HighSI/LowSI = %v/%v, want 0.75/0.25", raw[0].HighSI, raw[0].LowSI)
	}
}

// Without two parents, sites where both bulks are (near-)fixed for one allele
// do not segregate in the cross and are removed by the pooled-MAF floor.
func TestSegregationFloorForBulkModes(t *testing.T) {
	vcf := header("high", "low") +
		rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP", "1/1:1,39:40", "1/1:0,40:40") + // fixed ALT
		rec("chr1", 200, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/1:5,35:40", "0/1:25,15:40") // segregating
	cfg := bulksCfg()
	cfg.MinPooledMAF = -1 // auto
	cfg.Population = "F2"
	out := run(t, vcf, cfg, gatk(), "2b", []int{0, 1})
	if len(out.sites) != 1 || out.sites[0].Pos != 200 {
		t.Fatalf("auto floor: want only pos 200, got %+v", out.sites)
	}
	cfg.MinPooledMAF = 0 // disabled
	if out := run(t, vcf, cfg, gatk(), "2b", []int{0, 1}); len(out.sites) != 2 {
		t.Fatalf("floor disabled: want 2 sites, got %d", len(out.sites))
	}
	// Backcross: the expected minor AF is 0.25 (BC1), floor 0.1, so a pooled
	// ALT frequency of 0.85 (minor 0.15) is still segregating.
	bc := header("high", "low") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/1:6,34:40", "0/1:6,34:40")
	cfg.MinPooledMAF, cfg.Population = -1, "BC1H"
	if out := run(t, bc, cfg, gatk(), "2b", []int{0, 1}); len(out.sites) != 1 {
		t.Fatal("BC1 site with minor AF 0.15 removed by the segregation floor")
	}
}

func TestLoneParentMustBeHomozygous(t *testing.T) {
	het := header("hp", "hb", "lb") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP",
		"0/1:16,14:30", "0/1:20,20:40", "0/1:20,20:40")
	cfg := utils.AnalysisConfig{HighParentIdx: 0, HighBulkIdx: 1, LowBulkIdx: 2, HighParentDepth: 1, HighBulkDepth: 1, LowBulkDepth: 1}
	if out := run(t, het, cfg, gatk(), "hp2b", []int{0, 1, 2}); len(out.sites) != 0 {
		t.Fatal("heterozygous lone parent set polarity")
	}
}

func TestOneParentOneBulkNonParentAlleleIsOptIn(t *testing.T) {
	cfg := utils.AnalysisConfig{HighParentIdx: 0, HighBulkIdx: 1, HighParentDepth: 1, HighBulkDepth: 1}
	vcf := header("hp", "hb") +
		rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/0:30,0:30", "0/0:40,0:40") + // bulk fixed for parent allele
		rec("chr1", 200, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/0:30,0:30", "0/1:30,10:40") // other allele present
	if out := run(t, vcf, cfg, gatk(), "hphb", []int{0, 1}); len(out.sites) != 2 {
		t.Fatalf("default: want both sites kept, got %+v", out.sites)
	}
	cfg.RequireNonParentAllele = true
	out := run(t, vcf, cfg, gatk(), "hphb", []int{0, 1})
	if len(out.sites) != 1 || out.sites[0].Pos != 200 {
		t.Fatalf("--require-nonparent-allele: want only pos 200, got %+v", out.sites)
	}
}

// ── outputs ──────────────────────────────────────────────────────────────────

// The .tbi must locate records exactly: query regions and compare with a scan.
func TestTabixIndexIsQueryable(t *testing.T) {
	var b strings.Builder
	b.WriteString(header("high", "low"))
	const n = 40000 // many BGZF blocks per chromosome
	for _, c := range []string{"chr1", "chr2"} {
		for i := 1; i <= n; i++ {
			b.WriteString(rec(c, i*10, "A", "T", "60", "PASS", "QD=20;FS=1", "GT:AD:DP", "0/1:10,10:20", "0/1:10,10:20"))
		}
	}
	out := run(t, b.String(), bulksCfg(), gatk(), "2b", []int{0, 1})
	vcfPath := filepath.Join(out.dir, "stats", "GoBSAseq.2b.hardfiltered.vcf.gz")

	tf, err := os.Open(vcfPath + ".tbi")
	if err != nil {
		t.Fatalf("no index: %v", err)
	}
	defer tf.Close()
	tgz, err := gzip.NewReader(tf)
	if err != nil {
		t.Fatal(err)
	}
	idx, err := tabix.ReadFrom(tgz)
	if err != nil {
		t.Fatalf("read index: %v", err)
	}
	vf, err := os.Open(vcfPath)
	if err != nil {
		t.Fatal(err)
	}
	defer vf.Close()
	bz, err := bgzf.NewReader(vf, 1)
	if err != nil {
		t.Fatal(err)
	}

	for _, q := range []struct {
		chrom    string
		beg, end int // 0-based half-open
	}{{"chr1", 0, 1000}, {"chr1", 199_995, 205_000}, {"chr1", 399_000, 400_000}, {"chr2", 0, 100}, {"chr2", 123_450, 123_560}, {"chr2", 399_900, 400_000}} {
		chunks, err := idx.Chunks(q.chrom, q.beg, q.end)
		if err != nil {
			t.Fatalf("%v: %v", q, err)
		}
		cr, err := index.NewChunkReader(bz, chunks)
		if err != nil {
			t.Fatal(err)
		}
		got := 0
		sc := bufio.NewScanner(cr)
		for sc.Scan() {
			f := strings.SplitN(sc.Text(), "\t", 3)
			var pos int
			fmt.Sscan(f[1], &pos)
			if f[0] == q.chrom && pos-1 >= q.beg && pos-1 < q.end {
				got++
			}
		}
		want := 0
		for p := 10; p <= n*10; p += 10 {
			if p-1 >= q.beg && p-1 < q.end {
				want++
			}
		}
		if got != want {
			t.Errorf("%s:%d-%d: index returned %d records, want %d", q.chrom, q.beg, q.end, got, want)
		}
	}
}

func TestGzipAndPlainInputs(t *testing.T) {
	vcf := header("high", "low") + rec("chr1", 100, "A", "T", "60", "PASS", ".", "GT:AD:DP", "0/1:10,10:20", "0/1:10,10:20")
	dir := t.TempDir()
	gzPath := filepath.Join(dir, "in.vcf.gz")
	f, _ := os.Create(gzPath)
	gw := gzip.NewWriter(f) // plain gzip, not BGZF
	gw.Write([]byte(vcf))
	gw.Close()
	f.Close()
	cfg := bulksCfg()
	cfg.VCF, cfg.OutputDir = gzPath, dir
	sites, _, _, err := filter.HardFilterVcf(cfg, gatk(), "2b", []int{0, 1})
	if err != nil || len(sites) != 1 {
		t.Fatalf("plain-gzip input: sites=%d err=%v", len(sites), err)
	}
}

func TestMalformedADRejectsBulk(t *testing.T) {
	// AD has 2 values for 3 alleles → not usable → bulk fails.
	vcf := header("high", "low") + rec("chr1", 100, "A", "C,G", "60", "PASS", ".", "GT:AD:DP", "0/1:10,10:20", "0/1:10,10,0:20")
	if out := run(t, vcf, bulksCfg(), gatk(), "2b", []int{0, 1}); len(out.sites) != 0 {
		t.Fatalf("malformed AD accepted: %+v", out.sites)
	}
}
