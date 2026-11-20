package stats

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
)

func site2b(pos int64, hRef, hAlt, lRef, lAlt, hOther int32) filter.Site {
	s := filter.Site{Chrom: "chr1", Pos: pos, Ref: "A", Alt: "T", NSamples: 2}
	s.Samples[0] = filter.SampleAD{Ref: hRef, Alt: hAlt, Other: hOther, HasAD: true, GTRef: 1, GTAlt: 1}
	s.Samples[1] = filter.SampleAD{Ref: lRef, Alt: lAlt, HasAD: true, GTRef: 1, GTAlt: 1}
	return s
}

// Stored values are full precision: a small ΔSI gives a non-zero ED⁴, and the
// TSV prints it rather than 0.000000.
func TestRawStatsNotRounded(t *testing.T) {
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), RefParent: "high"}
	raw, err := RawStats(cfg, "2b", []int{0, 1}, []filter.Site{site2b(100, 101, 99, 99, 101, 3)})
	if err != nil {
		t.Fatal(err)
	}
	s := raw[0]
	if s.HighSI != 0.505 || s.ED <= 0 || s.ED != 4*s.DeltaSI*s.DeltaSI*s.DeltaSI*s.DeltaSI {
		t.Fatalf("HighSI=%v DeltaSI=%v ED=%v: want full-precision values, ED = 4·ΔSI⁴ > 0", s.HighSI, s.DeltaSI, s.ED)
	}
	tsv, err := os.ReadFile(filepath.Join(cfg.OutputDir, "stats", "GoBSAseq.2b.raw.tsv"))
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSpace(string(tsv)), "\n")
	hdr, row := strings.Split(lines[0], "\t"), strings.Split(lines[1], "\t")
	col := map[string]string{}
	for i, h := range hdr {
		col[h] = row[i]
	}
	if col["ED4"] != "4e-08" || col["HighBulkAD"] != "101,99" || col["HighBulkOther"] != "3" || col["LowBulkOther"] != "0" {
		t.Fatalf("TSV row %v", col)
	}
}

func TestRawStatsSampleCountMismatchIsError(t *testing.T) {
	s := site2b(100, 10, 10, 10, 10, 0)
	s.NSamples = 3
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir()}
	if _, err := RawStats(cfg, "2b", []int{0, 1}, []filter.Site{s}); err == nil {
		t.Fatal("site with the wrong sample count was not an error")
	}
}

// A failed write (here: a full device) must be reported, not leave a
// truncated TSV behind a nil error.
func TestWriteRawTSVReportsWriteErrors(t *testing.T) {
	if _, err := os.Stat("/dev/full"); err != nil {
		t.Skip("/dev/full not available")
	}
	sites := benchSites(2, 50000) // more output than the 1 MB bufio buffer
	calc := rawCalc{roles: siteRoles{hp: -1, lp: -1, hb: 0, lb: 1}, signedDSI: true}
	if _, err := computeAndWriteRaw("/dev/full", sites, calc, rawLayout{highBulk: true, lowBulk: true, signedDSI: true}); err == nil {
		t.Fatal("writing to a full device returned nil")
	}
}

// Known 2×2 tables: standard G (no pseudo-counts, 0·ln 0 = 0).
func TestGStatisticKnownValues(t *testing.T) {
	cases := []struct {
		a1, b1, a2, b2 int
		want           float64
	}{
		{30, 10, 10, 30, 20.9299}, // 2·[2·30·ln(30/20) + 2·10·ln(10/20)]
		{20, 20, 20, 20, 0},       // identical bulks
		{40, 0, 0, 40, 110.9035},  // fully separated: 2·(2·40·ln 2)
		{0, 0, 10, 10, 0},         // empty bulk: no information
		{10, 0, 10, 0, 0},         // monomorphic: no information
	}
	for _, c := range cases {
		g := GStatistic(c.a1, c.b1, c.a2, c.b2)
		if math.Abs(g-c.want) > 1e-3 {
			t.Errorf("G%v = %.4f, want %.4f", []int{c.a1, c.b1, c.a2, c.b2}, g, c.want)
		}
		if g2 := GStatistic(c.b1, c.a1, c.b2, c.a2); math.Abs(g-g2) > 1e-12 {
			t.Errorf("G not symmetric in allele order: %v vs %v", g, g2)
		}
	}
}

func TestED4(t *testing.T) {
	if got := euclideanDistance4(0.75, 0.25); got != 0.25 { // 4·0.5⁴
		t.Fatalf("ED4(0.75, 0.25) = %v, want 0.25", got)
	}
}

// Unequal bulk depths are kept separately; Depth is their harmonic mean.
func TestRawStatsBulkDepths(t *testing.T) {
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), RefParent: "high"}
	raw, err := RawStats(cfg, "2b", []int{0, 1}, []filter.Site{site2b(100, 40, 40, 12, 13, 0)})
	if err != nil {
		t.Fatal(err)
	}
	s := raw[0]
	if s.HighDepth != 80 || s.LowDepth != 25 || s.Depth != 38 { // 2·80·25/105 = 38.1
		t.Fatalf("HighDepth/LowDepth/Depth = %d/%d/%d, want 80/25/38", s.HighDepth, s.LowDepth, s.Depth)
	}
	if h, l := siteBulkDepths(SmoothedStats{Depth: 38, HighDepth: 80, LowDepth: 25}); h != 80 || l != 25 {
		t.Fatalf("threshold depth pair = %d/%d, want 80/25", h, l)
	}
}

// sitesHPHB builds one-parent one-bulk sites whose bulk shows a REF mapping
// bias: 55 REF / 45 ALT reads at every segregating site, half with the high
// allele on REF and half on ALT, plus nonSeg sites where the bulk is fixed.
func sitesHPHB(n, nonSeg int) []filter.Site {
	var out []filter.Site
	for i := 0; i < n+nonSeg; i++ {
		s := filter.Site{Chrom: "chr1", Pos: int64(i + 1), Ref: "A", Alt: "T", NSamples: 2, HighIsAlt: i%2 == 1}
		s.Samples[0] = filter.SampleAD{Ref: 30, HasAD: true, GTRef: 2}
		if s.HighIsAlt {
			s.Samples[0] = filter.SampleAD{Alt: 30, HasAD: true, GTAlt: 2}
		}
		s.Samples[1] = filter.SampleAD{Ref: 55, Alt: 45, HasAD: true, GTRef: 1, GTAlt: 1}
		if i >= n { // bulk fixed for the parent's allele: not segregating
			s.Samples[1] = filter.SampleAD{Ref: 100, HasAD: true, GTRef: 2}
			s.HighIsAlt = false
		}
		out = append(out, s)
	}
	return out
}

func TestOneBulkP0Empirical(t *testing.T) {
	sites := sitesHPHB(1000, 400)
	for _, tc := range []struct {
		mode       string
		p0Ref      float64
		p0Alt      float64
		afDevRef   float64
		afDevAltHi float64
	}{
		{"empirical", 0.55, 0.45, 0, 0},
		{"", 0.55, 0.45, 0, 0}, // unset = empirical
		{"theoretical", 0.5, 0.5, 0.05, -0.05},
	} {
		cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), Population: "F2", OneBulkP0Mode: tc.mode}
		raw, err := RawStats(cfg, "hphb", []int{0, 1}, sites)
		if err != nil {
			t.Fatal(err)
		}
		ref, alt := raw[0], raw[1] // site 0: high allele REF; site 1: ALT
		if math.Abs(ref.OneBulkP0-tc.p0Ref) > 1e-9 || math.Abs(alt.OneBulkP0-tc.p0Alt) > 1e-9 {
			t.Fatalf("mode %q: p0 = %v / %v, want %v / %v", tc.mode, ref.OneBulkP0, alt.OneBulkP0, tc.p0Ref, tc.p0Alt)
		}
		if math.Abs(ref.OneBulkAFDev-tc.afDevRef) > 1e-9 || math.Abs(alt.OneBulkAFDev-tc.afDevAltHi) > 1e-9 {
			t.Fatalf("mode %q: AFDev = %v / %v, want %v / %v", tc.mode, ref.OneBulkAFDev, alt.OneBulkAFDev, tc.afDevRef, tc.afDevAltHi)
		}
		if tc.mode != "theoretical" && ref.OneBulkGstat > 1e-9 {
			t.Fatalf("mode %q: G1 = %v at a site matching the empirical p0, want 0", tc.mode, ref.OneBulkGstat)
		}
	}
}

func TestOneBulkP0FallsBackWithFewSites(t *testing.T) {
	p0 := oneBulkP0ByPolarity("empirical", sitesHPHB(100, 0), 1, 0.5) // 50 sites per class
	if p0 != [2]float64{0.5, 0.5} {
		t.Fatalf("p0 = %v with too few sites, want theoretical 0.5", p0)
	}
}

// A shift common to both polarity classes (a QTL, or a wrong population
// setting) must not leak into p0; only the opposite-direction mapping bias does.
func TestOneBulkP0IgnoresCommonShift(t *testing.T) {
	sites := sitesHPHB(1000, 0)
	for i := range sites {
		if sites[i].HighIsAlt { // high-allele fraction 0.50 (ALT 50 / REF 50)
			sites[i].Samples[1] = filter.SampleAD{Ref: 50, Alt: 50, HasAD: true}
		} else { // high-allele fraction 0.60 (REF 60 / ALT 40)
			sites[i].Samples[1] = filter.SampleAD{Ref: 60, Alt: 40, HasAD: true}
		}
	}
	// Classes at 0.60 and 0.50: common shift +0.05, bias +0.05.
	p0 := oneBulkP0ByPolarity("empirical", sites, 1, 0.5)
	if math.Abs(p0[0]-0.55) > 1e-9 || math.Abs(p0[1]-0.45) > 1e-9 {
		t.Fatalf("p0 = %v, want [0.55 0.45]", p0)
	}
}
