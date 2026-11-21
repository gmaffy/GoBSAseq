package stats

import (
	"io"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
)

// Ends are interpolated where |value| − threshold crosses 0; a run reaching a
// chromosome end stops at that marker; the peak is the largest |value|.
func TestScanBRMBlocksInterpolatesEnds(t *testing.T) {
	sm := []SmoothedStats{{POS: 100}, {POS: 200}, {POS: 300}, {POS: 400}, {POS: 500}}
	vals := map[int64]float64{100: 0, 200: 0.5, 300: 2, 400: 1.5, 500: 0}
	blocks := scanBRMBlocks("c", sm, 1, func(s *SmoothedStats) (float64, float64) { return vals[s.POS], 1 })
	// margins −1, −0.5, 1, 0.5, −1: start 200 + 0.5/1.5·100, stop 400 + 0.5/1.5·100
	if len(blocks) != 1 || blocks[0].Start != 233 || blocks[0].Stop != 433 || blocks[0].PeakPos != 300 || blocks[0].Peak != 2 || blocks[0].PeakZ != 2 {
		t.Fatalf("blocks %+v", blocks)
	}
	vals[100], vals[500] = 3, -4 // runs at both chromosome ends
	blocks = scanBRMBlocks("c", sm, 1, func(s *SmoothedStats) (float64, float64) { return vals[s.POS], 1 })
	// margins 2, −0.5, 1, 0.5, 3: [100, 100 + 2/2.5·100] and [200 + 0.5/1.5·100, 500]
	if len(blocks) != 2 || blocks[0].Start != 100 || blocks[0].Stop != 180 || blocks[1].Start != 233 ||
		blocks[1].Stop != 500 || blocks[1].Peak != -4 || blocks[1].PeakPos != 500 {
		t.Fatalf("end runs %+v", blocks)
	}
}

// Low-bulk-only one-bulk modes test SmAFDev (SmHighSI is 0 there and must not
// be used): a null-sized deviation gives no block.
func TestBRMOneBulkUsesAFDev(t *testing.T) {
	cfg := utils.AnalysisConfig{Population: "F2", LowBulkSize: 20}
	m, err := newBRMModel(cfg, "2plb")
	if err != nil {
		t.Fatal(err)
	}
	sm := make([]SmoothedStats, 50)
	for i := range sm {
		sm[i] = SmoothedStats{CHROM: "c", POS: int64(1000 * (i + 1)), Depth: 60, OneBulkP0: 0.5, SmLowSI: 0.5, ReadVar0: 1e-4}
	}
	if b := m.blocks("c", sm, 0, 1.96); len(b) != 0 {
		t.Fatalf("no deviation, got blocks %+v", b)
	}
	sm[25].SmAFDev = 0.4
	if b := m.blocks("c", sm, 0, 1.96); len(b) != 1 || b[0].PeakPos != 26000 {
		t.Fatalf("one deviating marker, got %+v", b)
	}
}

func brmInput() []SmoothedStats {
	var sm []SmoothedStats
	for _, c := range []string{"chr9", "chr1", "chr3"} { // not sorted
		for i := range 40 {
			d := 0.0
			if i >= 15 && i < 25 {
				d = 0.4
			}
			sm = append(sm, SmoothedStats{CHROM: c, POS: int64(10_000 * (i + 1)), HighDepth: 50, LowDepth: 50, Depth: 50,
				SmHighSI: 0.5 + d/2, SmLowSI: 0.5 - d/2, SmDeltaSI: d, ReadVar0: 1e-4, ReadVar1: 1e-2})
		}
	}
	return sm
}

// Blocks come out in sorted chromosome order on every run.
func TestRunBRMSortedAndStable(t *testing.T) {
	color.Output = io.Discard
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), Population: "F2", HighBulkSize: 20, LowBulkSize: 20, BrmAlpha: 0.05, RefParent: "high"}
	var first []byte
	for run := range 5 {
		blocks, err := RunBRM(cfg, "2p2b", brmInput(), brmTestThresholds)
		if err != nil {
			t.Fatal(err)
		}
		if len(blocks) != 3 || blocks[0].Chrom != "chr1" || blocks[1].Chrom != "chr3" || blocks[2].Chrom != "chr9" {
			t.Fatalf("blocks %+v", blocks)
		}
		tsv, err := os.ReadFile(filepath.Join(cfg.OutputDir, "stats", "GoBSAseq.2p2b.brm_blocks.tsv"))
		if err != nil {
			t.Fatal(err)
		}
		if run == 0 {
			first = tsv
		} else if string(tsv) != string(first) {
			t.Fatal("TSV differs between runs")
		}
	}
}

// brmTestThresholds: a pointwise-α-like u for unit tests (no null run).
var brmTestThresholds = Thresholds{BRMU: 1.96, BRMAlpha: 0.05}

func TestRunBRMErrors(t *testing.T) {
	color.Output = io.Discard
	base := utils.AnalysisConfig{OutputDir: t.TempDir(), Population: "F2", HighBulkSize: 20, LowBulkSize: 20, BrmAlpha: 0.05, RefParent: "high"}
	for name, cfg := range map[string]utils.AnalysisConfig{
		"no bulk size": func() utils.AnalysisConfig { c := base; c.LowBulkSize = 0; return c }(),
		"bad pop":      func() utils.AnalysisConfig { c := base; c.Population = "XYZ"; return c }(),
		"unwritable":   func() utils.AnalysisConfig { c := base; c.OutputDir = "/dev/null/x"; return c }(),
	} {
		if _, err := RunBRM(cfg, "2p2b", brmInput(), brmTestThresholds); err == nil {
			t.Errorf("%s: expected an error", name)
		}
	}
	if _, err := RunBRM(base, "2p2b", brmInput(), Thresholds{}); err == nil {
		t.Error("no genome-wide threshold: expected an error")
	}
	noVar := brmInput()
	noVar[7].ReadVar0 = 0
	if _, err := RunBRM(base, "2p2b", noVar, brmTestThresholds); err == nil {
		t.Error("marker without read variance: expected an error")
	}
	for _, a := range []float64{0, 1, math.NaN()} {
		c := base
		c.BrmAlpha = a
		if _, err := CalculateThresholds(c, "2p2b", []BSAstats{{CHROM: "c", POS: 1, HighDepth: 30, LowDepth: 30}}); err == nil {
			t.Errorf("--brm-alpha %v: expected an error", a)
		}
	}
	if err := writeBRMBlocks("/dev/full", make([]BRMBlock, 100000)); err == nil {
		t.Error("write to /dev/full: expected an error")
	}
}

// ReadVar0 + ρ·ReadVar1 is the variance of the smoothed ΔSI from read
// sampling: checked against the empirical variance over replicates with
// beta-binomial reads around a fixed pool frequency (no pool variation, huge
// bulks), in the interior, near chromosome ends (local linear blend) and in
// a sparse stretch, in grid mode with duplicate positions.
func TestReadVarMatchesSmoothedVariance(t *testing.T) {
	color.Output = io.Discard
	const rho, reps = 0.01, 400
	cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 1_000_000, LowBulkSize: 1_000_000}
	weight, err := markerWeight(cfg, "2p2b")
	if err != nil {
		t.Fatal(err)
	}
	readVar, err := readVariance(cfg, "2p2b")
	if err != nil {
		t.Fatal(err)
	}
	bw := float64(cfg.WindowSize) / 2
	sc := &scorer{bw: bw, step: gridStep(0, bw), weight: weight, readVar: readVar,
		lay: smoothLayout{highBulk: true, lowBulk: true, signedDSI: true}}
	tmpl := smoothInput(11, 1200)
	order := newMarkerOrder(tmpl)
	calc := rawCalc{roles: siteRoles{hp: -1, lp: -1, hb: 0, lb: 1}, signedDSI: true}
	rng := rand.New(rand.NewSource(7))
	draw := func(d int) int32 {
		p := 0.5 + math.Sqrt(rho*0.25)*rng.NormFloat64()
		return int32(binomSample(rng, d, math.Min(math.Max(p, 0), 1)))
	}
	type acc struct{ ss, pred float64 }
	type key struct {
		chrom string
		pos   int64
	}
	byPos := map[key]*acc{}
	chromEnd := map[string][2]int64{}
	raw := make([]BSAstats, len(tmpl))
	for r := range reps {
		for i, s := range tmpl {
			site := filter.Site{Chrom: s.CHROM, Pos: s.POS, Ref: "A", Alt: "T", NSamples: 2}
			hH, hL := draw(s.HighDepth), draw(s.LowDepth)
			site.Samples[0] = filter.SampleAD{Ref: hH, Alt: int32(s.HighDepth) - hH, HasAD: true}
			site.Samples[1] = filter.SampleAD{Ref: hL, Alt: int32(s.LowDepth) - hL, HasAD: true}
			raw[i] = calc.site(&site)
		}
		sm, err := sc.scoreOrdered(raw, order, false, 1)
		if err != nil {
			t.Fatal(err)
		}
		for _, s := range sm {
			k := key{s.CHROM, s.POS}
			a := byPos[k]
			if a == nil {
				a = &acc{}
				byPos[k] = a
			}
			a.ss += s.SmDeltaSI * s.SmDeltaSI
			if r == 0 {
				a.pred = s.ReadVar0 + rho*s.ReadVar1
				e := chromEnd[s.CHROM]
				if e[1] == 0 || s.POS < e[0] {
					e[0] = s.POS
				}
				e[1] = max(e[1], s.POS)
				chromEnd[s.CHROM] = e
			}
		}
	}
	// classes: near an end (< 2σ), sparse stretch, interior
	var obs, pred [3]float64
	var n [3]int
	for _, idx := range order.chroms {
		for _, i := range idx {
			s := tmpl[i]
			a := byPos[key{s.CHROM, s.POS}]
			if a == nil || a.pred == 0 {
				continue
			}
			e := chromEnd[s.CHROM]
			cls := 2
			switch {
			case float64(s.POS-e[0]) < 2*bw || float64(e[1]-s.POS) < 2*bw:
				cls = 0
			case s.POS > 8_000_000 && s.POS < 14_000_000:
				cls = 1
			}
			obs[cls] += a.ss / reps
			pred[cls] += a.pred
			n[cls]++
		}
	}
	for cls, name := range []string{"chromosome ends", "sparse stretch", "interior"} {
		if n[cls] == 0 {
			t.Fatalf("%s: no markers", name)
		}
		ratio := obs[cls] / pred[cls]
		t.Logf("%s: %d markers, empirical/predicted variance %.3f", name, n[cls], ratio)
		if math.Abs(ratio-1) > 0.08 {
			t.Errorf("%s: empirical/predicted smoothed read variance %.3f, want 1 ± 0.08", name, ratio)
		}
	}
}

// Merging joins same-sign blocks within the gap and keeps the larger peak
// (with its z); opposite signs and wider gaps stay separate.
func TestMergeBRMBlocks(t *testing.T) {
	b := []BRMBlock{
		{Start: 100, Stop: 200, PeakPos: 150, Peak: 0.3, PeakZ: 3},
		{Start: 250, Stop: 300, PeakPos: 280, Peak: 0.5, PeakZ: 5},
		{Start: 320, Stop: 400, PeakPos: 350, Peak: -0.4, PeakZ: -4},
		{Start: 900, Stop: 950, PeakPos: 920, Peak: -0.6, PeakZ: -6},
	}
	got := mergeBRMBlocks(b, 100)
	if len(got) != 3 || got[0].Stop != 300 || got[0].PeakPos != 280 || got[0].PeakZ != 5 || got[1].Start != 320 || got[2].Start != 900 {
		t.Fatalf("merged %+v", got)
	}
}

// --min-qtl-markers drops BRM blocks with fewer markers beyond the line.
func TestRunBRMMinMarkers(t *testing.T) {
	color.Output = io.Discard
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), Population: "F2", HighBulkSize: 20, LowBulkSize: 20, BrmAlpha: 0.05, RefParent: "high"}
	for _, c := range []struct{ min, want int }{{10, 3}, {11, 0}} { // each block has 10 markers
		cfg.MinQTLMarkers = c.min
		blocks, err := RunBRM(cfg, "2p2b", brmInput(), brmTestThresholds)
		if err != nil {
			t.Fatal(err)
		}
		if len(blocks) != c.want || (c.want > 0 && blocks[0].NMarkers != 10) {
			t.Fatalf("min %d: %+v", c.min, blocks)
		}
	}
}
