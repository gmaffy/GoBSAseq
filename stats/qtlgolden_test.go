package stats

import (
	"bytes"
	"io"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// qtlGoldenInput: three chromosomes of smooth bumps (Gaussian in position):
// chr1 a "+" QTL supported by every statistic, chr2 a "−" QTL, chr3 an HMM
// LOD plateau with a small DsiGconsensusZ blip elsewhere (partial support).
func qtlGoldenInput() []SmoothedStats {
	bump := func(x, mu, sd float64) float64 { return math.Exp(-0.5 * (x - mu) * (x - mu) / (sd * sd)) }
	var sm []SmoothedStats
	for _, c := range []string{"chr1", "chr2", "chr3"} {
		for i := range 400 {
			x := float64(50_000 * (i + 1))
			s := SmoothedStats{CHROM: c, POS: int64(x), Depth: 60, HighDepth: 60, LowDepth: 60, NEff: 50,
				ReadVar0: 2e-4, ReadVar1: 2e-2, SmHighSI: 0.5, SmLowSI: 0.5}
			switch c {
			case "chr1":
				b := bump(x, 8e6, 1.5e6)
				s.SmDeltaSI, s.SmGstat, s.SmED, s.HMMLOD, s.DsiGconsensusZ = 0.5*b, 20*b, 0.3*b*b, 15*b, 25*b
			case "chr2":
				b := bump(x, 12e6, 1e6)
				s.SmDeltaSI, s.SmGstat, s.SmED, s.HMMLOD, s.DsiGconsensusZ = -0.45*b, 16*b, 0.2*b*b, 12*b, 18*b
			case "chr3":
				if x > 4e6 && x < 7e6 {
					s.HMMLOD = 9 + 0.2*math.Sin(x/1e5) // plateau
				}
				s.DsiGconsensusZ = 14 * bump(x, 15e6, 3e5)
			}
			s.SmHighSI, s.SmLowSI = 0.5+s.SmDeltaSI/2, 0.5-s.SmDeltaSI/2
			sm = append(sm, s)
		}
	}
	return sm
}

// TestQTLGolden pins the three QTL lists (TSV bytes) for a fixed input.
// Regenerate with -update-golden only for intended changes.
func TestQTLGolden(t *testing.T) {
	color.Output = io.Discard
	sm := qtlGoldenInput()
	var th Thresholds
	tb := &th.TwoBulk
	tb.DeltaSIP99, tb.DeltaSIP95, tb.DeltaSIMp99, tb.DeltaSIMp95 = 0.3, 0.25, -0.3, -0.25
	tb.GstatP99, tb.GstatP95, tb.ED4P99, tb.ED4P95 = 10, 8, 0.1, 0.08
	th.HMM.LODP99, th.HMM.LODP95 = 5, 4
	th.Consensus.P99, th.Consensus.P95 = 12, 9
	th.BRMU, th.BRMAlpha, th.Rho = 3, 0.05, 0.01
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), Population: "F2", HighBulkSize: 20, LowBulkSize: 20, BrmAlpha: 0.05,
		RefParent: "high", MinQTLMarkers: 5}
	gap, minM := int64(1_000_000), cfg.EffectiveMinQTLMarkers()

	peaks, ind := DetectIndividualQTLs(sm, th, "2p2b", true, gap, minM)
	cz := DetectDsiGconsensusZQTLs(sm, th, gap, minM)
	brm, err := RunBRM(cfg, "2p2b", sm, th)
	if err != nil {
		t.Fatal(err)
	}
	final := ConsensusQTLs(sm, th, brm, peaks, true, true, gap, minM)

	dir := t.TempDir()
	write := map[string]func(string) error{
		"individual":     func(p string) error { return WriteIndividualQTLs(ind, "2p2b", p) },
		"dsigconsensusz": func(p string) error { return WriteDsiGconsensusZQTLs(cz, p) },
		"final":          func(p string) error { return WriteFinalQTLs(final, p) },
	}
	for name, w := range write {
		if err := w(filepath.Join(dir, name+".xlsx")); err != nil {
			t.Fatal(err)
		}
		got, err := os.ReadFile(filepath.Join(dir, name+".tsv"))
		if err != nil {
			t.Fatal(err)
		}
		path := filepath.Join("testdata", "qtl_golden", name+".tsv")
		if *updateGolden {
			if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
				t.Fatal(err)
			}
			if err := os.WriteFile(path, got, 0o644); err != nil {
				t.Fatal(err)
			}
			continue
		}
		want, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("%v (run with -update-golden to create)", err)
		}
		if !bytes.Equal(got, bytes.ReplaceAll(want, []byte("\r\n"), []byte("\n"))) { // ignore CRLF from Windows checkouts
			t.Errorf("%s.tsv differs from golden:\n%s", name, got)
		}
	}
}

// Each QTL reports its support count and a lead: the HMM LOD support
// midpoint when present, else the statistic furthest beyond its threshold.
func TestIndividualLead(t *testing.T) {
	q := ConsolidatedQTL{Peaks: [nSlots]*Peak{
		slotGstat:   {Pos: 100, Value: 20, Threshold: 10, Stat: "Gstat"},        // 2×
		slotDeltaSI: {Pos: 150, Value: -0.9, Threshold: -0.3, Stat: "DeltaSI-"}, // 3×
	}}
	setLead(&q, nil)
	if q.NStats != 2 || q.LeadStat != "DeltaSI-" || q.LeadPos != 150 {
		t.Fatalf("no HMM: %+v", q)
	}
	sm := make([]SmoothedStats, 10)
	for i := range sm {
		sm[i] = SmoothedStats{POS: int64(100 * (i + 1)), HMMLOD: []float64{0, 1, 7, 8, 8.5, 8, 7.5, 2, 0, 0}[i]}
	}
	q.Peaks[slotHMMLOD] = &Peak{Pos: 500, Value: 8.5, Threshold: 5, Stat: "HMMLOD"}
	setLead(&q, sm)
	// support within 1.5 LOD of 8.5: 300–700, midpoint 500
	if q.NStats != 3 || q.LeadStat != "HMMLOD" || q.LeadPos != 500 {
		t.Fatalf("with HMM: %+v", q)
	}
}
