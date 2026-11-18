package stats

import (
	"testing"

	"github.com/gmaffy/GoBSAseq/utils"
)

func TestRegionFilter(t *testing.T) {
	var sm []SmoothedStats
	for p := int64(0); p <= 10_000_000; p += 100_000 { // one marker per 100 kb
		sm = append(sm, SmoothedStats{CHROM: "c1", POS: p})
	}
	cfg := utils.AnalysisConfig{WindowSize: 2_000_000, MinQTLMarkers: 5} // width 0 → σ = 1 Mb
	f := NewRegionFilter(cfg, sm)

	cases := []struct {
		start, end float64
		want       bool
	}{
		{1_000_000, 3_000_000, true},  // 2 Mb, 21 markers
		{1_000_000, 1_500_000, false}, // narrower than σ
		{1_000_000, 2_000_000, true},  // exactly σ, 11 markers
	}
	for _, c := range cases {
		if got := f.Keep("c1", c.start, c.end); got != c.want {
			t.Errorf("Keep(%v, %v) = %v, want %v", c.start, c.end, got, c.want)
		}
	}

	// Wide enough but too few markers (sparse region on another chromosome).
	sparse := append(sm, SmoothedStats{CHROM: "c2", POS: 0}, SmoothedStats{CHROM: "c2", POS: 3_000_000})
	if NewRegionFilter(cfg, sparse).Keep("c2", 0, 3_000_000) {
		t.Error("a 3 Mb region with 2 markers should be dropped (min 5)")
	}

	// Negative width and zero markers disable both criteria; a nil filter keeps all.
	off := NewRegionFilter(utils.AnalysisConfig{WindowSize: 2_000_000, MinQTLWidth: -1}, sm)
	if !off.Keep("c1", 1_000_000, 1_000_001) {
		t.Error("disabled filter should keep everything")
	}
	var none *RegionFilter
	if !none.Keep("c1", 0, 1) {
		t.Error("nil filter should keep everything")
	}
}

// chromosomes must share memory for sorted, contiguous input and still group
// correctly (in order of first appearance) when chromosomes are interleaved.
func TestChromosomes(t *testing.T) {
	sm := []SmoothedStats{{CHROM: "a", POS: 1}, {CHROM: "a", POS: 2}, {CHROM: "b", POS: 5}}
	g := chromosomes(sm)
	if len(g) != 2 || len(g[0]) != 2 || len(g[1]) != 1 || &g[0][0] != &sm[0] {
		t.Fatalf("contiguous input should be split without copying: %+v", g)
	}
	mixed := []SmoothedStats{{CHROM: "b", POS: 9}, {CHROM: "a", POS: 3}, {CHROM: "b", POS: 4}}
	g = chromosomes(mixed)
	if len(g) != 2 || g[0][0].CHROM != "b" || g[0][0].POS != 4 || g[0][1].POS != 9 || g[1][0].CHROM != "a" {
		t.Fatalf("interleaved input grouped wrongly: %+v", g)
	}
}
