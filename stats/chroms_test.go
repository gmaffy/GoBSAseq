package stats

import "testing"

func TestByChromosome(t *testing.T) {
	sm := []SmoothedStats{
		{CHROM: "chr2", POS: 1}, {CHROM: "chr2", POS: 5}, {CHROM: "chr2", POS: 5},
		{CHROM: "chr1", POS: 3}, {CHROM: "chr1", POS: 9},
	}
	got := byChromosome(sm)
	if len(got) != 2 || got[0].chrom != "chr1" || got[1].chrom != "chr2" {
		t.Fatalf("chromosomes %+v", got)
	}
	if &got[0].sm[0] != &sm[3] || &got[1].sm[0] != &sm[0] || len(got[1].sm) != 3 {
		t.Fatal("contiguous, sorted chromosomes must be subslices of the input, not copies")
	}
	// No record copies for sorted input: allocations do not grow with markers.
	big := make([]SmoothedStats, 200000)
	for i := range big {
		big[i] = SmoothedStats{CHROM: []string{"a", "b", "c"}[i*3/len(big)], POS: int64(i)}
	}
	if a := testing.AllocsPerRun(5, func() { byChromosome(big) }); a > 20 {
		t.Fatalf("%v allocations for sorted input", a)
	}

	// Split or unsorted chromosomes are gathered and sorted (stable).
	sm = []SmoothedStats{
		{CHROM: "c", POS: 7}, {CHROM: "d", POS: 1}, {CHROM: "c", POS: 2, REF: "x"}, {CHROM: "c", POS: 2, REF: "y"},
		{CHROM: "e", POS: 4}, {CHROM: "e", POS: 3},
	}
	got = byChromosome(sm)
	c, e := got[0].sm, got[2].sm
	if len(c) != 3 || c[0].REF != "x" || c[1].REF != "y" || c[2].POS != 7 || e[0].POS != 3 || e[1].POS != 4 {
		t.Fatalf("c %+v e %+v", c, e)
	}
	if &c[0] == &sm[2] {
		t.Fatal("a split chromosome must be a copy")
	}
}
