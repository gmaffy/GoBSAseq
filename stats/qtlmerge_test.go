package stats

import "testing"

type run struct {
	chrom      string
	start, end float64
	tag        string
	n          int
}

func runSpan(r *run) (string, float64, float64) { return r.chrom, r.start, r.end }

func joinRun(a *run, b run) {
	a.end = max(a.end, b.end)
	a.n += b.n
}

// mergeRuns joins items within gap on one chromosome, in sorted order, and
// respects canJoin between consecutive items.
func TestMergeRuns(t *testing.T) {
	in := []run{
		{"b", 0, 10, "x", 1}, {"a", 50, 60, "x", 1}, {"a", 0, 10, "x", 1}, {"a", 15, 30, "x", 1},
		{"a", 100, 110, "y", 1}, {"a", 112, 120, "x", 1},
	}
	got := mergeRuns(in, 5, runSpan, nil, joinRun)
	want := []run{{"a", 0, 30, "x", 2}, {"a", 50, 60, "x", 1}, {"a", 100, 120, "y", 2}, {"b", 0, 10, "x", 1}}
	if len(got) != len(want) {
		t.Fatalf("got %+v", got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("got %+v, want %+v", got, want)
		}
	}
	// canJoin: only same tag; consecutive items only (the y run blocks x–x)
	got = mergeRuns(in, 5, runSpan, func(a, b *run) bool { return a.tag == b.tag }, joinRun)
	if len(got) != 5 || got[2].tag != "y" || got[3].start != 112 {
		t.Fatalf("canJoin: %+v", got)
	}
	if len(in) != 6 || in[0].chrom != "b" {
		t.Fatal("input was modified")
	}
}
