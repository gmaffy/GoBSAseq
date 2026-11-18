package stats

import "testing"

// A triangular LOD profile peaking at 50 kb (LOD 10, falling 1 per 1 kb):
// a 1.5-LOD drop keeps markers with LOD ≥ 8.5, i.e. 49–51 kb.
func triangle() []SmoothedStats {
	var sm []SmoothedStats
	for i := 0; i <= 100; i++ {
		d := i - 50
		if d < 0 {
			d = -d
		}
		sm = append(sm, SmoothedStats{CHROM: "c", POS: int64(i * 1000), SmLOD: 10 - float64(d)})
	}
	return sm
}

func TestSupportInChromLODDrop(t *testing.T) {
	sm := triangle()
	lod := supportLOD("2p2b")

	got := supportInChrom(sm, lod, 0, 100000, 1.5)
	want := SupportInterval{Start: 49000, End: 51000, PeakPos: 50000, PeakLOD: 10}
	if got != want {
		t.Fatalf("got %+v, want %+v", got, want)
	}

	// The peak is searched within the region only, but the interval may
	// extend past the region's edge.
	got = supportInChrom(sm, lod, 20000, 49500, 3)
	want = SupportInterval{Start: 46000, End: 54000, PeakPos: 49000, PeakLOD: 9}
	if got != want {
		t.Fatalf("restricted region: got %+v, want %+v", got, want)
	}

	if s := supportInChrom(sm, lod, 200000, 300000, 1.5); s.Valid() {
		t.Fatalf("region with no markers should give no interval, got %+v", s)
	}
}

func TestSupportUnion(t *testing.T) {
	a := SupportInterval{Start: 10, End: 20, PeakPos: 15, PeakLOD: 5}
	b := SupportInterval{Start: 18, End: 40, PeakPos: 30, PeakLOD: 8}
	got := a.union(b)
	want := SupportInterval{Start: 10, End: 40, PeakPos: 30, PeakLOD: 8}
	if got != want {
		t.Fatalf("got %+v, want %+v", got, want)
	}
	if a.union(SupportInterval{}) != a || (SupportInterval{}).union(a) != a {
		t.Fatal("union with an empty interval should return the other one")
	}
}
