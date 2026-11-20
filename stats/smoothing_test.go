package stats

import (
	"io"
	"math"
	"os"
	"sort"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

func flatRaw(chrom string, n int, step int64, f func(i int) (si float64, depth int)) []BSAstats {
	raw := make([]BSAstats, n)
	for i := range raw {
		si, d := f(i)
		raw[i] = BSAstats{CHROM: chrom, POS: int64(i+1) * step, HighSI: si, LowSI: si, Depth: d, HighDepth: d, LowDepth: d}
	}
	return raw
}

func TestSmoothConstantStaysConstant(t *testing.T) {
	color.Output = io.Discard
	raw := flatRaw("chr1", 200, 10_000, func(i int) (float64, int) { return 0.3, 10 + i%50 })
	sm, err := gaussianSmooth(raw, 100_000, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	for _, s := range sm {
		if math.Abs(s.SmHighSI-0.3) > 1e-12 {
			t.Fatalf("constant 0.3 smoothed to %v at %d", s.SmHighSI, s.POS)
		}
	}
}

// A step from 0 to 1 at the middle of evenly spaced, equal-depth markers is
// smoothed symmetrically: values at mirror positions sum to 1.
func TestSmoothStepIsSymmetric(t *testing.T) {
	color.Output = io.Discard
	raw := flatRaw("chr1", 400, 10_000, func(i int) (float64, int) {
		if i >= 200 {
			return 1, 30
		}
		return 0, 30
	})
	sm, err := gaussianSmooth(raw, 200_000, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	for k := 0; k < 200; k++ {
		if s := sm[k].SmHighSI + sm[399-k].SmHighSI; math.Abs(s-1) > 1e-9 {
			t.Fatalf("asymmetric step at k=%d: %v + %v", k, sm[k].SmHighSI, sm[399-k].SmHighSI)
		}
	}
	if sm[199].SmHighSI >= 0.5 || sm[200].SmHighSI <= 0.5 {
		t.Fatalf("step midpoint: %v, %v", sm[199].SmHighSI, sm[200].SmHighSI)
	}
}

// depthWeight weights markers by depth, for kernel tests independent of
// markerWeight.
func depthWeight(s *BSAstats) float64 { return float64(s.Depth) }

// Marker weights enter the kernel: heavier markers dominate. (Markers sit a
// few bp apart, so their Gaussian distance weights are equal to ~1e-12.)
func TestSmoothMarkerWeighting(t *testing.T) {
	color.Output = io.Discard
	var raw []BSAstats
	for i := 0; i < 5; i++ {
		raw = append(raw, BSAstats{CHROM: "c", POS: int64(1000 + 2*i), HighSI: 1, Depth: 90}, BSAstats{CHROM: "c", POS: int64(1001 + 2*i), HighSI: 0, Depth: 10})
	}
	sm, err := gaussianSmooth(raw, 1e6, depthWeight, true, false, false, true, false)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(sm[0].SmHighSI-0.9) > 1e-9 {
		t.Fatalf("depth-weighted mean = %v, want 0.9", sm[0].SmHighSI)
	}
}

// Chromosomes are smoothed independently, and a contig with fewer than
// minEffectiveMarkers markers is dropped.
func TestSmoothChromosomesIndependentAndSparseDropped(t *testing.T) {
	color.Output = io.Discard
	raw := append(flatRaw("a", 50, 1000, func(int) (float64, int) { return 0, 20 }),
		flatRaw("b", 50, 1000, func(int) (float64, int) { return 1, 20 })...)
	raw = append(raw, flatRaw("tiny", int(minEffectiveMarkers)-1, 1000, func(int) (float64, int) { return 0.5, 20 })...)
	sm, err := gaussianSmooth(raw, 1e6, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	if len(sm) != 100 {
		t.Fatalf("got %d smoothed markers, want 100 (tiny contig dropped)", len(sm))
	}
	for _, s := range sm {
		if want := map[string]float64{"a": 0, "b": 1}[s.CHROM]; s.SmHighSI != want {
			t.Fatalf("%s:%d smoothed to %v, want %v (leak across chromosomes)", s.CHROM, s.POS, s.SmHighSI, want)
		}
	}
}

func TestRobustZScores(t *testing.T) {
	// median 3, |dev| = 2 1 0 1 2 → MAD 1 → scale 1.4826.
	z := robustZScores([]float64{1, 2, 3, 4, 5})
	for i, want := range []float64{-2 / 1.4826, -1 / 1.4826, 0, 1 / 1.4826, 2 / 1.4826} {
		if math.Abs(z[i]-want) > 1e-6 {
			t.Fatalf("z[%d] = %v, want %v", i, z[i], want)
		}
	}
	// MAD = 0 → classical Z (population SD).
	z = robustZScores([]float64{0, 0, 0, 0, 4})
	if math.Abs(z[4]-2) > 1e-6 || math.Abs(z[0]+0.5) > 1e-6 {
		t.Fatalf("MAD=0 fallback = %v, want [-0.5 … 2]", z)
	}
	if z := robustZScores([]float64{7, 7, 7}); z[0] != 0 || z[2] != 0 {
		t.Fatalf("constant input = %v, want zeros", z)
	}
}

// countZStats (used by the CompositeZ null) must match the number of Z-scores
// consolidate actually combines. With every component Z = 1, CompositeZ = √k.
func TestCompositeZStatCountMatchesConsolidate(t *testing.T) {
	for _, m := range allModes {
		hb, lb, both, one := BulkFlags(m.mode)
		for _, hmm := range []bool{true, false} {
			sm := []SmoothedStats{{ZDeltaSI: 1, ZGstat: 1, ZAFDev: 1, ZOneBulkG: 1, ZHMMLOD: 1}}
			consolidate(sm, hb, lb, both, one, true, hmm)
			k := int(math.Round(sm[0].CompositeZ * sm[0].CompositeZ))
			if k != countZStats(m.mode, hmm) {
				t.Errorf("%s hmm=%v: consolidate combines %d Z-scores, countZStats says %d", m.mode, hmm, k, countZStats(m.mode, hmm))
			}
		}
	}
}

func TestWriteSmoothedTSVReportsWriteErrors(t *testing.T) {
	if _, err := os.Stat("/dev/full"); err != nil {
		t.Skip("/dev/full not available")
	}
	sm := make([]SmoothedStats, 20000) // more output than the 1 MB buffer
	if err := writeSmoothedTSV("/dev/full", sm, smoothLayout{highBulk: true, lowBulk: true, signedDSI: true}); err == nil {
		t.Fatal("writing to a full device returned nil")
	}
}

func TestBulkFlagsExactModes(t *testing.T) {
	for _, m := range allModes {
		hb, lb, _, _ := BulkFlags(m.mode)
		wantHB := m.mode == "2p2b" || m.mode == "hp2b" || m.mode == "lp2b" || m.mode == "2b" || m.mode[len(m.mode)-2:] == "hb"
		wantLB := m.mode == "2p2b" || m.mode == "hp2b" || m.mode == "lp2b" || m.mode == "2b" || m.mode[len(m.mode)-2:] == "lb"
		if hb != wantHB || lb != wantLB {
			t.Errorf("%s: high/low bulk = %v/%v, want %v/%v", m.mode, hb, lb, wantHB, wantLB)
		}
	}
	if hb, lb, _, _ := BulkFlags("2b_typo_hb"); hb || lb {
		t.Error("unknown mode matched by substring")
	}
}

// Inverse-variance weights (F2, 20-plant bulks): pool variance per bulk is
// 1/(2·20), so a 200x marker is only 1.67× as precise as a 40x one.
func TestMarkerWeightInverseVariance(t *testing.T) {
	cfg := utils.AnalysisConfig{Population: "F2", HighBulkSize: 20, LowBulkSize: 20}
	two, err := markerWeight(cfg, "2p2b")
	if err != nil {
		t.Fatal(err)
	}
	r := two(&BSAstats{HighDepth: 200, LowDepth: 200}) / two(&BSAstats{HighDepth: 40, LowDepth: 40})
	if math.Abs(r-5.0/3) > 1e-9 {
		t.Fatalf("two-bulk weight ratio 200x/40x = %v, want 5/3", r)
	}
	// Unequal depths: the shallow bulk limits precision.
	if a, b := two(&BSAstats{HighDepth: 200, LowDepth: 10}), two(&BSAstats{HighDepth: 40, LowDepth: 40}); a >= b {
		t.Fatalf("200x/10x weight %v should be below 40x/40x weight %v", a, b)
	}
	one, err := markerWeight(cfg, "hphb")
	if err != nil {
		t.Fatal(err)
	}
	if r := one(&BSAstats{Depth: 200}) / one(&BSAstats{Depth: 40}); math.Abs(r-5.0/3) > 1e-9 {
		t.Fatalf("one-bulk weight ratio = %v, want 5/3", r)
	}
	if _, err := markerWeight(utils.AnalysisConfig{Population: "F2"}, "2p2b"); err == nil {
		t.Fatal("missing bulk sizes not reported")
	}
}

// An isolated marker 2.8σ from a cluster has 31 markers within ±3σ, but its
// effective number of markers is ~2.5, so it is dropped rather than smoothed
// to (essentially) its own value.
func TestSmoothDropsIsolatedMarkerByNEff(t *testing.T) {
	color.Output = io.Discard
	raw := flatRaw("c", 30, 1000, func(int) (float64, int) { return 0, 20 })
	raw = append(raw, BSAstats{CHROM: "c", POS: 30_000 + 2_800_000, HighSI: 1, LowSI: 1, Depth: 20})
	sm, err := gaussianSmooth(raw, 1e6, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	if len(sm) != 30 {
		t.Fatalf("got %d markers, want 30 (isolated marker dropped)", len(sm))
	}
	if sm[0].NEff < 29 || sm[0].NEff > 31 {
		t.Fatalf("cluster NEff = %v, want ~30", sm[0].NEff)
	}
}

// SI Z-scores enter MaxAbsZ only in one-bulk modes; two-bulk modes do not
// normalise SI (ΔSI carries the contrast).
func TestMaxAbsZSIOnlyInOneBulk(t *testing.T) {
	sm := []SmoothedStats{{ZHighSI: 10, ZLowSI: -9, ZDeltaSI: 1, ZGstat: 1}}
	for _, signed := range []bool{true, false} {
		consolidate(sm, true, true, true, false, signed, true)
		if sm[0].MaxAbsZ != 1 {
			t.Fatalf("two-bulk MaxAbsZ = %v (signed %v), want 1", sm[0].MaxAbsZ, signed)
		}
	}
	one := []SmoothedStats{{ZHighSI: 10, ZAFDev: 1}}
	consolidate(one, true, false, false, true, true, true)
	if one[0].MaxAbsZ != 10 {
		t.Fatalf("one-bulk MaxAbsZ = %v, want 10", one[0].MaxAbsZ)
	}
}

// Loess (local quadratic) reproduces a quadratic exactly, including at the
// chromosome ends, where a local average is biased.
func TestLoessReproducesQuadratic(t *testing.T) {
	color.Output = io.Discard
	n := 300
	raw := make([]BSAstats, n)
	f := func(p float64) float64 { x := p / 1e6; return 0.05 + 0.01*x + 0.002*x*x }
	for i := range raw {
		p := int64(i+1) * 10_000
		raw[i] = BSAstats{CHROM: "c", POS: p, Depth: 30, HighDepth: 30, LowDepth: 30, ED: f(float64(p))}
	}
	sm, err := gaussianSmooth(raw, 250_000, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	for _, s := range sm {
		if want := f(float64(s.POS)); math.Abs(s.SmED-want) > 1e-9 {
			t.Fatalf("loess at %d = %v, want %v", s.POS, s.SmED, want)
		}
	}
}

// G′: tricube-weighted mean over ±h (h = 2·bw) with marker weights; checked
// by hand at one marker.
func TestGPrimeTricube(t *testing.T) {
	color.Output = io.Discard
	bw := 1000.0 // h = 2000
	var raw []BSAstats
	for i, g := range []float64{1, 2, 3, 4, 5, 6, 7} {
		raw = append(raw, BSAstats{CHROM: "c", POS: int64(1000 + 500*i), Gstat: g, Depth: 10 + i})
	}
	sm, err := gaussianSmooth(raw, bw, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	// Query at POS 2500 (index 3): neighbours within |Δ| < 2000.
	num, den := 0.0, 0.0
	for _, r := range raw {
		d := math.Abs(float64(r.POS-2500)) / 2000
		if d >= 1 {
			continue
		}
		w := float64(r.Depth) * math.Pow(1-d*d*d, 3)
		num += w * r.Gstat
		den += w
	}
	var got *SmoothedStats
	for i := range sm { // edge markers may be dropped by the NEff rule
		if sm[i].POS == 2500 {
			got = &sm[i]
		}
	}
	if got == nil {
		t.Fatal("marker at 2500 was dropped")
	}
	if math.Abs(got.SmGstat-num/den) > 1e-12 {
		t.Fatalf("G′ = %v, want %v", got.SmGstat, num/den)
	}
}

// The loess fit falls back gracefully and smoothed ED⁴ is never negative.
func TestLoessFallbackAndClamp(t *testing.T) {
	if v := loessFit(1, 0, 0, 0, 0, 0.3, 0, 0, 1); math.Abs(v-0.3) > 1e-12 {
		t.Fatalf("single point: %v, want 0.3", v)
	}
	if v := loessFit(2, 0, 2, 0, 2, 1, 0, 1, 2); math.IsNaN(v) || math.IsInf(v, 0) {
		t.Fatalf("two points: %v", v)
	}
	color.Output = io.Discard
	// A sharp dip to 0 next to large values: a quadratic fit undershoots.
	var raw []BSAstats
	for i := range 60 {
		ed := 1.0
		if i >= 28 && i <= 32 {
			ed = 0
		}
		raw = append(raw, BSAstats{CHROM: "c", POS: int64(1 + 1000*i), ED: ed, Depth: 30, HighDepth: 30, LowDepth: 30})
	}
	sm, err := gaussianSmooth(raw, 3000, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	for _, s := range sm {
		if s.SmED < 0 {
			t.Fatalf("negative smoothed ED⁴ %v at %d", s.SmED, s.POS)
		}
	}
}

// robustZField (in place, one buffer) gives exactly robustZScores' values.
func TestRobustZFieldMatchesRobustZScores(t *testing.T) {
	for _, vals := range [][]float64{{1, 5, 2, 8, 3, 3, 9, 0.5}, {0, 0, 0, 0, 4}, {7, 7, 7}} {
		sm := make([]SmoothedStats, len(vals))
		for i, v := range vals {
			sm[i].SmGstat = v
		}
		robustZField(sm, zField{func(s *SmoothedStats) float64 { return s.SmGstat }, func(s *SmoothedStats, v float64) { s.ZGstat = v }}, make([]float64, len(sm)))
		want := robustZScores(vals)
		for i := range sm {
			if sm[i].ZGstat != want[i] {
				t.Fatalf("%v: z[%d] = %v, want %v", vals, i, sm[i].ZGstat, want[i])
			}
		}
	}
}

// Gaussian tracks switch to a local linear fit where a chromosome end cuts
// the kernel off: a linear trend is reproduced exactly within 2σ of each end
// (a local average is pulled toward the interior there), and to within
// 0.1% of the trend in the 2–3σ blend zone and interior (the old local
// average was off by ~8% of the trend at the chromosome end).
func TestLocalLinearNoEdgeBias(t *testing.T) {
	color.Output = io.Discard
	n, bw := 400, 200_000.0
	raw := make([]BSAstats, n)
	for i := range raw {
		si := 0.2 + 0.6*float64(i)/float64(n-1) // 0.2 → 0.8 along the chromosome
		raw[i] = BSAstats{CHROM: "c", POS: int64(i+1) * 5_000, HighSI: si, LowSI: si, Depth: 30, HighDepth: 30, LowDepth: 30}
	}
	sm, err := gaussianSmooth(raw, bw, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	first, last := float64(raw[0].POS), float64(raw[n-1].POS)
	for k, s := range sm {
		err := math.Abs(s.SmHighSI - raw[k].HighSI)
		edge := math.Min(float64(s.POS)-first, last-float64(s.POS)) <= 2*bw
		if (edge && err > 1e-9) || err > 0.001*0.6 {
			t.Fatalf("SI at %d = %v, want %v (error %.2g, edge zone %v)", s.POS, s.SmHighSI, raw[k].HighSI, err, edge)
		}
	}
	if edgeBlend(0, bw) != 1 || edgeBlend(3*bw, bw) != 0 || edgeBlend(2.5*bw, bw) != 0.5 {
		t.Fatal("edgeBlend ramp")
	}
}

// Grid evaluation on a dense chromosome (grid points, interpolation, and
// neighbours aggregated into σ/100 bins) agrees with exact evaluation at every
// marker to within 1% of each track's range. The input is pure per-site
// noise, the worst case for aggregation; smooth signals agree far more closely.
func TestGridMatchesExact(t *testing.T) {
	color.Output = io.Discard
	raw := smoothInput(11, 20000) // dense, with a sparse stretch at 8–14 Mb
	o := smoothOpts{bw: 1e6, weight: depthWeight, hasHighBulk: true, hasLowBulk: true, signedDSI: true, every: 1}
	exact, err := gaussianSmoothOpts(raw, o)
	if err != nil {
		t.Fatal(err)
	}
	o.step = gridStep(100_000, o.bw)
	grid, err := gaussianSmoothOpts(raw, o)
	if err != nil {
		t.Fatal(err)
	}
	if len(exact) != len(grid) {
		t.Fatalf("%d vs %d markers", len(exact), len(grid))
	}
	tracks := map[string]func(s SmoothedStats) float64{
		"HighSI": func(s SmoothedStats) float64 { return s.SmHighSI }, "DeltaSI": func(s SmoothedStats) float64 { return s.SmDeltaSI },
		"G′": func(s SmoothedStats) float64 { return s.SmGstat }, "ED4": func(s SmoothedStats) float64 { return s.SmED },
	}
	for name, get := range tracks {
		lo, hi, worst := math.Inf(1), math.Inf(-1), 0.0
		for i := range exact {
			v := get(exact[i])
			lo, hi = math.Min(lo, v), math.Max(hi, v)
			worst = math.Max(worst, math.Abs(v-get(grid[i])))
		}
		if worst > 0.01*(hi-lo) {
			t.Errorf("%s: max |grid − exact| %.3g exceeds 1%% of range %.3g", name, worst, hi-lo)
		}
	}
}

// Markers sharing a position (split multi-allelic views) count as one
// marker: two views with the same value smooth exactly like one marker.
func TestDuplicatePositionsShareWeight(t *testing.T) {
	color.Output = io.Discard
	base := flatRaw("c", 60, 20_000, func(i int) (float64, int) { return 0, 30 })
	one := append([]BSAstats(nil), base...)
	one[30].HighSI, one[30].LowSI = 1, 1
	two := append(append([]BSAstats(nil), one...), one[30]) // second view at the same position
	a, err := gaussianSmooth(one, 200_000, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	b, err := gaussianSmooth(two, 200_000, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	for i := range a {
		if math.Abs(a[i].SmHighSI-b[i].SmHighSI) > 1e-12 || math.Abs(a[i].NEff-b[i].NEff) > 1e-9 {
			t.Fatalf("marker %d: one view %v (NEff %v), two views %v (NEff %v)", a[i].POS, a[i].SmHighSI, a[i].NEff, b[i].SmHighSI, b[i].NEff)
		}
	}
}

// Contigs shorter than two kernel widths (6σ) keep the local average: no
// local linear fit near their ends.
func TestShortContigNoEdgeFit(t *testing.T) {
	color.Output = io.Discard
	bw := 200_000.0 // 6σ = 1.2 Mb; the contig spans 1 Mb
	var raw []BSAstats
	for i := range 101 {
		si := float64(i) / 100
		raw = append(raw, BSAstats{CHROM: "short", POS: int64(1 + 10_000*i), HighSI: si, LowSI: si, Depth: 30, HighDepth: 30, LowDepth: 30})
	}
	sm, err := gaussianSmooth(raw, bw, depthWeight, true, true, true, false, false)
	if err != nil {
		t.Fatal(err)
	}
	// At the first marker the local average of a rising trend is pulled up.
	num, den := 0.0, 0.0
	for _, r := range raw {
		d := float64(r.POS-raw[0].POS) / bw
		if d > 3 {
			continue
		}
		w := 30 * math.Exp(-0.5*d*d)
		num += w * r.HighSI
		den += w
	}
	if math.Abs(sm[0].SmHighSI-num/den) > 1e-12 {
		t.Fatalf("short contig first marker = %v, want local average %v", sm[0].SmHighSI, num/den)
	}
}

// robustZScores is an independent reference implementation of the robust
// Z-score (copies and sorts its own data), used to check robustZField.
func robustZScores(vals []float64) []float64 {
	n := len(vals)
	if n == 0 {
		return nil
	}
	median := func(v []float64) float64 {
		s := append([]float64(nil), v...)
		sort.Float64s(s)
		if len(s)%2 == 1 {
			return s[len(s)/2]
		}
		return (s[len(s)/2-1] + s[len(s)/2]) / 2
	}
	med := median(vals)
	dev := make([]float64, n)
	for i, v := range vals {
		dev[i] = math.Abs(v - med)
	}
	mad := median(dev) * 1.4826
	zs := make([]float64, n)
	if mad > 0 {
		for i, v := range vals {
			zs[i] = (v - med) / mad
		}
		return zs
	}
	mean, variance := 0.0, 0.0
	for _, v := range vals {
		mean += v
	}
	mean /= float64(n)
	for _, v := range vals {
		variance += (v - mean) * (v - mean)
	}
	if sd := math.Sqrt(variance / float64(n)); sd > 0 {
		for i, v := range vals {
			zs[i] = (v - mean) / sd
		}
	}
	return zs
}

// Markers all at one position (a multi-allelic site split on a tiny contig)
// with a grid step set: no crash. They are one site (NEff 1), so all are
// dropped as too sparse.
func TestSmoothSamePositionMarkers(t *testing.T) {
	color.Output = io.Discard
	var raw []BSAstats
	for range 6 {
		raw = append(raw, BSAstats{CHROM: "tiny", POS: 5000, HighSI: 0.5, LowSI: 0.5, Depth: 30, HighDepth: 30, LowDepth: 30})
	}
	sm, err := gaussianSmoothOpts(raw, smoothOpts{bw: 1e6, step: 50_000, weight: depthWeight, hasHighBulk: true, hasLowBulk: true, signedDSI: true, every: 1})
	if err != nil {
		t.Fatal(err)
	}
	if len(sm) != 0 {
		t.Fatalf("got %d markers, want 0 (one site, NEff 1)", len(sm))
	}
}

// Compact null output (every > 1) holds exactly the sampled markers with the
// same values as the full output.
func TestSmoothCompactOutput(t *testing.T) {
	color.Output = io.Discard
	raw := smoothInput(5, 3000)
	o := smoothOpts{bw: 1e6, step: 50_000, weight: depthWeight, hasHighBulk: true, hasLowBulk: true, signedDSI: true, every: 1}
	full, err := gaussianSmoothOpts(raw, o)
	if err != nil {
		t.Fatal(err)
	}
	byKey := map[[2]int64]SmoothedStats{}
	for _, s := range full {
		byKey[[2]int64{int64(len(s.CHROM)), s.POS}] = s // CHROM lengths differ ("chr2" vs "chr10")
	}
	o.every = 7
	compact, err := gaussianSmoothOpts(raw, o)
	if err != nil {
		t.Fatal(err)
	}
	if cap(compact) >= len(raw) || len(compact) == 0 {
		t.Fatalf("compact output len %d cap %d, raw %d", len(compact), cap(compact), len(raw))
	}
	for _, s := range compact {
		f, ok := byKey[[2]int64{int64(len(s.CHROM)), s.POS}]
		if !ok || f.SmDeltaSI != s.SmDeltaSI || f.SmGstat != s.SmGstat || f.NEff != s.NEff {
			t.Fatalf("compact marker %s:%d differs from full output", s.CHROM, s.POS)
		}
	}
}
