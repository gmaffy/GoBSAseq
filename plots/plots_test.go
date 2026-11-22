package plots

import (
	"io"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
)

// interp is the simplified line's value at x (pts sorted by position).
func interp(pts [][2]float64, x float64) float64 {
	for i := 1; i < len(pts); i++ {
		if x <= pts[i][0] {
			a, b := pts[i-1], pts[i]
			if b[0] == a[0] {
				return b[1]
			}
			return a[1] + (b[1]-a[1])*(x-a[0])/(b[0]-a[0])
		}
	}
	return pts[len(pts)-1][1]
}

// Every marker lies within the tolerance of the plotted line; the ends and
// the extremes are kept exactly; non-finite values are dropped; a smooth
// track needs far fewer points than a jagged one.
func TestSimplify(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	n := 20000
	x := make([]int64, n)
	for i := range x {
		x[i] = int64(1000 + 1000*i + rng.Intn(500))
	}
	walk, smooth := make([]float64, n), make([]float64, n)
	for i := range walk {
		if i > 0 {
			walk[i] = walk[i-1] + rng.NormFloat64()
		}
		smooth[i] = math.Sin(float64(x[i]) / 3e6)
	}
	walk[777] = math.NaN()
	sizes := map[string]int{}
	for name, y := range map[string][]float64{"walk": walk, "smooth": smooth} {
		pts := simplify(x, y, 5)
		sizes[name] = len(pts)
		lo, hi := math.Inf(1), math.Inf(-1)
		for _, v := range y {
			if finite(v) {
				lo, hi = math.Min(lo, v), math.Max(hi, v)
			}
		}
		tol := plotTolerance*(hi-lo) + 1e-5 // + rounding to 5 decimals
		for i := 1; i < len(pts); i++ {
			if pts[i][0] < pts[i-1][0] {
				t.Fatalf("%s: positions not sorted at %d", name, i)
			}
		}
		for i, v := range y {
			if !finite(v) {
				continue
			}
			if d := math.Abs(interp(pts, float64(x[i])) - v); d > tol {
				t.Fatalf("%s: marker %d is %.3g from the plotted line (tolerance %.3g)", name, i, d, tol)
			}
		}
		if pts[0][0] != float64(x[0]) || pts[len(pts)-1][0] != float64(x[n-1]) {
			t.Fatalf("%s: ends not kept", name)
		}
		var pmin, pmax float64 = math.Inf(1), math.Inf(-1)
		for _, p := range pts {
			pmin, pmax = math.Min(pmin, p[1]), math.Max(pmax, p[1])
			if math.IsNaN(p[1]) {
				t.Fatalf("%s: NaN plotted", name)
			}
		}
		if math.Abs(pmax-hi) > 1e-5 || math.Abs(pmin-lo) > 1e-5 {
			t.Fatalf("%s: extremes %v..%v plotted as %v..%v", name, lo, hi, pmin, pmax)
		}
	}
	if !(sizes["smooth"] < n/50 && sizes["smooth"] < sizes["walk"]) {
		t.Fatalf("points kept: %v of %d", sizes, n)
	}
	if simplify(nil, nil, 5) != nil {
		t.Fatal("empty track should plot nothing")
	}
}

// The pages use a base-pair value axis, straight segments, thresholds as mark
// lines (both tails), JavaScript without escaped double quotes, and stay small.
func TestGeneratePlotsPages(t *testing.T) {
	color.Output = io.Discard
	var sm []stats.SmoothedStats
	for _, c := range []string{"chr1", "chr2"} {
		for i := range 5000 {
			p := float64(i) / 5000
			sm = append(sm, stats.SmoothedStats{CHROM: c, POS: int64(1000 + 4000*i),
				SmHighSI: 0.5 + 0.2*math.Sin(6*p), SmLowSI: 0.5 - 0.2*math.Sin(6*p), SmDeltaSI: 0.4 * math.Sin(6*p),
				SmGstat: 10 + 5*math.Sin(6*p), SmED: 0.01, HMMLOD: 3 + 3*math.Sin(6*p), DsiGconsensusZ: 2 * math.Sin(6*p)})
		}
	}
	th := stats.Thresholds{HMM: stats.HMMThresholds{LODP99: 5, LODP95: 4}, Consensus: stats.ConsensusThresholds{P99: 4, P95: 3, R: 0.9}}
	th.TwoBulk = stats.TwoBulkThresholds{DeltaSIP99: 0.3, DeltaSIP95: 0.25, DeltaSIMp99: -0.3, DeltaSIMp95: -0.25,
		HighSIP99: 0.7, HighSIP95: 0.65, HighSIMp99: 0.3, HighSIMp95: 0.35, LowSIP99: 0.7, LowSIP95: 0.65, LowSIMp99: 0.3, LowSIMp95: 0.35,
		GstatP99: 14, GstatP95: 12, ED4P99: 0.05, ED4P95: 0.04}
	cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), RefParent: "high", WindowSize: 2_000_000}
	brm := []stats.BRMBlock{{Chrom: "chr1", Start: 2_000_000, Stop: 3_000_000, Peak: 0.4}}
	if err := GeneratePlots(cfg, "2p2b", sm, th, brm); err != nil {
		t.Fatal(err)
	}
	for _, name := range []string{"individual_plots", "dsigconsensusz"} {
		b, err := os.ReadFile(filepath.Join(cfg.OutputDir, "plots", "GoBSAseq.2p2b."+name+".html"))
		if err != nil {
			t.Fatal(err)
		}
		h := string(b)
		if strings.Contains(h, `"smooth":true`) || !strings.Contains(h, `"type":"value"`) || !strings.Contains(h, `"min":"dataMin"`) {
			t.Errorf("%s: want a linear value axis and no spline", name)
		}
		if !strings.Contains(h, `"markLine":{"data":[{"name":"p99`) {
			t.Errorf("%s: no p99 mark line", name)
		}
		if name == "individual_plots" && !strings.Contains(h, `"name":"p99 valley"`) {
			t.Errorf("%s: no lower-tail line", name)
		}
		if regexp.MustCompile(`=== \\"`).MatchString(h) {
			t.Errorf("%s: escaped double quotes in JavaScript", name)
		}
		if len(b) > 400_000 {
			t.Errorf("%s: %d bytes for 10,000 markers", name, len(b))
		}
	}
}
