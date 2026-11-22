package stats

import (
	"io"
	"math"
	"math/rand"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"gonum.org/v1/gonum/stat/distuv"
)

// normalQuantile matches gonum in the body and inverts the exact upper tail
// (math.Erfc) far beyond where 1 − p is representable.
func TestNormalQuantile(t *testing.T) {
	for _, p := range []float64{1e-6, 0.001, 0.02, 0.2, 0.5, 0.7, 0.975, 0.999} {
		if got, want := normalQuantile(p), distuv.UnitNormal.Quantile(p); math.Abs(got-want) > 1e-9*math.Max(1, math.Abs(want)) {
			t.Fatalf("Φ⁻¹(%v) = %v, want %v", p, got, want)
		}
	}
	for _, z := range []float64{5, 10, 20, 35} {
		q := 0.5 * math.Erfc(z/math.Sqrt2) // P(Z > z)
		if got := -normalQuantile(q); math.Abs(got-z) > 1e-9*z {
			t.Fatalf("upper tail %g: z = %v, want %v", q, got, z)
		}
		if got := normalUpperLog(-math.Log(q)); math.Abs(got-z) > 1e-9*z {
			t.Fatalf("normalUpperLog: z = %v, want %v", got, z)
		}
	}
	// beyond float64: the Mills-ratio expansion continues the tail
	if a, b := normalUpperLog(699), normalUpperLog(701); !(b > a && b-a < 0.1) {
		t.Fatalf("tail not continuous at the switch: %v, %v", a, b)
	}
}

// NullCDF scores of null samples are N(0, 1); a mass of ties gets its
// midpoint; the exponential tail is continuous and increasing.
func TestNullCDFScores(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	sample := func(n int, draw func() float64) []float64 {
		x := make([]float64, n)
		for i := range x {
			x[i] = draw()
		}
		return x
	}
	for name, draw := range map[string]func() float64{
		"chi2_1": func() float64 { z := rng.NormFloat64(); return z * z },
		"skewed": func() float64 { return rng.ExpFloat64() - 1 },
		"mass0": func() float64 { // 30% exactly 0, as the polarity-free √ΔSI²c
			if rng.Float64() < 0.3 {
				return 0
			}
			return rng.ExpFloat64()
		},
	} {
		fit := fitNullCDF(sample(200_000, draw))
		fresh := sample(200_000, draw)
		var m, v float64
		for _, x := range fresh {
			s := fit.score(x)
			m += s
			v += s * s
		}
		m /= float64(len(fresh))
		sd := math.Sqrt(v/float64(len(fresh)) - m*m)
		// A mass of ties shares its midpoint score, so the scores are not
		// exactly N(0, 1): for 30% at 0 the mean is 0.3·Φ⁻¹(0.15) + φ(Φ⁻¹(0.3))
		// ≈ 0.037 and the SD below 1.
		wantMean := 0.0
		if name == "mass0" {
			wantMean = 0.3*normalQuantile(0.15) + math.Exp(-0.5*math.Pow(normalQuantile(0.3), 2))/math.Sqrt(2*math.Pi)
		}
		if math.Abs(m-wantMean) > 0.01 || (name != "mass0" && math.Abs(sd-1) > 0.02) || sd > 1.02 {
			t.Fatalf("%s: null scores mean %.3f, SD %.3f", name, m, sd)
		}
		last := fit.Knots[cdfKnots-1]
		if a, b := fit.score(last), fit.score(last+1e-9); math.Abs(a-b) > 1e-6 {
			t.Fatalf("%s: tail not continuous: %v, %v", name, a, b)
		}
		if !(fit.score(last+1) > fit.score(last) && fit.score(last+10) > fit.score(last+1)) {
			t.Fatalf("%s: tail not increasing", name)
		}
	}
	mass := fitNullCDF(sample(100_000, func() float64 {
		if rng.Float64() < 0.3 {
			return 0
		}
		return 1 + rng.ExpFloat64()
	}))
	if got, want := mass.score(0), normalQuantile(0.15); math.Abs(got-want) > 0.01 {
		t.Fatalf("tied mass at 0: score %v, want the midpoint Φ⁻¹(0.15) = %v", got, want)
	}
}

// smallNullRaw is a small two-chromosome marker set for null tests.
func smallNullRaw() []BSAstats {
	var raw []BSAstats
	for _, c := range []string{"c1", "c2"} {
		for i := range 1500 {
			d := 30 + i%25
			raw = append(raw, BSAstats{CHROM: c, POS: int64(1000 + 10_000*i), HighDepth: d, LowDepth: 70 - d, Depth: effectiveDepth(d, 70-d), OneBulkP0: 0.5})
		}
	}
	return raw
}

// On fresh null genomes (not those that fitted the model) both components
// are close to N(0, 1), their correlation is close to the model's r, and
// DsiGconsensusZ has unit variance.
func TestDsiGconsensusNullCalibrated(t *testing.T) {
	color.Output = io.Discard
	raw := smallNullRaw()
	for _, mode := range []struct{ mode, ref string }{{"2p2b", "high"}, {"2b", "none"}, {"hphb", ""}} {
		cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: mode.ref,
			NullReps: 60, BrmAlpha: 0.01}
		th, err := genomeNullThresholds(cfg, mode.mode, raw)
		if err != nil {
			t.Fatal(err)
		}
		m := th.Consensus.Model
		sc, _ := newScorer(cfg, mode.mode)
		sim, _ := newNullSim(cfg, mode.mode, raw)
		sim.rho = th.Rho
		order := newMarkerOrder(raw)
		spec, _ := newConsensusSpec(cfg, mode.mode, th.Rho)
		pool := newNullPool(sim, sc, order, 1)
		var s1, s2, c []float64
		err = pool.run(40, 99, false, func(_ int, sm []SmoothedStats, _ float64) error {
			for i := range sm {
				mom := m.Moments[sm[i].CHROM]
				mean, sd := mom.at(sm[i].POS)
				nzd, nzg, cons := m.scores(spec, &sm[i], mean, sd)
				s1, s2, c = append(s1, m.CDF1.score(nzd*nzd)), append(s2, nzg), append(c, cons)
			}
			return nil
		})
		if err != nil {
			t.Fatal(err)
		}
		mv := func(x []float64) (float64, float64) {
			var a, b float64
			for _, v := range x {
				a += v
				b += v * v
			}
			a /= float64(len(x))
			return a, math.Sqrt(b/float64(len(x)) - a*a)
		}
		// Null genomes are a few linked blocks each, so 40 of them pin the
		// moments only roughly.
		for name, x := range map[string][]float64{"ΔSI component": s1, "G component": s2} {
			if mean, sd := mv(x); math.Abs(mean) > 0.25 || math.Abs(sd-1) > 0.25 {
				t.Errorf("%s %s: fresh null mean %.3f SD %.3f", mode.mode, name, mean, sd)
			}
		}
		if r := pearson(s1, s2); math.Abs(r-m.R) > 0.15 || !(m.R > 0 && m.R < 1) {
			t.Errorf("%s: model r = %.3f, fresh null r = %.3f", mode.mode, m.R, r)
		}
		if _, sd := mv(c); math.Abs(sd-1) > 0.25 {
			t.Errorf("%s: DsiGconsensusZ null SD %.3f", mode.mode, sd)
		}
		if !(th.Consensus.P99 > th.Consensus.P95 && th.Consensus.P95 > 1 && th.Consensus.P99 < 10) {
			t.Errorf("%s: DsiGconsensusZ thresholds p95 %v p99 %v", mode.mode, th.Consensus.P95, th.Consensus.P99)
		}
	}
}

// The null's results do not depend on how many replicates run at once.
func TestNullThresholdsIndependentOfWorkers(t *testing.T) {
	color.Output = io.Discard
	raw := smallNullRaw()
	cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high",
		NullReps: 24, BrmAlpha: 0.01}
	run := func(outer int) Thresholds {
		nullOuterOverride = outer
		defer func() { nullOuterOverride = 0 }()
		th, err := genomeNullThresholds(cfg, "2p2b", raw)
		if err != nil {
			t.Fatal(err)
		}
		return th
	}
	a, b := run(1), run(5)
	if !reflect.DeepEqual(a, b) {
		t.Fatalf("thresholds differ between 1 and 5 concurrent replicates:\n%+v\n%+v", a, b)
	}
}

// maxHMMLOD (with the lodUpper pruning) equals the maximum of the full
// per-site LOD, in two-bulk and one-bulk modes.
func TestMaxHMMLODExact(t *testing.T) {
	rng := rand.New(rand.NewSource(5))
	var raw []BSAstats
	for _, c := range []string{"a", "b", "c"} {
		f := 0.5
		for i := range 3000 {
			f = math.Min(0.95, math.Max(0.05, f+0.01*rng.NormFloat64()))
			if c == "b" && i > 1200 && i < 1500 {
				f = 0.85 // a QTL-like stretch
			}
			dh, dl := 20+rng.Intn(40), 20+rng.Intn(40)
			hh, hl := distuv.Binomial{N: float64(dh), P: f, Src: rng}.Rand(), distuv.Binomial{N: float64(dl), P: 1 - f, Src: rng}.Rand()
			raw = append(raw, BSAstats{CHROM: c, POS: int64(1000 + 3000*i), HighBulkH: int(hh), HighBulkL: dh - int(hh),
				LowBulkH: int(hl), LowBulkL: dl - int(hl), HighDepth: dh, LowDepth: dl, OneBulkP0: 0.5})
		}
	}
	for _, mode := range []string{"2p2b", "hphb"} {
		cfg := utils.AnalysisConfig{Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high"}
		setup, err := newHMMSetup(cfg, mode)
		if err != nil {
			t.Fatal(err)
		}
		order := newMarkerOrder(raw)
		cp := append([]BSAstats(nil), raw...)
		computeHMMLOD(setup, cp, order)
		want := 0.0
		for i := range cp {
			want = math.Max(want, cp[i].HMMLOD)
		}
		for _, w := range []int{1, 3} {
			if got := setup.maxHMMLODWorkers(raw, order, w); got != want {
				t.Fatalf("%s, %d workers: max LOD %v, full computation %v", mode, w, got, want)
			}
		}
		if mode == "2p2b" && want < 5 {
			t.Fatalf("test data has no clear peak (max LOD %v)", want)
		}
	}
}

// The null cache round-trips, and its key changes with the data and with
// the settings the null depends on.
func TestNullCache(t *testing.T) {
	color.Output = io.Discard
	raw := smallNullRaw()
	cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high",
		NullReps: 20, BrmAlpha: 0.01}
	th, err := genomeNullThresholds(cfg, "2p2b", raw)
	if err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "stats", "null_cache.gob")
	key := nullCacheKey(cfg, "2p2b", raw)
	if err := saveNullCache(path, key, th); err != nil {
		t.Fatal(err)
	}
	got, ok := loadNullCache(path, key)
	if !ok || !reflect.DeepEqual(got, th) {
		t.Fatalf("cache round trip: ok %v\n%+v\n%+v", ok, got, th)
	}
	if _, ok := loadNullCache(path, key+"x"); ok {
		t.Fatal("cache hit with another key")
	}
	changed := func(f func(c *utils.AnalysisConfig, r []BSAstats)) bool {
		c, r := cfg, append([]BSAstats(nil), raw...)
		f(&c, r)
		return nullCacheKey(c, "2p2b", r) != key
	}
	for name, f := range map[string]func(c *utils.AnalysisConfig, r []BSAstats){
		"null reps":   func(c *utils.AnalysisConfig, _ []BSAstats) { c.NullReps = 21 },
		"bulk size":   func(c *utils.AnalysisConfig, _ []BSAstats) { c.LowBulkSize = 21 },
		"window":      func(c *utils.AnalysisConfig, _ []BSAstats) { c.WindowSize++ },
		"population":  func(c *utils.AnalysisConfig, _ []BSAstats) { c.Population = "F3" },
		"a depth":     func(_ *utils.AnalysisConfig, r []BSAstats) { r[7].HighDepth++ },
		"an SI value": func(_ *utils.AnalysisConfig, r []BSAstats) { r[9].LowSI += 1e-9 },
	} {
		if !changed(f) {
			t.Errorf("cache key unchanged when the %s changes", name)
		}
	}
	if changed(func(c *utils.AnalysisConfig, _ []BSAstats) { c.MinQTLMarkers, c.MergeDistance = 9, 123 }) {
		t.Error("cache key changed with QTL-calling options the null does not use")
	}
}

// Without a null model the NZ columns and DsiGconsensusZ stay zero; with one,
// G moments are interpolated between the positions the null scored.
func TestApplyDsiGConsensus(t *testing.T) {
	cfg := utils.AnalysisConfig{Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high"}
	sm := []SmoothedStats{{CHROM: "c", POS: 150, SmDeltaSI: 0.2, SmGstat: 5, ReadVar0: 1e-4}}
	if err := applyDsiGConsensus(cfg, "2p2b", sm, Thresholds{}); err != nil || sm[0].DsiGconsensusZ != 0 || sm[0].NZDeltaSI != 0 {
		t.Fatalf("no model: %+v, %v", sm[0], err)
	}
	g := &GMoments{Pos: []int64{100, 200}, Mean: []float64{2, 4}, SD: []float64{1, 3}}
	if mean, sd := g.at(150); mean != 3 || sd != 2 {
		t.Fatalf("interpolated moments %v, %v", mean, sd)
	}
	if mean, sd := g.at(50); mean != 2 || sd != 1 {
		t.Fatalf("moments before the first position %v, %v", mean, sd)
	}
	std := fitNullCDF([]float64{-2, -1, 0, 1, 2})
	m := &ConsensusModel{DSIScale: 1, Moments: map[string]*GMoments{"c": g}, CDF1: std, CDF2: std, R: 0.5}
	th := Thresholds{Consensus: ConsensusThresholds{Model: m}}
	if err := applyDsiGConsensus(cfg, "2p2b", sm, th); err != nil {
		t.Fatal(err)
	}
	spec, _ := newConsensusSpec(cfg, "2p2b", 0)
	if want := spec.z(&sm[0]); math.Abs(sm[0].NZDeltaSI-want) > 1e-12 || want <= 0 {
		t.Fatalf("NZDeltaSI %v, want z = %v", sm[0].NZDeltaSI, want)
	}
	if want := std.score((5.0 - 3) / 2); sm[0].NZGstat != want {
		t.Fatalf("NZGstat %v, want %v (G standardised with the interpolated moments)", sm[0].NZGstat, want)
	}
	if want := (std.score(sm[0].NZDeltaSI*sm[0].NZDeltaSI) + sm[0].NZGstat) / math.Sqrt(3); math.Abs(sm[0].DsiGconsensusZ-want) > 1e-12 {
		t.Fatalf("DsiGconsensusZ %v, want (s1 + s2)/√(2 + 2r) = %v", sm[0].DsiGconsensusZ, want)
	}
}
