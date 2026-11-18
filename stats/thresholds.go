package stats

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"hash/fnv"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"gonum.org/v1/gonum/stat"
)

// ── Threshold structs ────────────────────────────────────────────────────────

type Thresholds struct {
	TwoBulk TwoBulkThresholds
	OneBulk OneBulkThresholds
	Z       ZThresholds
	BRM     BRMThreshold

	// SupportDrop is the LOD drop (smoothed-LOD units) that ends a support
	// interval: SupportZDrop × the null LOD scale (see support.go).
	SupportDrop float64
}

// BRMThreshold is BRM's genome-wide multiplier: a BRM block is called where
// |ΔSI| (or |SI − p0|) ≥ U·σ (see "BRM significance" in brm.go). U is the
// (1 − Alpha) quantile of the per-genome maximum of |ΔSI|/σ in null genomes.
type BRMThreshold struct {
	U     float64
	Alpha float64
}

type TwoBulkThresholds struct {
	HighSIP99, HighSIP95     float64
	HighSIMp99, HighSIMp95   float64
	LowSIP99, LowSIP95       float64
	LowSIMp99, LowSIMp95     float64
	DeltaSIP99, DeltaSIP95   float64
	DeltaSIMp99, DeltaSIMp95 float64

	GstatP99, GstatP95     float64
	ED4P99, ED4P95         float64
	LODP99, LODP95         float64
	BBLogBFP99, BBLogBFP95 float64
}

type OneBulkThresholds struct {
	AFDevP99, AFDevP95                   float64
	AFDevMp99, AFDevMp95                 float64 // lower tail
	OneBulkGstatP99, OneBulkGstatP95     float64
	OneBulkLODP99, OneBulkLODP95         float64
	OneBulkBBLogBFP99, OneBulkBBLogBFP95 float64
}

type ZThresholds struct {
	ZP99, ZP95                   float64
	CompositeZP99, CompositeZP95 float64
	NumStats                     int
}

// deriveSeed deterministically mixes a base seed with a simulation's identity
// (its parameter key), so every simulation gets a distinct but reproducible
// RNG stream no matter which goroutine runs it or in what order.
func deriveSeed(base int64, key string) int64 {
	h := fnv.New64a()
	_, _ = h.Write([]byte(key))
	return base ^ int64(h.Sum64())
}

// countZStats returns k, the number of Z-scores summed by consolidate() into
// CompositeZ.  It MUST mirror the czs slice built in consolidate() exactly;
// any change to the statistic set there requires a matching change here.
//
// Current mapping (see consolidate() in smoothing.go for rationale):
//
//	Two-bulk  → k = 3  (ZDeltaSI, ZGstat, ZLOD)
//	One-bulk  → k = 3  (ZAFDev, ZOneBulkG, ZOneBulkLOD)
//
// Note: HighSI, LowSI, ED4, BBLogBF, and OneBulkBBLogBF are intentionally
// excluded from CompositeZ (they remain available in MaxAbsZ and individual
// Z-score columns).
func countZStats(bsaType string) int {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
	n := 0
	if hasBothBulks {
		n += 3 // ZDeltaSI, ZGstat, ZLOD
	}
	if hasOneBulk {
		n += 3 // ZAFDev, ZOneBulkG, ZOneBulkLOD
	}
	return n
}

// ── Public entry point ───────────────────────────────────────────────────────

// CalculateThresholds derives genome-wide significance thresholds by
// simulating cfg.Rep null genomes (see nullsim.go). They are genome-wide, so
// one Thresholds value applies to every marker (the thresholds TSV still
// lists it per marker for convenience).
//
// Interpretation: a "P99" threshold is exceeded anywhere in the genome by only
// 1% of null genomes (genome-wide false-positive rate 1%); "P95" by 5%. Signed
// tracks (SI, ΔSI, AFDev) split the error rate across both tails (99.5th
// percentile of the per-genome maximum, 0.5th of the minimum).
func CalculateThresholds(
	cfg utils.AnalysisConfig,
	bsaType string,
	raw []BSAstats,
	smoothed []SmoothedStats,
) (Thresholds, error) {
	var none Thresholds

	if len(smoothed) == 0 || len(raw) == 0 {
		return none, fmt.Errorf("no data to calculate thresholds from")
	}
	if cfg.Rep <= 0 {
		return none, fmt.Errorf("--rep must be > 0, got %d", cfg.Rep)
	}
	if cfg.WindowSize <= 0 {
		return none, fmt.Errorf("WindowSize must be > 0, got %d", cfg.WindowSize)
	}
	p0, err := ExpectedAF(cfg.Population)
	if err != nil {
		return none, fmt.Errorf("cannot determine null AF: %w", err)
	}
	design, err := parseCrossDesign(cfg.Population)
	if err != nil {
		return none, err
	}

	hasHigh, hasLow, hasBoth, hasOne := BulkFlags(bsaType)
	highSize, lowSize, oneSize := bulkSizes(cfg, bsaType)
	if !(cfg.BrmAlpha > 0 && cfg.BrmAlpha < 1) {
		return none, fmt.Errorf("--brm-alpha must be between 0 and 1, got %v", cfg.BrmAlpha)
	}
	if (hasHigh && !hasOne && highSize <= 0) || (hasLow && !hasOne && lowSize <= 0) || (hasOne && oneSize <= 0) {
		return none, fmt.Errorf("bulk sizes must be > 0 (got high=%d low=%d); set them with --bulk-sizes", highSize, lowSize)
	}

	sigma := float64(cfg.WindowSize) / 2
	layout := buildNullLayout(raw, smoothed, bsaType, sigma, contigLengths(cfg), design, p0, highSize, lowSize, oneSize)
	layout.brmScale = PopulationVarianceScale(cfg.Population)

	color.Cyan(
		"\n============================ Calculating Thresholds (%s, p₀=%.4f, %d null genomes) ============================\n",
		bsaType, p0, cfg.Rep,
	)
	color.Blue("Null model: %s bulks of %s, %.0f cM per chromosome (assumed), observed marker positions and depths, seed %d",
		strings.ToUpper(cfg.Population), bulkSizeLabel(hasBoth, highSize, lowSize, oneSize), morgansPerChromosome*100, cfg.Seed)

	// Z-score scales from a separate batch of null genomes; the real data and
	// every threshold genome are normalised with the same fixed transform.
	layout.scales, layout.scalesSet = layout.estimateScales(cfg.Rep, cfg.Seed), true
	normaliseToNull(smoothed, layout.scales, hasHigh, hasLow, hasBoth, hasOne)
	if err := writeSmoothedOutput(cfg, bsaType, smoothed); err != nil {
		return none, err
	}
	if err := writeNullScalesTSV(filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.null_scales.tsv", bsaType)), layout); err != nil {
		return none, fmt.Errorf("failed to write null scales TSV: %w", err)
	}

	ext := simulateNullExtremes(layout, cfg.Rep, cfg.Seed)

	col := func(f func(nullExtremes) float64) []float64 {
		xs := make([]float64, len(ext))
		for i, e := range ext {
			xs[i] = f(e)
		}
		sort.Float64s(xs)
		return xs
	}
	q := func(xs []float64, p float64) float64 { return r6(stat.Quantile(p, stat.Empirical, xs, nil)) }
	// upper/lower for signed tracks: both tails share the error rate.
	upper := func(f func(nullExtremes) float64) (p99, p95 float64) {
		xs := col(f)
		return q(xs, 0.995), q(xs, 0.975)
	}
	lower := func(f func(nullExtremes) float64) (p99, p95 float64) {
		xs := col(f)
		return q(xs, 0.005), q(xs, 0.025)
	}
	oneSided := func(f func(nullExtremes) float64) (p99, p95 float64) {
		xs := col(f)
		return q(xs, 0.99), q(xs, 0.95)
	}

	var t Thresholds
	tb := &t.TwoBulk
	if hasHigh {
		tb.HighSIP99, tb.HighSIP95 = upper(func(e nullExtremes) float64 { return e.HighSIMax })
		tb.HighSIMp99, tb.HighSIMp95 = lower(func(e nullExtremes) float64 { return e.HighSIMin })
	}
	if hasLow {
		tb.LowSIP99, tb.LowSIP95 = upper(func(e nullExtremes) float64 { return e.LowSIMax })
		tb.LowSIMp99, tb.LowSIMp95 = lower(func(e nullExtremes) float64 { return e.LowSIMin })
	}
	if hasBoth {
		tb.DeltaSIP99, tb.DeltaSIP95 = upper(func(e nullExtremes) float64 { return e.DeltaSIMax })
		tb.DeltaSIMp99, tb.DeltaSIMp95 = lower(func(e nullExtremes) float64 { return e.DeltaSIMin })
		tb.GstatP99, tb.GstatP95 = oneSided(func(e nullExtremes) float64 { return e.GstatMax })
		tb.ED4P99, tb.ED4P95 = oneSided(func(e nullExtremes) float64 { return e.ED4Max })
		tb.LODP99, tb.LODP95 = oneSided(func(e nullExtremes) float64 { return e.LODMax })
		tb.BBLogBFP99, tb.BBLogBFP95 = oneSided(func(e nullExtremes) float64 { return e.BBLogBFMax })
	}
	if hasOne {
		ob := &t.OneBulk
		ob.AFDevP99, ob.AFDevP95 = upper(func(e nullExtremes) float64 { return e.AFDevMax })
		ob.AFDevMp99, ob.AFDevMp95 = lower(func(e nullExtremes) float64 { return e.AFDevMin })
		ob.OneBulkGstatP99, ob.OneBulkGstatP95 = oneSided(func(e nullExtremes) float64 { return e.OneGMax })
		ob.OneBulkLODP99, ob.OneBulkLODP95 = oneSided(func(e nullExtremes) float64 { return e.OneLODMax })
		ob.OneBulkBBLogBFP99, ob.OneBulkBBLogBFP95 = oneSided(func(e nullExtremes) float64 { return e.OneBBMax })
	}
	czP99, czP95 := oneSided(func(e nullExtremes) float64 { return e.CompositeZMax })
	t.Z = ZThresholds{
		ZP99: 2.576, ZP95: 1.960, // nominal per-stat Z guides for plots only
		CompositeZP99: czP99, CompositeZP95: czP95,
		NumStats: countZStats(bsaType),
	}

	brmZ := col(func(e nullExtremes) float64 { return e.BRMZMax })
	t.BRM = BRMThreshold{U: q(brmZ, 1-cfg.BrmAlpha), Alpha: cfg.BrmAlpha}
	if tail := cfg.BrmAlpha * float64(cfg.Rep); tail < 10 {
		color.Yellow("Note: the BRM threshold (alpha %.3g) rests on %.0f simulated genome(s) in the tail; increase --rep for a stable value.", cfg.BrmAlpha, tail)
	}

	lodTrack := tLOD
	if !hasBoth {
		lodTrack = tLOD1
	}
	t.SupportDrop = SupportZDrop * layout.scales[lodTrack].Scale

	color.Green("✔ Genome-wide thresholds: CompositeZ P99=%.3f P95=%.3f; BRM u=%.3f (genome-wide alpha %.3g); support-interval drop %.3f LOD", czP99, czP95, t.BRM.U, t.BRM.Alpha, t.SupportDrop)

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := writeNullExtremesTSV(filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.null_genome_extremes.tsv", bsaType)), ext, hasBoth, hasOne); err != nil {
		return t, fmt.Errorf("failed to write null extremes TSV: %w", err)
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.thresholds.tsv", bsaType))
	if err := writeThresholdsTSV(outPath, smoothed, t, bsaType); err != nil {
		return t, fmt.Errorf("failed to write thresholds TSV: %w", err)
	}
	return t, nil
}

func bulkSizeLabel(hasBoth bool, high, low, one int) string {
	if hasBoth {
		return fmt.Sprintf("%d (high) and %d (low) individuals", high, low)
	}
	return fmt.Sprintf("%d individuals", one)
}

// contigLengths returns contig lengths (bp) from the VCF header when present.
func contigLengths(cfg utils.AnalysisConfig) map[string]float64 {
	out := map[string]float64{}
	if cfg.Rdr == nil || cfg.Rdr.Header == nil {
		return out
	}
	for _, c := range cfg.Rdr.Header.Contigs {
		var n float64
		if _, err := fmt.Sscan(c["length"], &n); err == nil && n > 0 {
			out[c["ID"]] = n
		}
	}
	return out
}

// writeNullScalesTSV writes the null median and scale used for each Z-score.
func writeNullScalesTSV(path string, L *nullLayout) error {
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		return err
	}
	names := [nTracks]string{
		tHighSI: "HighSI", tLowSI: "LowSI", tDeltaSI: "DeltaSI", tG: "Gstat", tED4: "ED4",
		tLOD: "LOD", tBB: "BBLogBF", tAFDev: "AFDev", tG1: "OneBulkGstat", tLOD1: "OneBulkLOD", tBB1: "OneBulkBBLogBF",
	}
	var b strings.Builder
	b.WriteString("Statistic\tNullMedian\tNullScale\n")
	for _, tr := range L.activeTracks() {
		fmt.Fprintf(&b, "%s\t%.6g\t%.6g\n", names[tr], L.scales[tr].Med, L.scales[tr].Scale)
	}
	return os.WriteFile(path, []byte(b.String()), 0644)
}

// writeNullExtremesTSV writes the per-genome extremes behind the thresholds,
// so calibration can be inspected or re-derived at other error rates.
func writeNullExtremesTSV(path string, ext []nullExtremes, hasBoth, hasOne bool) error {
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	w := bufio.NewWriter(f)
	cols := []string{"Replicate"}
	if hasBoth {
		cols = append(cols, "DeltaSI_max", "DeltaSI_min", "Gstat_max", "ED4_max", "LOD_max", "BBLogBF_max")
	}
	if hasOne {
		cols = append(cols, "AFDev_max", "AFDev_min", "OneBulkGstat_max", "OneBulkLOD_max", "OneBulkBBLogBF_max")
	}
	cols = append(cols, "CompositeZ_max", "BRM_maxAbsZ")
	fmt.Fprintln(w, strings.Join(cols, "\t"))
	f6 := func(v float64) string { return fmt.Sprintf("%.6f", v) }
	for i, e := range ext {
		row := []string{fmt.Sprint(i)}
		if hasBoth {
			row = append(row, f6(e.DeltaSIMax), f6(e.DeltaSIMin), f6(e.GstatMax), f6(e.ED4Max), f6(e.LODMax), f6(e.BBLogBFMax))
		}
		if hasOne {
			row = append(row, f6(e.AFDevMax), f6(e.AFDevMin), f6(e.OneGMax), f6(e.OneLODMax), f6(e.OneBBMax))
		}
		row = append(row, f6(e.CompositeZMax), f6(e.BRMZMax))
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}
	if err := w.Flush(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// ── TSV output ───────────────────────────────────────────────────────────────

// writeThresholdsTSV writes one row per variant with its empirical thresholds.
// Column order within each statistic group is always P99 before P95, and upper
// tail before lower tail, matching the struct field order.
func writeThresholdsTSV(
	path string,
	smoothed []SmoothedStats,
	t Thresholds,
	bsaType string,
) error {
	if err := os.MkdirAll(filepath.Dir(path), 0755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()

	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()

	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	// ── header ──────────────────────────────────────────────────────────────
	header := []string{"CHROM", "POS", "DEPTH"}

	if hasBothBulks {
		header = append(header,
			"HighSI_P99", "HighSI_P95", "HighSI_Mp99", "HighSI_Mp95",
			"LowSI_P99", "LowSI_P95", "LowSI_Mp99", "LowSI_Mp95",
			"DeltaSI_P99", "DeltaSI_P95", "DeltaSI_Mp99", "DeltaSI_Mp95",
			"Gstat_P99", "Gstat_P95",
			"ED4_P99", "ED4_P95",
			"LOD_P99", "LOD_P95", // P99 before P95 — consistent throughout
			"BBLogBF_P99", "BBLogBF_P95",
		)
	}
	if hasOneBulk {
		header = append(header,
			"AFDev_P99", "AFDev_P95", "AFDev_Mp99", "AFDev_Mp95",
			"OneBulkGstat_P99", "OneBulkGstat_P95",
			"OneBulkLOD_P99", "OneBulkLOD_P95",
			"OneBulkBBLogBF_P99", "OneBulkBBLogBF_P95",
		)
	}
	header = append(header,
		"Z_P99", "Z_P95",
		"CompositeZ_P99", "CompositeZ_P95",
		"BRM_U",
	)

	if err := w.Write(header); err != nil {
		return err
	}

	// ── rows ─────────────────────────────────────────────────────────────────
	f4 := func(v float64) string { return fmt.Sprintf("%.4f", v) }

	for _, sm := range smoothed {
		row := []string{
			sm.CHROM,
			fmt.Sprintf("%d", sm.POS),
			fmt.Sprintf("%d", sm.Depth),
		}

		if hasBothBulks {
			tb := t.TwoBulk
			row = append(row,
				f4(tb.HighSIP99), f4(tb.HighSIP95), f4(tb.HighSIMp99), f4(tb.HighSIMp95),
				f4(tb.LowSIP99), f4(tb.LowSIP95), f4(tb.LowSIMp99), f4(tb.LowSIMp95),
				f4(tb.DeltaSIP99), f4(tb.DeltaSIP95), f4(tb.DeltaSIMp99), f4(tb.DeltaSIMp95),
				f4(tb.GstatP99), f4(tb.GstatP95),
				f4(tb.ED4P99), f4(tb.ED4P95),
				f4(tb.LODP99), f4(tb.LODP95),
				f4(tb.BBLogBFP99), f4(tb.BBLogBFP95),
			)
		}
		if hasOneBulk {
			ob := t.OneBulk
			row = append(row,
				f4(ob.AFDevP99), f4(ob.AFDevP95), f4(ob.AFDevMp99), f4(ob.AFDevMp95),
				f4(ob.OneBulkGstatP99), f4(ob.OneBulkGstatP95),
				f4(ob.OneBulkLODP99), f4(ob.OneBulkLODP95),
				f4(ob.OneBulkBBLogBFP99), f4(ob.OneBulkBBLogBFP95),
			)
		}

		z := t.Z
		row = append(row,
			f4(z.ZP99), f4(z.ZP95),
			f4(z.CompositeZP99), f4(z.CompositeZP95),
			f4(t.BRM.U),
		)

		if err := w.Write(row); err != nil {
			return err
		}
	}

	color.Green("✔ Thresholds written to %s\n", path)
	return nil
}

// ── Helper ───────────────────────────────────────────────────────────────────

// r6 rounds a float64 to 6 decimal places.
func r6(v float64) float64 {
	return math.Round(v*1e6) / 1e6
}
