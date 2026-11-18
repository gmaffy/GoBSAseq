package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

const minVariantsInKernel = 5

const kernelCutoffSigmas = 3.0

// compositeZStatCount is the number of Z-scores that consolidate() adds to the
// czs slice for each bsaType mode.  It must equal what countZStats() in
// thresholds.go returns; validateCompositeZStatCount() checks this at runtime.
//
// Two-bulk: ZDeltaSI + ZGstat + ZLOD = 3
// One-bulk: ZAFDev + ZOneBulkG + ZOneBulkLOD = 3
const compositeZStatCountPerMode = 3

// validateCompositeZStatCount panics if the number of statistics that consolidate()
// would place into CompositeZ for the given bsaType does not match the value
// returned by countZStats().  This catches any future divergence between the two
// functions at the earliest possible point in the pipeline rather than silently
// producing wrong CompositeZ values and miscalibrated thresholds.
//
// countZStats is passed as a function argument to avoid an import cycle
// (smoothing.go is in package stats; thresholds.go is also in package stats,
// so direct reference is fine — but accepting it as a func keeps the coupling
// explicit and testable).
func validateCompositeZStatCount(bsaType string, countZStats func(string) int) {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	expected := 0
	if hasBothBulks {
		expected += compositeZStatCountPerMode
	}
	if hasOneBulk {
		expected += compositeZStatCountPerMode
	}

	got := countZStats(bsaType)
	if got != expected {
		panic(fmt.Sprintf(
			"CompositeZ stat-count mismatch for bsaType %q: "+
				"consolidate() uses %d statistics but countZStats() returns %d. "+
				"Update countZStats() in thresholds.go to match consolidate() in smoothing.go.",
			bsaType, expected, got,
		))
	}
}

type SmoothedStats struct {
	CHROM string
	POS   int64
	REF   string
	ALT   string
	Depth int

	// ── smoothed raw values ─────────────────────────────────────────────
	// Two-bulk
	SmHighSI  float64
	SmLowSI   float64
	SmDeltaSI float64
	SmGstat   float64
	SmED      float64
	SmLOD     float64
	SmBBLogBF float64
	// One-bulk
	OneBulkP0        float64 // carries the per-site prior (not smoothed; taken from last contributor)
	SmAFDev          float64
	SmOneBulkG       float64
	SmOneBulkLOD     float64
	SmOneBulkBBLogBF float64

	// ── robust Z-scores of smoothed values ──────────────────────────────
	ZHighSI         float64
	ZLowSI          float64
	ZDeltaSI        float64
	ZGstat          float64
	ZED             float64
	ZLOD            float64
	ZBBLogBF        float64
	ZAFDev          float64
	ZOneBulkG       float64
	ZOneBulkLOD     float64
	ZOneBulkBBLogBF float64

	// ── consolidated scores ──────────────────────────────────────────────
	// Stouffer's Z: (|ZDeltaSI| + ZGstat + ZLOD) / √k — always ≥ 0.
	// Direction is handled by the individual-stat tracks; CompositeZ
	// measures signal strength only, enabling a single upper-tail threshold.
	CompositeZ float64
	// max|Z| across all valid statistics (unsigned)
	MaxAbsZ float64
}

func SmoothAndNormalise(cfg utils.AnalysisConfig, bsaType string, rawStats []BSAstats) ([]SmoothedStats, error) {

	// Fix 6: Verify that consolidate() and countZStats() agree on k before doing
	// any work.  Panics immediately if they diverge (e.g. after editing one but
	// not the other), preventing silently wrong CompositeZ values and thresholds.
	validateCompositeZStatCount(bsaType, countZStats)

	if len(rawStats) == 0 {
		return nil, fmt.Errorf("no raw stats supplied to SmoothAndNormalise")
	}
	if cfg.WindowSize <= 0 {
		return nil, fmt.Errorf("WindowSize must be > 0, got %d", cfg.WindowSize)
	}

	bw := float64(cfg.WindowSize) / 2.0

	color.Cyan("\n============================ Smoothing & Normalising (%s) ============================\n\n", bsaType)
	color.Cyan("Gaussian kernel σ = %.0f bp  (WindowSize=%d / 2)\n\n", bw, cfg.WindowSize)

	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	// ── 1. Gaussian kernel smoothing per chromosome ──────────────────────
	smoothed, err := gaussianSmooth(rawStats, bw, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk)
	if err != nil {
		return nil, err
	}

	// Z-scores, CompositeZ and MaxAbsZ need the null scales from the
	// null-genome simulation, so CalculateThresholds fills them (via
	// normaliseToNull) and writes the smoothed TSV.
	color.Green("Smoothed %d variants\n", len(smoothed))
	return smoothed, nil
}

// writeSmoothedOutput writes GoBSAseq.<type>.smoothed_and_normalised.tsv.
func writeSmoothedOutput(cfg utils.AnalysisConfig, bsaType string, smoothed []SmoothedStats) error {
	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return err
	}
	outPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.smoothed_and_normalised.tsv", bsaType))
	if err := writeSmoothedTSV(outPath, smoothed, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk); err != nil {
		return err
	}
	color.Green("Smoothed & normalised stats written to %s (%d variants)\n", outPath, len(smoothed))
	return nil
}

// ─── smoothing ─────────────────────────────────────────────────────────────────

// gaussianSmooth applies depth-weighted Gaussian kernel smoothing to each
// numeric statistic independently, grouped by chromosome.
// bw is the Gaussian σ in base-pairs (= cfg.WindowSize / 2).
// Variants further than kernelCutoffSigmas × bw bp from the query position
// are excluded, reducing per-chromosome complexity from O(n²) to O(n × m).

func gaussianSmooth(raw []BSAstats, bw float64, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) ([]SmoothedStats, error) {

	cutoff := kernelCutoffSigmas * bw // hard distance cutoff in bp

	// Sort a working index by CHROM then POS.
	idx := make([]int, len(raw))
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(a, b int) bool {
		sa, sb := raw[idx[a]], raw[idx[b]]
		if sa.CHROM != sb.CHROM {
			return sa.CHROM < sb.CHROM
		}
		return sa.POS < sb.POS
	})

	out := make([]SmoothedStats, len(raw))

	i := 0
	for i < len(idx) {
		// Find the full extent of the current chromosome.
		chrom := raw[idx[i]].CHROM
		j := i
		for j < len(idx) && raw[idx[j]].CHROM == chrom {
			j++
		}
		chromIdx := idx[i:j]
		n := len(chromIdx)

		// Pre-extract position and depth for this chromosome.
		pos := make([]float64, n)
		depth := make([]float64, n)
		for k, ci := range chromIdx {
			pos[k] = float64(raw[ci].POS)
			depth[k] = math.Max(float64(raw[ci].Depth), 1)
		}

		// Collect (raw values, setter) pairs for every active statistic.
		type statPair struct {
			vals   []float64
			setter func(sm *SmoothedStats, v float64)
		}
		var pairs []statPair

		if hasHighBulk {
			hsi := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.HighSI })
			pairs = append(pairs, statPair{hsi, func(sm *SmoothedStats, v float64) { sm.SmHighSI = v }})
		}
		if hasLowBulk {
			lsi := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.LowSI })
			pairs = append(pairs, statPair{lsi, func(sm *SmoothedStats, v float64) { sm.SmLowSI = v }})
		}
		if hasBothBulks {
			dsi := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.DeltaSI })
			gst := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.Gstat })
			ed := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.ED })
			lod := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.LOD })
			bb := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.BBLogBF })
			pairs = append(pairs,
				statPair{dsi, func(sm *SmoothedStats, v float64) { sm.SmDeltaSI = v }},
				statPair{gst, func(sm *SmoothedStats, v float64) { sm.SmGstat = v }},
				statPair{ed, func(sm *SmoothedStats, v float64) { sm.SmED = v }},
				statPair{lod, func(sm *SmoothedStats, v float64) { sm.SmLOD = v }},
				statPair{bb, func(sm *SmoothedStats, v float64) { sm.SmBBLogBF = v }},
			)
		}
		if hasOneBulk {
			afd := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.OneBulkAFDev })
			og := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.OneBulkGstat })
			ol := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.OneBulkLOD })
			ob := extractFloat(raw, chromIdx, func(s BSAstats) float64 { return s.OneBulkBBLogBF })
			pairs = append(pairs,
				statPair{afd, func(sm *SmoothedStats, v float64) { sm.SmAFDev = v }},
				statPair{og, func(sm *SmoothedStats, v float64) { sm.SmOneBulkG = v }},
				statPair{ol, func(sm *SmoothedStats, v float64) { sm.SmOneBulkLOD = v }},
				statPair{ob, func(sm *SmoothedStats, v float64) { sm.SmOneBulkBBLogBF = v }},
			)
		}

		// Smooth every statistic in one pass: the Gaussian × depth weight of
		// each (query, neighbour) pair is computed once and applied to all
		// statistics, and query positions are split across all cores. Each
		// statistic accumulates in the same order as a per-statistic pass would,
		// so results are bit-for-bit identical to smoothing them one at a time.
		vals := make([][]float64, len(pairs))
		for p, sp := range pairs {
			vals[p] = sp.vals
		}
		smoothedVals, nContrib := smoothAll(pos, depth, vals, bw, cutoff)

		// Assign to output; skip positions with insufficient kernel coverage.
		for k, ci := range chromIdx {
			if nContrib != nil && nContrib[k] < minVariantsInKernel {
				// Mark as invalid by leaving POS = 0; filtered out in caller.
				// We still need identity fields for the skip to be detectable.
				out[ci].CHROM = raw[ci].CHROM
				out[ci].POS = -1 // sentinel: too sparse
				continue
			}
			s := &out[ci]
			s.CHROM = raw[ci].CHROM
			s.POS = raw[ci].POS
			s.REF = raw[ci].REF
			s.ALT = raw[ci].ALT
			s.Depth = raw[ci].Depth
			if hasOneBulk {
				s.OneBulkP0 = raw[ci].OneBulkP0
			}
			for p, sp := range pairs {
				sp.setter(s, smoothedVals[p][k])
			}
		}

		i = j
	}

	// Filter out sentinel (sparse) positions.
	valid := out[:0]
	for _, s := range out {
		if s.POS >= 0 {
			valid = append(valid, s)
		}
	}

	color.Cyan("Smoothed %d variants (%d removed — fewer than %d kernel contributors)\n",
		len(valid), len(out)-len(valid), minVariantsInKernel)

	return valid, nil
}

// smoothAll computes depth-weighted Gaussian smoothed values of every series
// in vals at every position in pos (sorted). It returns the smoothed series
// and, per position, the number of contributing neighbours within cutoff.
func smoothAll(pos, depth []float64, vals [][]float64, bw, cutoff float64) ([][]float64, []int) {
	n := len(pos)
	nv := len(vals)
	out := make([][]float64, nv)
	for p := range out {
		out[p] = make([]float64, n)
	}
	nContrib := make([]int, n)

	workers := runtime.GOMAXPROCS(0)
	chunk := (n + workers - 1) / workers
	if chunk < 256 {
		chunk = 256
	}
	var wg sync.WaitGroup
	for start := 0; start < n; start += chunk {
		end := min(start+chunk, n)
		wg.Add(1)
		go func(start, end int) {
			defer wg.Done()
			wv := make([]float64, nv)
			for k := start; k < end; k++ {
				pK := pos[k]
				lo := sort.Search(n, func(m int) bool { return pos[m] >= pK-cutoff })
				hi := sort.Search(n, func(m int) bool { return pos[m] > pK+cutoff })

				for p := range wv {
					wv[p] = 0
				}
				wSum := 0.0
				for m := lo; m < hi; m++ {
					d := (pos[m] - pK) / bw
					w := depth[m] * math.Exp(-0.5*d*d) // depth × Gaussian
					if w <= 0 {
						continue
					}
					wSum += w
					for p := range wv {
						wv[p] += w * vals[p][m]
					}
					nContrib[k]++
				}
				if wSum > 0 {
					for p := range wv {
						out[p][k] = wv[p] / wSum
					}
				}
			}
		}(start, end)
	}
	wg.Wait()
	return out, nContrib
}

// extractFloat builds a float64 slice for the variants referenced by chromIdx.
func extractFloat(raw []BSAstats, chromIdx []int, fn func(BSAstats) float64) []float64 {
	vals := make([]float64, len(chromIdx))
	for k, ci := range chromIdx {
		vals[k] = fn(raw[ci])
	}
	return vals
}

// ─── normalisation ─────────────────────────────────────────────────────────────

// normaliseToNull fills the Z-score fields with z = (x − median) / scale,
// where median and scale (1.4826 × MAD) come from simulated null genomes
// (NullScales, see nullsim.go), then fills CompositeZ and MaxAbsZ. The same
// fixed transform is applied to every simulated genome, so CompositeZ and
// its threshold are on the same scale.
func normaliseToNull(sm []SmoothedStats, sc NullScales, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) {
	type field struct {
		track int
		get   func(*SmoothedStats) float64
		set   func(*SmoothedStats, float64)
	}
	var fields []field
	if hasHighBulk {
		fields = append(fields,
			field{tHighSI, func(s *SmoothedStats) float64 { return s.SmHighSI }, func(s *SmoothedStats, v float64) { s.ZHighSI = v }},
		)
	}
	if hasLowBulk {
		fields = append(fields,
			field{tLowSI, func(s *SmoothedStats) float64 { return s.SmLowSI }, func(s *SmoothedStats, v float64) { s.ZLowSI = v }},
		)
	}
	if hasBothBulks {
		fields = append(fields,
			field{tDeltaSI, func(s *SmoothedStats) float64 { return s.SmDeltaSI }, func(s *SmoothedStats, v float64) { s.ZDeltaSI = v }},
			field{tG, func(s *SmoothedStats) float64 { return s.SmGstat }, func(s *SmoothedStats, v float64) { s.ZGstat = v }},
			field{tED4, func(s *SmoothedStats) float64 { return s.SmED }, func(s *SmoothedStats, v float64) { s.ZED = v }},
			field{tLOD, func(s *SmoothedStats) float64 { return s.SmLOD }, func(s *SmoothedStats, v float64) { s.ZLOD = v }},
			field{tBB, func(s *SmoothedStats) float64 { return s.SmBBLogBF }, func(s *SmoothedStats, v float64) { s.ZBBLogBF = v }},
		)
	}
	if hasOneBulk {
		fields = append(fields,
			field{tAFDev, func(s *SmoothedStats) float64 { return s.SmAFDev }, func(s *SmoothedStats, v float64) { s.ZAFDev = v }},
			field{tG1, func(s *SmoothedStats) float64 { return s.SmOneBulkG }, func(s *SmoothedStats, v float64) { s.ZOneBulkG = v }},
			field{tLOD1, func(s *SmoothedStats) float64 { return s.SmOneBulkLOD }, func(s *SmoothedStats, v float64) { s.ZOneBulkLOD = v }},
			field{tBB1, func(s *SmoothedStats) float64 { return s.SmOneBulkBBLogBF }, func(s *SmoothedStats, v float64) { s.ZOneBulkBBLogBF = v }},
		)
	}
	for _, f := range fields {
		for i := range sm {
			f.set(&sm[i], sc[f.track].z(f.get(&sm[i])))
		}
	}
	consolidate(sm, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk)
}

// ─── consolidation ─────────────────────────────────────────────────────────────

// consolidate fills CompositeZ (Stouffer's Z) and MaxAbsZ for every variant.
//
// # Statistic selection
//
// CompositeZ combines only a minimal, low-redundancy set of Z-scores:
//
//	Two-bulk: |ZDeltaSI| + ZGstat + ZLOD  (k = 3)
//	One-bulk: |ZAFDev|  + ZOneBulkG + ZOneBulkLOD  (k = 3)
//
// |ZDeltaSI| / |ZAFDev| are used rather than the signed values because
// ZGstat and ZLOD are always ≥ 0 (they measure magnitude of allele frequency
// separation, not direction).  Using the signed ΔSI with unsigned statistics
// produces a right-skewed null under H0, requiring a two-tailed threshold.
// Taking |ZDeltaSI| makes all three components non-negative and the null
// symmetric around zero, so a single upper-tail threshold is both correct
// and sufficient.  Direction is fully captured by the individual-stat tracks
// (ZDeltaSI / ZAFDev retain their signs there).
//
// Excluded and why:
//   - ZHighSI / ZLowSI — subsumed by ZDeltaSI (which is their difference);
//     including them counts the same contrast up to three times.
//   - ZED — a monotone function of |ΔSI|; provides no independent information
//     and inflates the composite for every locus, not just QTLs.
//   - ZBBLogBF / ZOneBulkBBLogBF — strongly correlated with ZGstat and ZLOD
//     (all three derive from the same 2×2 allele-count table); omitting them
//     reduces the spurious null-variance inflation from correlated components
//     without meaningful loss of power at true QTLs.
//
// Stouffer's method assumes (approximately) independent inputs.  Restricting k
// to three statistics with lower pairwise correlation makes the √k normalisation
// more valid and improves calibration of the empirical null.
//
// # MaxAbsZ
//
// MaxAbsZ is computed over ALL available Z-scores (not just the three above) and
// is retained as a secondary diagnostic column for exploratory use.
//
// IMPORTANT: countZStats() in thresholds.go must return the same k = 3 value for
// each bsaType mode so that the CompositeZ null simulation uses the correct
// denominator.  Any change here must be mirrored there.
func consolidate(sm []SmoothedStats, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) {
	for i := range sm {
		s := &sm[i]

		// ── CompositeZ: minimal non-redundant subset ─────────────────────────
		var czs []float64
		if hasBothBulks {
			// |ZDeltaSI| folds both tails into a single positive signal-strength
			// score, matching the always-positive ZGstat and ZLOD components so
			// the Stouffer sum and its null distribution are symmetric.
			czs = append(czs, math.Abs(s.ZDeltaSI), s.ZGstat, s.ZLOD)
		}
		if hasOneBulk {
			czs = append(czs, math.Abs(s.ZAFDev), s.ZOneBulkG, s.ZOneBulkLOD)
		}

		// ── MaxAbsZ: full envelope across all available Z-scores ─────────────
		var allZs []float64
		if hasHighBulk {
			allZs = append(allZs, s.ZHighSI)
		}
		if hasLowBulk {
			allZs = append(allZs, s.ZLowSI)
		}
		if hasBothBulks {
			allZs = append(allZs, s.ZDeltaSI, s.ZGstat, s.ZED, s.ZLOD, s.ZBBLogBF)
		}
		if hasOneBulk {
			allZs = append(allZs, s.ZAFDev, s.ZOneBulkG, s.ZOneBulkLOD, s.ZOneBulkBBLogBF)
		}

		if len(czs) == 0 {
			continue
		}

		s.CompositeZ = stoufferSum(czs)

		maxAbs := 0.0
		for _, z := range allZs {
			if a := math.Abs(z); a > maxAbs {
				maxAbs = a
			}
		}
		s.MaxAbsZ = math.Round(maxAbs*1e6) / 1e6
	}
}

// compositeZ is the CompositeZ of one marker from its signed-difference
// Z-score (ZDeltaSI or ZAFDev) and its two magnitude Z-scores (G and LOD).
// The null simulation (nullsim.go) calls this so the simulated statistic is
// computed by exactly the same formula as the observed one.
func compositeZ(zSigned, zG, zLOD float64) float64 {
	return stoufferSum([]float64{math.Abs(zSigned), zG, zLOD})
}

// stoufferSum returns Σz/√k rounded to 6 decimals.
func stoufferSum(czs []float64) float64 {
	sum := 0.0
	for _, z := range czs {
		sum += z
	}
	return math.Round((sum/math.Sqrt(float64(len(czs))))*1e6) / 1e6
}

// ─── TSV output ────────────────────────────────────────────────────────────────

func writeSmoothedTSV(filename string, sm []SmoothedStats, hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	w := bufio.NewWriter(f)
	defer w.Flush()

	// Build header.
	header := []string{"CHROM", "POS", "REF", "ALT", "Depth"}
	if hasHighBulk {
		header = append(header, "Sm_HighSI", "Z_HighSI")
	}
	if hasLowBulk {
		header = append(header, "Sm_LowSI", "Z_LowSI")
	}
	if hasBothBulks {
		header = append(header,
			"Sm_DeltaSI", "Z_DeltaSI",
			"Sm_Gstat", "Z_Gstat",
			"Sm_ED4", "Z_ED4",
			"Sm_LOD", "Z_LOD",
			"Sm_BBLogBF", "Z_BBLogBF",
		)
	}
	if hasOneBulk {
		header = append(header,
			"P0",
			"Sm_AFDev", "Z_AFDev",
			"Sm_OneBulkGstat", "Z_OneBulkGstat",
			"Sm_OneBulkLOD", "Z_OneBulkLOD",
			"Sm_OneBulkBBLogBF", "Z_OneBulkBBLogBF",
		)
	}
	header = append(header, "CompositeZ", "MaxAbsZ")
	fmt.Fprintln(w, strings.Join(header, "\t"))

	// Write rows.
	f6 := func(v float64) string { return fmt.Sprintf("%.6f", v) }
	for _, s := range sm {
		row := []string{
			s.CHROM,
			fmt.Sprintf("%d", s.POS),
			s.REF,
			s.ALT,
			fmt.Sprintf("%d", s.Depth),
		}
		if hasHighBulk {
			row = append(row, f6(s.SmHighSI), f6(s.ZHighSI))
		}
		if hasLowBulk {
			row = append(row, f6(s.SmLowSI), f6(s.ZLowSI))
		}
		if hasBothBulks {
			row = append(row,
				f6(s.SmDeltaSI), f6(s.ZDeltaSI),
				f6(s.SmGstat), f6(s.ZGstat),
				f6(s.SmED), f6(s.ZED),
				f6(s.SmLOD), f6(s.ZLOD),
				f6(s.SmBBLogBF), f6(s.ZBBLogBF),
			)
		}
		if hasOneBulk {
			row = append(row,
				f6(s.OneBulkP0),
				f6(s.SmAFDev), f6(s.ZAFDev),
				f6(s.SmOneBulkG), f6(s.ZOneBulkG),
				f6(s.SmOneBulkLOD), f6(s.ZOneBulkLOD),
				f6(s.SmOneBulkBBLogBF), f6(s.ZOneBulkBBLogBF),
			)
		}
		row = append(row, f6(s.CompositeZ), f6(s.MaxAbsZ))
		fmt.Fprintln(w, strings.Join(row, "\t"))
	}
	return nil
}

// ─── helpers ───────────────────────────────────────────────────────────────────

func BulkFlags(bsaType string) (hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) {
	hasHighBulk = strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLowBulk = strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	hasBothBulks = hasHighBulk && hasLowBulk
	hasOneBulk = hasHighBulk != hasLowBulk
	return
}
