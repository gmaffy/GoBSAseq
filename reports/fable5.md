# GoBSAseq Code Review — Fable 5

**Reviewer:** Claude Fable 5 (independent read of the source tree)
**Date:** 2026-07-15
**Scope reviewed:** `main.go`, `cmd/root.go`, `run/run.go`, `filter/filter.go`, `stats/{stats,smoothing,thresholds,detect,brm,genespace}.go`, `utils/{config,utils}.go`, plus the `filter` and `stats` test files. `plots/plots.go` skimmed at the interface level.
**Method:** Traced the runtime call graph in `run/run.go` end to end, checked each statistical routine against the biology it claims to model, and verified build/vet/test status locally. This is an independent read, not a synthesis of the other reports in `reports/`.

Build, `go vet`, and `go test ./...` all pass cleanly (two packages have tests; the rest have none).

---

## Executive summary

GoBSAseq is a capably engineered BSA-seq pipeline. The Go is idiomatic and concurrent where it matters, the VCF-handling is defensively written against real-world malformed data, and the statistical toolkit (SNP-index, ΔSI, G-statistic, ED⁴, LOD, beta-binomial Bayes factor, Stouffer's composite Z, BRM) is broad and mostly correct. Two design choices stand out as genuinely good: the **two-stage Monte-Carlo null** (sample pool allele frequency, *then* sample reads) and **population-structure-aware variance scaling**.

The most important problem is a **calibration mismatch in the core QTL-calling path**: significance thresholds are simulated from *single-locus, unsmoothed* null draws, but they are applied to *Gaussian-smoothed* statistics. Smoothing pools neighbouring markers and shrinks variance, so a threshold calibrated on unsmoothed data is systematically too strict — the tool loses power at true QTLs and the reported significance levels (p99/p95) do not mean what they say. This affects both the individual-stat and CompositeZ tracks, though CompositeZ is partially self-correcting because it re-normalises.

Beyond that, the read-alignment input path is advertised but not implemented (it silently succeeds without doing anything), output directories no longer carry timestamps (re-runs overwrite prior results), and test coverage is thin over exactly the statistically hardest code.

**Overall: a solid, above-average bioinformatics codebase with one real correctness gap in its headline deliverable and a couple of rough edges around unfinished features.**

---

## What is done well

**Architecture and mode dispatch.** The ten analysis modes (`2p2b`, `2phb`, `2plb`, `hp2b`, `lp2b`, `hphb`, `hplb`, `lphb`, `lplb`, `2b`) are an exhaustive boolean case-split over which of {high parent, low parent, high bulk, low bulk} are present (`run/run.go:87`). The same mode string threads consistently through filtering, stats, thresholds, TSV headers, and Excel output via the `BulkFlags` idiom (`stats/smoothing.go:625`), and the BAM-input path mirrors it (`bsaSeqTypeBam`).

**Defensive VCF parsing.** `safeRefDepthSample` / `safeAltDepthsSample` (`stats/stats.go:59`, `filter/filter.go:51`) hand-parse `AD`/`RO`/`AO` FORMAT fields to sidestep a known `vcfgo` panic on malformed `AD`, and the reader loop tolerates `.`-as-missing parse errors by clearing and continuing (`filter/filter.go:668`). This is the kind of hardening real VCF data demands.

**LightFilter reasoning.** Dropping FS/SOR/MQRankSum/ReadPosRankSum for pooled bulks (`filter/filter.go:255`) is biologically correct and well-argued in the code comment: those strand/position tests assume a single diploid genotype, and a bulk of many individuals is *expected* to show skewed ref/alt ratios that GATK's single-sample filters would misread as artifacts.

**Two-stage null model.** `simulateTwoBulk` / `simulateOneBulk` (`stats/thresholds.go:160`) draw the realised pool allele frequency from a population-appropriate genotype sampler (`samplePoolAF`), then draw observed reads `Binomial(depth, realized_af)`. This correctly decomposes pooled-sequencing variance into its two real sources (finite individuals, finite reads); a naive single-stage `Binomial(depth, p0)` understates variance at high depth.

**Population genetics.** `ExpectedAF` (`stats/stats.go:484`) and `PopulationVarianceScale` (`stats/brm.go:68`) encode Mendelian expectations for F2/F3/F4/RIL and generalise the backcross closed forms for arbitrary BCₙ. This level of population-structure rigor is uncommon in comparable open-source tools.

**Concurrency and hygiene.** Worker pools are sized by `GOMAXPROCS`; `singleflight.Group` deduplicates simulation work per (depth, bulk-size) key (`stats/thresholds.go:255`); the filter pipeline is a genuine reader→workers→writer design that restores coordinate order after out-of-order completion (`filter/filter.go:768`), and this is actually exercised by `TestHardFilterPreservesInputOrder`. `validateCompositeZStatCount` (`stats/smoothing.go:38`) panics early if `consolidate()` and `countZStats()` ever disagree on *k* — a good guard against silent drift between two coupled functions.

**Numerical care.** LOD and Bayes-factor code works in log space, uses `math.Lgamma` for the Beta function, and clamps probabilities away from 0/1 (`stats/stats.go:561`).

---

## Correctness issues

### 1. (Major) Thresholds are calibrated on unsmoothed data but applied to smoothed data

The pipeline order in `bsaseq()` is:

```
RawStats → SmoothAndNormalise → CalculateThresholds → RunBRM
         → DetectIndividualQTLs → DetectCompositeZQTLs → MergeCompositeBRM
```

`SmoothAndNormalise` applies depth-weighted **Gaussian kernel smoothing** (σ = `WindowSize/2`) to every statistic (`stats/smoothing.go:162`). `DetectIndividualQTLs` then compares the *smoothed* values `s.SmGstat`, `s.SmDeltaSI`, `s.SmLOD`, … (`stats/detect.go:393`) against the thresholds from `CalculateThresholds`.

But those thresholds come from `simulateTwoBulk`, which computes each statistic from a **single simulated locus** — one `Binomial(depth, realized_af)` read draw — and takes the 99.5/95th percentile of that per-locus distribution (`stats/thresholds.go:201-238`). No smoothing is applied to the simulated null.

Smoothing a statistic across a kernel of many markers dramatically reduces its variance (roughly by a factor related to the effective number of markers in the window). Comparing a low-variance smoothed value against a threshold derived from the high-variance unsmoothed null means the threshold is **systematically too high**. The practical consequences:

- The pipeline is conservative and loses power — real QTLs whose smoothed peak sits below the unsmoothed p99 are missed.
- The p99/p95 labels are not calibrated: the false-positive rate at "p99" is not ~1%.

The CompositeZ track has the same structural issue but is *partially* self-correcting: both the real CompositeZ (`consolidate`, `stats/smoothing.go:491`) and the simulated null (`simulateNullCompositeZTwoBulk`, `stats/thresholds.go:454`) robust-Z-normalise before combining, so each is expressed on its own dispersion scale. That cancels the first-order variance shrinkage, but not the change in tail shape — smoothing induces autocorrelation and pulls the smoothed distribution toward Gaussian with lighter tails than the unsmoothed null, so a residual conservative bias remains.

The correct fix is to smooth the simulated null the same way the real data is smoothed (e.g. simulate a chromosome-length null track, run it through `gaussianSmooth`, and take percentiles of the smoothed simulated statistic), or equivalently to derive thresholds from a permutation/circular-shift null of the already-smoothed genome-wide track. At minimum this mismatch should be documented as a known limitation.

Note: the *active* `FindPeakIntersections` (`stats/detect.go:210`) does correctly use per-variant, depth-dependent thresholds `threshFn(thresholds[i])` — the older flat chromosome-average version is commented out above it (`stats/detect.go:71`). So the depth-dependence the Monte-Carlo apparatus produces *is* preserved into detection; the problem is purely the smoothed-vs-unsmoothed scale mismatch.

### 2. (Major) The `reads` input path is advertised but silently does nothing

`getRunType` recognises a `reads` run type, and `README.md` documents `--hp-reads`/`--lp-reads`/`--hb-reads`/`--lb-reads`. But:

```go
if runType == "reads" {
    fmt.Printf("Running BSAseq with read-based analysis\n")
    return nil
}
```
(`run/run.go:475`)

It prints a message and returns `nil` — success — without aligning anything or producing output. A user who supplies FASTQ inputs gets a clean exit and no results, with no error to indicate the feature isn't wired up. Either implement the alignment path (the BAM path at `run/run.go:480` shows the intended shape) or return an explicit "not yet implemented" error so the run fails loudly.

### 3. (Moderate) Output directory no longer timestamps — re-runs overwrite

`CreateResultsDir` (`utils/utils.go:28`) has its timestamping logic commented out and now returns a fixed `goBSAseqResults` path. Every run writes into the same directory, so a second run silently overwrites the first's `stats/`, `qtls/`, and `plots/`. If this is intentional (idempotent output location) it's defensible, but the dead timestamp code and the `<timestamp>` reference in the README's quick-start (`ls results/goBSAseqResults/<timestamp>/`) suggest it's an unfinished revert. Reconcile the code, the comment, and the README.

### 4. (Minor) Two detection tracks use inconsistent null models

QTL calls come from two independent nulls that are never reconciled: the Monte-Carlo empirical percentiles (individual stats + CompositeZ) and BRM's *analytic* variance threshold `uAlpha·√variance` with a law-of-total-variance depth term (`stats/brm.go:157`). BRM includes a depth term but no simulation; the MC path simulates but (per issue 1) doesn't smooth. When `MergeCompositeBRM` unions the two (`stats/detect.go:886`), a QTL's support can come from tracks calibrated on different assumptions. Worth a note in the docs about how to interpret `Source = "CompositeZ"` vs `"BRM"` vs `"CompositeZ+BRM"`.

### 5. (Minor) Interactive prompts block non-interactive runs in one case

Sample selection is skipped when names are supplied, which is good. But when gene-space parameters are partially missing, `Run` unconditionally prompts `Continue without gene space analysis? [y/N]` via `fmt.Scan` (`run/run.go:444`). In a scripted/pipeline context with no stdin this will error or hang. Consider a `--yes`/`--no-gene-space` flag to make the run fully non-interactive.

---

## Code-quality notes (non-blocking)

- **Dead code.** `minI64`/`maxI64` (`stats/detect.go:997`) are unused. `OneParentIdx`/`OneParentName` in `AnalysisConfig` (`utils/config.go:60`) are never read. The large commented-out `FindPeakIntersections` block (`stats/detect.go:71-208`) should be deleted — git history preserves it.
- **Duplication.** `safeRefDepthSample`/`safeAltDepthsSample` exist in both `filter` and `stats` verbatim. Hoist to `utils` alongside `GetFloat`.
- **Dead assignments.** `_ = highTotal; _ = lowTotal` in `RawStats` (`stats/stats.go:332`) compute values only to discard them.
- **`cmd/root.go` indentation.** A block of flag reads and `AnalysisConfig` fields (lines 214-217, 382-387) is tab-misaligned relative to its neighbours — cosmetic but `gofmt` would catch it; worth confirming `gofmt -l` is clean in CI.
- **Config-file vs flag paths diverge.** The config-file branch in `cmd/root.go` sets gene-space params but not `cds`/`protein`, while the flag branch sets all of them. Minor inconsistency worth aligning.

---

## Test coverage

Two test files exist (`filter/filter_test.go`, `stats/regression_test.go`) and they cover real invariants: ALT selection flowing into raw stats, coordinate-sorted output, caller-FAILED rejection, local-threshold peak detection, QTL merge/gap behaviour, and the diploid/RIL pool sampler. That's a good start on the highest-risk plumbing.

The gap is the **statistically hardest code**, which is entirely untested:
- `simulateTwoBulk`/`simulateOneBulk` and the CompositeZ null (no test that thresholds have the intended tail probability against a known-null dataset — which is exactly what would have surfaced issue 1).
- `gaussianSmooth` (kernel weights, cutoff, min-contributor sentinel filtering).
- `calculateBRMBlocksTwoBulk`/`OneBulk` and the variance formula.
- `robustZScores` MAD→SD fallback path.

A property test asserting that at a simulated true null, ~1% of smoothed loci exceed the p99 threshold would be the single most valuable test to add.

---

## Recommendations, in priority order

1. **Fix the smoothed-vs-unsmoothed threshold calibration** (issue 1) — smooth the simulated null, or switch to a permutation null on the smoothed track. This is the one change that affects scientific validity.
2. **Make the `reads` path either work or fail loudly** (issue 2).
3. **Add a null-calibration test** proving the p99/p95 false-positive rates hold for smoothed data.
4. **Resolve the output-directory timestamp question** and align the README (issue 3).
5. **Clean up** dead code, duplicated depth parsers, and `gofmt`; add a non-interactive escape for the gene-space prompt.
6. **Document** the two-null-model interpretation for merged QTLs.

The engineering foundation here is strong; the fixes above are about closing the gap between what the statistical machinery is designed to test and what it actually tests.
