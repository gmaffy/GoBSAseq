# GoBSAseq Code Review — Correctness & Biological Accuracy

**Reviewer role:** Computational biologist / Go pipeline developer
**Scope:** `filter/filter.go`, `stats/stats.go`, `stats/smoothing.go`, `stats/thresholds.go`, `stats/detect.go`, `stats/brm.go`, `stats/genespace.go`, `run/run.go`, `cmd/root.go`, `plots/plots.go`, `filter/filter_test.go`
**Method:** Independent read of the source (not a synthesis of `copilot.md`/`cursor.md`/`gemini35flash.md`/`geminipro.md`/`maicode_report.md`/`vibe_report.md`), tracing the actual runtime call graph in `run/run.go` end to end, and checking every statistical claim in `README.md` against the code that implements it.

---

## Executive summary

GoBSAseq is a well-engineered piece of Go software wrapped around a mostly-standard BSA-seq statistical toolkit (SNP-index, ΔSI, G-statistic, LOD, beta-binomial Bayes factor, Stouffer's Z), with two genuinely good ideas layered on top: a two-stage (bulk-size → depth) Monte Carlo null model, and population-structure-aware variance scaling (`PopulationVarianceScale`). The concurrency, VCF-edge-case handling, and code documentation are all above the bar for a bioinformatics side project.

However, tracing the actual pipeline (`run/run.go`, steps 2→8) turns up a **statistical calibration bug that sits at the center of the tool's stated purpose**: the significance thresholds used for QTL calling are simulated from *unsmoothed, single-locus* null data, but they are compared against *Gaussian-smoothed* statistics. Smoothing pools information across many neighboring markers and mechanically shrinks variance; a threshold calibrated on the wider, unsmoothed null will be systematically too strict for smoothed data. Compounding this, `FindPeakIntersections` (the function underlying every one of the "individual-stat" and CompositeZ QTL calls) collapses the per-variant, depth-aware thresholds the whole Monte Carlo apparatus was built to produce into **one flat, chromosome-wide average**, discarding the depth-dependence that is the documented reason for running Monte Carlo simulation per depth bin in the first place. Neither of these is cosmetic — together they mean the "primary" QTL-calling method (`*.qtls.tsv`, CompositeZ) is not testing what its own design doc says it tests.

This is not a fatal flaw (the pipeline will not crash or produce nonsense; it will mostly just be conservative and lose power, and in `brm.go` the opposite direction — no depth term at all — creates a different miscalibration), but it is a real correctness problem in the code's most important deliverable, and none of the five prior AI-generated reviews in this repo (`copilot.md`, `cursor.md`, `gemini35flash.md`, `geminipro.md`, `maicode_report.md`, `vibe_report.md`) caught it — they praised the Monte Carlo model as "state of the art" and "excellently calibrated" without checking whether the calibration target (raw stats) matched the tested quantity (smoothed stats). I flag it first because it changes the overall score materially relative to those reports.

**Score: 6.8 / 10** — a genuinely capable engineering effort with a real statistical-calibration gap in its core detection logic, thin test coverage over exactly the parts that need it most, and several smaller biological simplifications that are defensible as v1 choices but are not documented as such.

---

## 1. What is correct and well done

### 1.1 Pipeline architecture and mode dispatch (`run/run.go`)
The 10-mode dispatch (`2p2b`, `2phb`, `2plb`, `hp2b`, `lp2b`, `hphb`, `hplb`, `lphb`, `lplb`, `2b`) in `bsaseqType()` is a clean, exhaustive boolean case-split on which of {high parent, low parent, high bulk, low bulk} are supplied, and it is threaded consistently through filtering, stats, and TSV headers via the repeated `hasHighParent`/`hasBothBulks`/`hasOneBulk` idiom (`BulkFlags` in `smoothing.go`). This correctly answers the user's "1-bulk vs 2-bulk" requirement, and the same switch is duplicated correctly for the BAM-input path (`bsaSeqTypeBam`).

### 1.2 Hard filtering (`filter/filter.go`)
- The SNP/INDEL split (`classifyVariant`) correctly buckets same-length multi-nucleotide substitutions (MNPs) with SNPs rather than silently dropping them — a real, cited fix (line 236-241 comment).
- `LightFilter` deliberately drops FS/SOR/MQRankSum/ReadPosRankSum for pooled bulk samples, which is the right call: those tests assume a single diploid genotype, and a bulk of dozens of individuals is *expected* to show non-Mendelian, skewed ref/alt ratios that these filters would misclassify as strand/mapping artifacts. This is a real, non-obvious improvement over blindly applying GATK best-practice defaults to pooled data, and it is documented in the code with the correct reasoning (lines 249-260).
- `safeRefDepthSample`/`safeAltDepthsSample` (duplicated in `filter.go` and `stats.go`) sidestep a real `vcfgo` panic on malformed `AD` fields by hand-parsing the field and falling back to `RO`/`AO`. This is a sensible, defensive pattern given how much production VCF data has malformed FORMAT fields.
- `parentAllele()` uses a purity threshold (`parentAllelePurity = 0.85`) plus a small-count escape hatch (`nonRef <= 1 && refDepth > nonRef`) rather than requiring perfect homozygosity, which is realistic — inbred lines routinely retain 1-2% residual heterozygosity or index-hopping contamination.
- Multi-allelic sites are handled by testing every real ALT and keeping the first one that satisfies the segregation filter (`BsaSeqTargetAlt`/`realAltIndices`), rather than the common shortcut of discarding multi-allelic records outright.
- The filter pipeline is a genuine three-stage concurrent producer/worker/writer design with sequence numbers used to restore coordinate order after out-of-order worker completion (`pending map[int]variantResult`), which is correct and is actually exercised by `TestHardFilterPreservesInputOrder`.

### 1.3 Raw statistics (`stats/stats.go`)
- `determineHighAllele` inferring polarity from parental *read depth* rather than trusting `GT[0]` is a real robustness improvement — a parent's called genotype is one bit of information, its full allele-depth profile is more.
- `ExpectedAF` and `PopulationVarianceScale` (in `brm.go`) correctly encode Mendelian expectations for F2/F3/F4/RIL and generalize the backcross formula (`p₀ = 1 − 0.5^(n+1)`) for arbitrary BCₙH/BCₙL — this is the correct closed-form expectation for repeated backcrossing to one parent and is a genuine point of biological rigor rare in comparable open-source tools.
- LOD and Bayes-factor implementations are numerically stable (log-space, `math.Lgamma` for the Beta function, clamped probabilities to avoid `log(0)`).

### 1.4 Two-stage Monte Carlo null (`stats/thresholds.go`)
Sampling `Binomial(bulk_size, p₀)` for the realized pool allele count and then `Binomial(depth, realized_af)` for the observed reads is the statistically correct way to decompose variance in a pooled-sequencing experiment into its two real sources (finite number of individuals, finite number of reads). A naive single-stage `Binomial(depth, p₀)` model — which is what several published BSA tools still use — understates variance at high depth because it silently assumes every read comes from a different individual. This part of the design is sound *in isolation*.

### 1.5 Concurrency and engineering hygiene
- `runtime.GOMAXPROCS(0)`-sized worker pools are used consistently in `filter.go`, `stats.go`, and the threshold cache warm-ups.
- `singleflight.Group` correctly prevents redundant simulation work when many goroutines request the same (depth, bulk-size) cache key simultaneously.
- `validateCompositeZStatCount` is a genuinely good defensive pattern: it panics at runtime if `consolidate()` (which builds the CompositeZ inputs) and `countZStats()` (which tells the null simulator how many components to expect) ever diverge, which is exactly the kind of two-function coupling that silently rots in real projects.

---

## 2. Correctness / biological-accuracy problems

### 2.1 (Major) The QTL-calling null model does not match the QTL-calling statistic

Trace the actual call graph in `run/run.go`:

```
RawStats → SmoothAndNormalise → CalculateThresholds → RunBRM
        → DetectIndividualQTLs → DetectCompositeZQTLs → MergeCompositeBRM
```

`SmoothAndNormalise` (`smoothing.go`) applies depth-weighted **Gaussian kernel smoothing** (`gaussianSmooth`, σ = `WindowSize/2`) to every raw statistic before anything downstream ever sees it. `DetectIndividualQTLs` and `DetectCompositeZQTLs` then compare the *smoothed* values (`s.SmGstat`, `s.SmDeltaSI`, `s.CompositeZ`, …) against thresholds from `CalculateThresholds`.

But `CalculateThresholds` computes those thresholds by calling `simulateTwoBulk`/`simulateOneBulk` and `simulateNullCompositeZTwoBulk`/`...OneBulk`, all of which simulate exactly **one** locus's worth of reads per replicate — there is no kernel step, no neighboring markers, nothing that mimics the smoothing the real data went through (`thresholds.go`, e.g. lines 107-206, 431-499). The 99th-percentile envelope produced is therefore the null distribution of a **single unsmoothed SNP's** statistic.

Gaussian smoothing over a window containing *m* effectively-independent markers shrinks variance roughly by a factor related to *m* (exactly, for the depth-weighted Gaussian kernel used here, an effective-sample-size argument applies). Comparing a shrunk statistic to a threshold calibrated on the un-shrunk distribution means:

- In high-marker-density regions, the smoothed null is far tighter than the simulated null, so real signal has to be disproportionately strong before `SmGstat`/`CompositeZ` reaches the (too-wide) threshold — **the pipeline is more conservative than advertised, and that conservatism varies by local marker density**, which is exactly the kind of confound a genome scan should not have.
- In low-marker-density regions (fewer neighbors to average over), the smoothed statistic is closer to the raw one, so detection there is *comparatively* more sensitive than in dense regions — i.e., **statistical power is not comparable across the genome**, which will bias which QTLs are found toward the sparsely-genotyped parts of the assembly.

This is not a hypothetical: `simulateNullCompositeZTwoBulk` documents in its own comment (lines 409-430) that it is meant to "mirror the genome-wide `normalise()` step," which it does for the *robust-Z normalization*, but it never mirrors the *smoothing* step that comes before normalization. The one-stage bypass is real and verifiable by reading the two functions side by side.

**Fix:** either (a) run the null simulation through the same `gaussianSmooth` step using the *actual* marker spacing/density observed on each chromosome (a block/window permutation of the real positions is the standard way to do this and also naturally handles edge effects), or (b) derive an analytic variance-deflation factor for the Gaussian kernel as a function of local marker density and depth, and apply it to the simulated single-locus null before taking quantiles. Option (a) is more defensible biologically and is not expensive given the caching infrastructure already built.

### 2.2 (Major) Per-variant thresholds are collapsed into one chromosome-wide constant

Independent of §2.1, `FindPeakIntersections` in `detect.go` (lines 210-221) does this:

```go
var threshSum float64
for _, t := range thresholds {
    threshSum += threshFn(t)
}
avgThresh := threshSum / float64(len(thresholds))
```

Every one of `twoBulkStatConfigs`/`oneBulkStatConfigs` (`Gstat`, `DeltaSI±`, `ED4`, `LOD`, `BBLogBF`) and the CompositeZ detector (`DetectCompositeZQTLs`) go through `FindPeakIntersections`/`findPeaksWithFallback`, so **every QTL-calling routine in the codebase averages the per-variant, depth-specific Monte Carlo threshold across the entire chromosome before using it**, rather than comparing each smoothed value to the threshold computed for *that variant's own depth*.

This directly contradicts the design intent stated in the README ("Monte Carlo detection provides ... per-variant thresholds" / "MC QTL Detection — QTLs using per-variant MC thresholds") and undoes essentially all of the depth-conditioning work the two-stage Monte Carlo model and the per-depth cache (`warmUpTwoBulkCache`, keyed by exact depth) were built for. A chromosome with a mix of 20× and 150× coverage regions will get one blended threshold that is too strict for the 20× regions and too lax for the 150× regions.

**Fix:** thread the per-variant threshold (already computed and stored in `thresholds[i]`) through the peak-finding boundary-crossing logic directly (`d1 = y1 - threshFn(th[i])`, not a single scalar `avgThresh`), the same way `y1`/`y2` are already indexed per-variant. This is a localized, mechanical fix — the infrastructure to do it correctly already exists in the `Thresholds` slice that's passed in but currently unused per-element.

### 2.3 (Moderate) BRM threshold ignores sequencing depth entirely

`calculateBRMBlocksTwoBulk` (`brm.go`, lines 114-185) computes

```
σ² = (n1 + n2) / (V_scale · n1 · n2) · p(1-p)
```

using only bulk sizes (`n1`, `n2`) and population variance scale — there is no depth term anywhere in this formula, even though the input (`s.SmDeltaSI`) is a statistic estimated from finite sequencing depth. This is the opposite problem from the two-stage Monte Carlo model elsewhere in the codebase, which correctly recognizes that both bulk-size variance *and* depth variance contribute to the noise in an observed SI. BRM's analytic threshold only models the first component. At low depth (the exact "shallow sequencing" scenario BRM would plausibly be used for as a fast analytic alternative to Monte Carlo), this will understate the true variance of `SmDeltaSI` and over-call BRM blocks. The README's claim that BRM "incorporates bulk size" is true but incomplete — it does not incorporate depth, unlike every other threshold system in the pipeline.

**Fix:** add a depth term, e.g. `σ² = (n1+n2)/(V_scale·n1·n2)·p(1-p) + p(1-p)·(1/depth_high + 1/depth_low)/4` (or derive the exact combined-variance formula via the law of total variance over the same two-stage model used elsewhere), so BRM and the Monte Carlo/CompositeZ thresholds are internally consistent about what generates noise in `ΔSI`.

### 2.4 (Moderate) `GStatistic` applies a +0.5 pseudocount unconditionally, not just for zero cells

```go
hba := float64(highBulkAlt) + 0.5
hbr := float64(highBulkRef) + 0.5
lba := float64(lowBulkAlt) + 0.5
lbr := float64(lowBulkRef) + 0.5
```

Standard treatments of the BSA-seq G-statistic (Magwene et al. 2011; Sokal & Rohlf) add a continuity correction only to cells that are exactly zero, to avoid `0·log(0)` — the `if hba > 0 { ... }` guards later in the same function suggest the zero-handling was already intended to be conditional. Here, the +0.5 is applied to *every* cell regardless of whether it's zero, which systematically shrinks the G-statistic at low depth (where +0.5 is a large fraction of the count) and has negligible effect at high depth. Because this same `GStatistic` function is reused inside the Monte Carlo null simulator (`simulateTwoBulk`, `simulateNullCompositeZTwoBulk`), the null and the real statistic are at least computed consistently with each other — but the choice still deviates from the textbook definition without being called out as a deliberate design decision anywhere in the code or README, and it will shrink power specifically at the low end of the depth range where power is already limited.

**Fix:** either document this explicitly as a deliberate small-sample regularization (and justify the magnitude), or restrict the pseudocount to genuinely zero cells as the literature formula does.

### 2.5 (Minor/documentation) README overstates a fixed CompositeZ ≥ 3.0 threshold

The README states plainly: *"Regions with |CompositeZ| ≥ 3.0 (primary method)"* and lists this as the "Primary" detection rule in the QTL-methods table. The actual code (`calculateEmpiricalCompositeZThresholds`) computes a **dynamic, Monte-Carlo-derived** `CompositeZP99`/`CompositeZP95` per run and falls back to fixed z=2.576/1.960 (not 3.0) only if the simulation fails. Neither branch ever uses literally "3.0". This is a real documentation/implementation mismatch that will mislead a user trying to reason about output stringency from the README alone (and it also means §2.1/§2.2 above bite silently — the user has no way to know from the docs that the actual operative threshold is both dynamically simulated *and* averaged across the chromosome).

### 2.6 (Minor, but should be documented) Bulks-only ("2b") mode has no basis for allele polarity

In `determineHighAllele` (`stats.go`), when neither parent is supplied (`highParentIdx < 0 && lowParentIdx < 0`), the function falls through to `return 0` — i.e., the reference allele is *always* treated as the "high" allele. This is a defensible fallback (many BSA tools without parental genotypes compute SNP-index relative to the reference genome and rely on the *sign/enrichment direction* of ΔSI rather than allele identity to interpret QTLs), but it is silent in the code and undocumented in the README's `2b` row, which just says "Bulks only, no parents" without noting that "HighSI"/"LowSI" in that mode are computed relative to the reference allele, not a phenotype-linked allele, and that a user cannot say "the high-bulk allele" means anything about which parent contributed it. This should be one sentence in the README so users don't over-interpret `HighSI`/`LowSI` column semantics in that mode.

### 2.7 (Minor) No multiple-testing / genome-wide error-rate control
Thresholds are per-statistic p95/p99 quantiles applied independently at every marker (or, per §2.2, every chromosome-average). With millions of markers genome-wide, and no FDR/Bonferroni-style correction and no permutation-based genome-wide significance level, the nominal per-marker α does not translate to a controlled genome-wide false-positive rate. This is a common simplification in BSA-seq tools generally (QTLseqr and MutMap have similar issues), so it's not unique to GoBSAseq, but given how much statistical machinery already exists here (Monte Carlo everything else), a permutation-based genome-wide threshold (label-swap the two bulks, rerun the whole smoothing+CompositeZ pipeline, take the max CompositeZ across the genome per permutation, then take the empirical 95th/99th percentile of *that*) would be a natural, and correct, extension, and would also sidestep §2.1/§2.2 by construction since it would run the real smoothing step on real marker spacing.

### 2.8 (Minor) Test coverage does not reach the statistical core
`filter/filter_test.go` contains exactly two tests, both about `filter.go` (multi-allelic ALT selection and coordinate-sorted output preservation). There are **zero automated tests** for `stats.go` (G-statistic, LOD, Bayes factor, `ExpectedAF`, `determineHighAllele`), `smoothing.go` (Gaussian kernel correctness, robust Z, Stouffer combination), `thresholds.go` (Monte Carlo simulators, empirical quantiles), `detect.go` (peak-finding, boundary interpolation, QTL merging), or `brm.go` (variance-scale formulas, block detection). This matters more than usual here because the bugs in §2.1/§2.2 are exactly the kind that a single unit test — "simulate a smoothed constant-signal chromosome, assert the per-variant threshold used in `FindPeakIntersections` equals `thresholds[i]`, not a chromosome average" — would have caught immediately. `go test ./...` passing (as the other reports note) currently only proves the filtering/I-O layer works, not that the statistics or QTL calls are correct.

---

## 3. Scorecard

| Component | Score | Rationale |
|---|---|---|
| CLI / mode dispatch | 9/10 | Clean, exhaustive, correctly threaded through the whole pipeline. |
| Hard filtering (`filter.go`) | 8.5/10 | Thoughtful pooled-data-aware design (`LightFilter`), robust parsing, correct concurrency. Minor: multi-allelic handling still treats each ALT independently rather than jointly modeling segregation across all alleles. |
| Raw statistics (`stats.go`) | 7/10 | Correct core math; unconditional G-stat pseudocount (§2.4) and undocumented reference-allele polarity default in `2b` mode (§2.6) are real but fixable issues. |
| Smoothing (`smoothing.go`) | 7.5/10 | Correct, efficient (binary-search windowing) Gaussian kernel and robust-Z implementation; fixed physical window size rather than marker-density-adaptive bandwidth is a reasonable v1 limitation. |
| Thresholds (`thresholds.go`) | 6/10 | The two-stage Monte Carlo *idea* is sound and well implemented for a single locus, but it is never actually calibrated against the smoothed statistic it's meant to threshold (§2.1) — a genuine miscalibration of the pipeline's central method. |
| QTL detection (`detect.go`) | 5.5/10 | The chromosome-average-threshold bug (§2.2) undermines the depth-aware design of the whole threshold system; peak interpolation and interval-merging logic themselves are implemented correctly. |
| BRM (`brm.go`) | 6.5/10 | Analytically elegant population-variance scaling, but omits sequencing-depth variance entirely (§2.3), making it inconsistent with the rest of the pipeline's noise model. |
| Documentation accuracy | 6/10 | README overstates a fixed CompositeZ≥3.0 rule that the code doesn't implement (§2.5), and doesn't document the `2b`-mode polarity default (§2.6) or that "per-variant" MC thresholds are currently chromosome-averaged (§2.2). |
| Test coverage | 4/10 | Two tests total, both on the least statistically risky module; the modules with real bugs in this review have zero coverage. |
| Engineering (concurrency, error handling) | 9/10 | Genuinely strong: worker pools, `singleflight` dedup, sequence-preserving pipelines, defensive VCF parsing, a real runtime consistency check (`validateCompositeZStatCount`). |

**Overall: 6.8 / 10.**

---

## 4. Roadmap to 10/10

**Statistical correctness (highest priority — these are the difference between "the tool works" and "the tool's numbers mean what they claim to mean"):**
1. Fix the null-vs-statistic mismatch (§2.1): simulate the null through the actual smoothing step (permutation-based, using real marker spacing) rather than through a single unsmoothed locus.
2. Fix per-variant threshold averaging (§2.2): compare each smoothed statistic to its own per-variant threshold in `FindPeakIntersections`, not a chromosome-wide mean.
3. Add a depth term to the BRM variance formula (§2.3) so it's consistent with the two-stage model used everywhere else.
4. Either restrict the G-statistic's +0.5 correction to zero-count cells or explicitly document/justify the always-on pseudocount (§2.4).
5. Add genome-wide error-rate control via permutation (bulk label-swapping), which also structurally fixes §2.1 for free.

**Documentation fidelity:**
6. Correct the README's "CompositeZ ≥ 3.0" claim to describe the actual dynamic/empirical threshold, and note that it is currently a per-chromosome average pending fix #2.
7. Document the `2b`-mode reference-allele polarity convention explicitly (§2.6).

**Test coverage:**
8. Unit tests for `stats.go` (G-stat/LOD/BF against hand-computed values), `smoothing.go` (Gaussian kernel against a closed-form example), `thresholds.go` (quantile sanity checks, cache-key correctness), and — most importantly — a regression test asserting `FindPeakIntersections` uses per-variant, not averaged, thresholds.

**Scalability (secondary, but real for genome-scale data as other reviews correctly note):**
9. Streaming/chunked (per-chromosome) processing instead of holding all variants/smoothed stats/thresholds in memory simultaneously.
10. Depth-binning in the Monte Carlo cache for very deep or highly variable-depth datasets to bound cache size.
11. Plot data-thinning for chromosomes with millions of markers so the HTML outputs stay browser-renderable.

Addressing items 1-5 is what would move this from a well-built pipeline with a real calibration bug in its central method to a tool whose QTL calls can be trusted at face value; items 6-11 are what would make it publication- and production-grade on top of that.
