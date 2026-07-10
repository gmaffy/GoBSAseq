# GoBSAseq technical and biological review

**Reviewer:** Computational biologist and Go pipeline developer  
**Date:** 2026-07-10  
**Scope:** Complete Go pipeline and the existing repository Markdown reviews.

## Verdict

**Overall score: 6.0 / 10.**

GoBSAseq has a promising, well-organised Go implementation: the mode dispatch is clear, VCF handling is unusually defensive, parallel filtering preserves order, and the raw-statistic calculations are transparent. `go test ./...` passes (with `GOCACHE` redirected to a writable temporary directory); however, only the `filter` package has substantive tests.

The score is held down by three issues that affect reported biological discoveries rather than just presentation: (1) finite-pool simulation treats diploid individuals as one Bernoulli allele draw, (2) thresholds are generated for unsmoothed single markers but are applied to smoothed tracks, and (3) final merging reduces all CompositeZ and BRM intervals on each chromosome to one broad QTL. Consequently, the current final QTL workbook should be treated as exploratory, not a calibrated genome-wide inference result.

## What works well

- **Pipeline structure and role selection.** VCF modes cover two-bulk and one-bulk designs, and mode selection rejects designs without a bulk ([run/run.go](run/run.go:87)). The stages are ordered sensibly: filter, raw statistics, smoothing, thresholds, QTL calling, plots, and optional annotation ([run/run.go](run/run.go:232)).
- **VCF robustness.** The filter reads AD/RO/AO defensively, carries the selected ALT allele with a multiallelic record, preserves source ordering despite worker concurrency, and creates indexed filtered VCF output ([filter/filter.go](filter/filter.go:334), [filter/filter.go](filter/filter.go:760)).
- **Useful core statistics.** SNP-index/allele frequency, delta SNP-index, G statistic, likelihood-ratio LOD, and beta-binomial Bayes-factor-like scores are computed from read counts with reasonable numerical guards ([stats/stats.go](stats/stats.go:324), [stats/stats.go](stats/stats.go:445)).
- **Efficient smoothing implementation.** Per-chromosome Gaussian smoothing uses a bounded kernel and binary search; sparse kernels are discarded rather than silently extrapolated ([stats/smoothing.go](stats/smoothing.go:162)).
- **Traceable outputs.** Raw, smoothed, threshold, BRM, Excel, and interactive HTML outputs make it possible to audit a call. The plots correctly keep threshold records aligned with smoothed records before rendering ([plots/plots.go](plots/plots.go:63)).

## Critical issues — correct before biological use

### 1. Finite-pool Monte Carlo uses the wrong biological sampling unit

The two-stage null simulator samples `Binomial(N=bulkSize, p=p0)` and describes `bulkSize` as a count of individuals ([stats/thresholds.go](stats/thresholds.go:114), [stats/thresholds.go](stats/thresholds.go:138)). For a diploid F2 or backcross pool, allele frequencies arise from **2N alleles**, not N Bernoulli draws. Using N inflates the pool-sampling variance approximately twofold for a simple F2 allele-frequency model, making thresholds too permissive and increasing false positives. It is only naturally defensible for fully homozygous lines such as a RIL when `bulkSize` represents independent fixed genotypes.

This is not solved by `PopulationVarianceScale`: that correction appears only in BRM, not in the Monte Carlo thresholds. More fundamentally, a binomial allele draw cannot represent every listed population. F2, F3, advanced selfed generations, RILs, and backcrosses differ in genotype/haplotype variance and linkage structure.

**Required correction:** define `bulk-size` explicitly as individuals; simulate individual genotype dosage (0/1/2) or parental haplotypes under the selected cross, then sample reads from the realised pool allele frequency. Validate the simulator empirically against known F2, BC, and RIL variance formulas. If a simpler model is retained, use 2N for diploid allele pools and clearly limit the supported population assumptions.

### 2. The threshold null is not the statistic being tested

Individual thresholds are simulated for an **unsmoothed one-site** statistic at a single depth ([stats/thresholds.go](stats/thresholds.go:107)), but QTL detection compares them to a Gaussian, depth-weighted average over nearby markers ([stats/smoothing.go](stats/smoothing.go:200)). Smoothing changes both variance and correlation substantially. Assigning the threshold from the focal record's `Depth` does not model the effective depth, marker spacing, or number of contributors in that kernel ([stats/thresholds.go](stats/thresholds.go:861)).

CompositeZ has the same central issue. Its null uses one representative (median) depth and normalises independent simulated loci, rather than simulating the actual chromosome positions, coverages, smoothing, and genome-wide robust-Z operation ([stats/thresholds.go](stats/thresholds.go:614)). Its empirical quantile is therefore not a genome-wide threshold for the plotted CompositeZ track.

**Required correction:** generate null data across the actual marker grid (or phenotype/bulk-label permutations where appropriate), apply the identical smoothing and normalisation code, and retain the maximum statistic per genome for each replicate. Use the resulting max-statistic percentile for family-wise error control. This resolves smoothing calibration and multiple testing together.

### 3. Per-site thresholding is discarded during individual QTL detection

`CalculateThresholds` constructs one threshold object per record, but `FindPeakIntersections` first averages them over an entire chromosome and uses only that constant thereafter ([stats/detect.go](stats/detect.go:216), [stats/detect.go](stats/detect.go:260)). This defeats depth-aware thresholding. Low-depth sites can be called too readily while high-depth sites can be penalised, and boundary interpolation is also wrong where depth changes.

**Required correction:** use `t1 := threshFn(thresholds[i])` and `t2 := threshFn(thresholds[i+1])`, including in the first-marker and end-marker cases. Add a regression test with two different depths and a signal that crosses only one of the two thresholds.

### 4. Final QTL merging collapses independent loci on the same chromosome

`bestCompositePeakByChrom` retains one CompositeZ peak per chromosome, and `summarizeBRMByChrom` expands all BRM blocks on that chromosome into one span ([stats/detect.go](stats/detect.go:796), [stats/detect.go](stats/detect.go:814)). `MergeCompositeBRM` then returns at most one final record for that chromosome ([stats/detect.go](stats/detect.go:845)). Two genuine QTLs on chromosome 1 will therefore be reported as one interval potentially spanning most of the chromosome.

**Required correction:** interval-join CompositeZ peaks and BRM blocks by reciprocal/any overlap (with an explicit user-configurable gap), retaining every disconnected component. Do not summarise by chromosome.

## Important biological and statistical limitations

### One-bulk BSA-seq is underpowered and needs stronger framing

The one-bulk modes compare one selected pool to a Mendelian background expectation (`p0`) ([stats/stats.go](stats/stats.go:336)). This can be useful, but it lacks the matched opposing bulk that protects a two-bulk analysis against locus-specific mapping bias, reference bias, population structure, and stochastic founder distortion. It should be labelled as an exploratory/single-tail design, with stronger requirements for parental informativeness and a recommendation to validate candidates in an independent bulk or individual genotyping panel.

### BRM thresholds omit sequencing and smoothing uncertainty

BRM uses only bulk-size/population variance and a normal quantile ([stats/brm.go](stats/brm.go:104)); sequencing depth, overdispersion, and kernel smoothing are absent. At shallow coverage this can be anti-conservative; at deep coverage it is inconsistent with the two-stage model. BRM should either use the same calibrated simulation framework as the main test or be presented as a heuristic support track rather than merged as equivalent evidence.

### The composite is not an independent Stouffer test

Delta-SNP-index, G, and LOD are all transformations of the same two-by-two read-count table. Limiting the composite to three metrics reduces redundancy, but it does not make the inputs independent ([stats/smoothing.go](stats/smoothing.go:451)). An empirical null can account for correlation only if it also mirrors the real spatial smoothing and marker properties. Until then, report CompositeZ as a ranking score, not a p-value or independently calibrated Z statistic.

### Hard filtering can retain known failed or unannotated calls

The hard filter does not require VCF `FILTER=PASS`; it decides from selected INFO values only ([filter/filter.go](filter/filter.go:246)). Missing QD, MQ, FS, SOR, and rank-sum fields generally pass because checks execute only when a value exists. This is permissive for caller-specific VCFs but unsafe as a silent default.

The light-filter rationale also says strand and rank-sum artefact tests are invalid because a bulk is not 50:50 ([filter/filter.go](filter/filter.go:249)). Strand-bias tests assess whether the *strand distribution* differs by allele; they do not require equal reference and alternative counts. Skewed biological allele frequency does not inherently justify ignoring severe strand bias.

**Required correction:** make `FILTER=PASS` the default requirement (with an opt-out), report missing annotation fields, expose caller-specific profiles, and retain or carefully relax artefact filters based on calibrated validation rather than expected segregation ratio alone.

### Multiallelic handling can distort frequencies

The pipeline selects one ALT at a multiallelic site, permits up to 20% support for other ALT alleles, then calculates the focal frequency as focal ALT versus reference while excluding those other reads ([filter/filter.go](filter/filter.go:125), [stats/stats.go](stats/stats.go:280)). This creates a denominator smaller than the actual coverage and can inflate a focal SNP index. Decompose multiallelic records into valid biallelic records upstream, or exclude loci with appreciable non-target ALT depth.

### Parent polarity is convenient but not proof of linkage phase

High-allele direction is inferred from the predominant parental read allele and otherwise defaults to reference ([stats/stats.go](stats/stats.go:265)). This is reasonable for clean, inbred parents, but partial-parent and bulks-only modes have no robust phase validation. In `2b`, “high allele” is simply the reference allele, so the sign is reference-relative rather than phenotype-parent-relative. Document this prominently and add per-site parent-purity/phase diagnostics.

### Filtering and modelling do not address overdispersion or replicates

Read counts are treated as binomial after a finite-pool draw. PCR duplicates, unequal DNA contributions, allele-specific mapping, local CNV, and bulk composition commonly create beta-binomial overdispersion. The tool accepts one sample per bulk and has no replicate-aware model. Add optional technical/biological bulk replicates, beta-binomial or Dirichlet-multinomial dispersion estimation, mapping-quality/base-quality/read-position exclusions, duplicate handling, and reference-bias checks.

## Engineering, usability, and reproducibility findings

| Area | Assessment | Needed improvement |
|---|---|---|
| Tests | `go test ./...` passes, but only `filter/filter_test.go` exercises substantive behaviour. | Add hand-calculated unit tests for G/LOD/BF, population expectations, smoothing, threshold assignment, peak boundaries, multi-QTL merge, and simulations. Add synthetic VCF integration tests with planted QTLs and null data. |
| Window semantics | Documentation calls `--window-size` Gaussian sigma, while code uses `WindowSize / 2` as sigma ([stats/smoothing.go](stats/smoothing.go:158)). | Pick one definition, rename/help-text it precisely, and record the effective sigma and cutoff in output metadata. |
| Step size | `StepSize` is parsed and stored but is not used in the analysis. | Remove it, or implement regular genomic-grid/window summaries and document their relationship to marker-centred smoothing. |
| QTL width | Single points and very short runs can be emitted as QTLs. | Require a minimum physical width and/or effective-marker count, with a transparent relaxation option. |
| Scalability | Raw variants, smoothed values, and thresholds all reside in memory. | Stream chromosome-by-chromosome, persist intermediate tables, and downsample plot payloads for dense genomes. |
| Reproducibility | Simulation seeds are implicit and not saved with results. | Add `--seed`, print it, and write a machine-readable run manifest with version, command, VCF header/sample mapping, thresholds, and parameters. |
| Output naming | README descriptions include outputs/method claims not fully aligned with code (for example fixed CompositeZ cut-offs and MC-QTL wording). | Generate documentation from tested behaviour and include a method/assumptions report with every run. |

## Prioritised roadmap to 10/10

1. **Repair inference first:** replace the finite-pool simulator with population-appropriate genotype/haplotype sampling and validate it analytically and by simulation.
2. **Calibrate the actual statistic:** use whole-genome permutation/null simulations through the same smoothing and normalisation path; control genome-wide error by max statistic.
3. **Fix calls and merging:** compare individual tracks with local thresholds and preserve all non-overlapping QTL intervals in the final merge.
4. **Harden VCF QC:** require/pass through FILTER state, handle missing INFO explicitly, biallelicise multiallelic sites, and add allele-bias/overdispersion safeguards.
5. **Make modes scientifically explicit:** distinguish diploid/F2, backcross, RIL, and one-bulk assumptions; document reference-relative polarity in bulks-only mode.
6. **Add a validation suite:** simulated null datasets should meet nominal false-positive rates, and planted QTL datasets should demonstrate power/localisation across coverage, bulk size, and cross type.
7. **Improve production readiness:** streaming per chromosome, deterministic run manifests, meaningful `--step-size`, and plot thinning.

## Final assessment

This is a strong engineering foundation with an especially good start on robust VCF processing and traceable outputs. It is not yet biologically calibrated enough to be called production-grade BSA-seq software. After the first three roadmap items are implemented and validated on simulated plus real benchmark crosses, the project could credibly reach **9–10/10**: the architecture is already suitable for that evolution.
