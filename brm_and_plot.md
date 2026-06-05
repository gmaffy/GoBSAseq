Plan: Implement BRM, Fisher test, and plotting (based on ../v1 and literature)

Summary

This plan describes how to implement BRM (block-regression style interval detection), per-SNP Fisher exact tests, and the plotting pipeline in the current v2 codebase by reusing the logic and insights from ../v1 and accepted BSAseq literature (e.g., QTL-seq/ΔSNP-index, Magwene G-statistic, PyBSASeq, BRM-like block thresholds). The goal is: correct, readable, minimally invasive code that reproduces v1 behavior and adds clean hooks for future permutation-based thresholds.

Key references / findings

- Fisher two-sided exact test on 2×2 allele counts (hypergeometric approach) is standard for per-SNP enrichment (PyBSASeq uses this). Two-sided p-values are computed by summing hypergeometric probabilities as equal-or-more-extreme.
- BRM (block thresholding) in v1: compute per-window expected pooled allele frequency (AFP = mean of bulk SIs), clamp AFP to avoid near-fixation zero variance, compute variance scale = (n1+n2)/(popScale * n1 * n2) where popScale = 2^(popLevel), threshold = uAlpha * sqrt(varianceScale * AFP*(1-AFP)). Flag contiguous windows where |DeltaSI| >= threshold, emit blocks using peak DeltaSI and block boundaries.
- Plotting (v1): uses go-echarts to build interactive HTML charts with: raw smoothed values, p95/p99 horizontal lines, valley (negative-side) lines for two-sided thresholds, and shaded BRM blocks.

Implementation steps (high-level)

1) Stats: Fisher exact
   - Add FisherExact2x2(a,b,c,d int) -> float64 in stats package (copy/adapt v1 hypergeom implementation).
   - Add unit tests adapted from ../v1/crefactor/stats tests.
   - Add a convenience wrapper to compute Fisher p-value per raw SNP (use allele depth counts already collected in RawStats). Option: only run for SNPs that pass depth/MAF filters.

2) BRM block detection
   - Add calculateBRMBlocks (two-bulk) and calculateBRMBlocksOneBulk (one-bulk) to stats (or new brm.go). Implementation will closely follow v1 logic:
     - Compute varianceScale from bulk sizes and popLevel.
     - Compute AFP floor clamp to avoid near-zero variance.
     - Use uAlpha = NormalQuantile(1 - alpha/2) where alpha is user-configurable (cfg.BrmAlpha or command flag).
     - Walk smoothed windows, mark contiguous windows with |stat| >= threshold, record peak position, peak stat and threshold.
   - Output BRM blocks TSV with columns: CHROM, START, STOP, PEAK_POS, PEAK_STAT, BRM_THRESHOLD (and optional VALIDATION/notes).

3) Fisher summary / aggregation
   - Decide where Fisher is useful: per-SNP p-values are useful for prefilter/annotation and as an aggregate per-window measure.
   - Implementation options:
     a) Per-SNP Fisher + BH q-values written as a TSV alongside raw stats (fast, independent).
     b) Window-wise Fisher: for each window, sum allele counts across SNPs and compute Fisher-like test or chi-square on pooled counts (less standard). Prefer (a) with possible per-window aggregation later.
   - Add BH correction implementation (v1 has BenjaminiHochberg in confidence.go) to compute q-values.

4) Plotting
   - Copy/port plotting helpers from v1 into a plotting package or stats/plotting.go:
     - Line charts for smoothed values with p95/p99 reference lines and negative-side thresholds where applicable.
     - Composite & robust-Z overlays (max |Z| across stats) and shaded BRM blocks.
   - Keep external dependency: go-echarts (already used in repo). Produce interactive HTML pages and, optionally, a static PNG export approach (documented). Provide a simple CreatePlots(cfg, stats, thresholds, brmBlocks) entrypoint.

5) Integration into run.bsaseq
   - After smoothing and normalization, call in order:
     a) compute Fisher p-values on raw SNPs (write per-SNP TSV)
     b) compute BRM blocks from smoothed windows (two- and one-bulk branches)
     c) compute and write BRM TSVs and thresholds (already created thresholds.go)
     d) invoke plotting (create HTML in <output>/plots or <output>/stats). Make plotting optional via flag.

6) Tests & validation
   - Port unit tests from v1 (crefactor tests) for Fisher, Wilson CI, Benjamini-Hochberg and BRM block detection.
   - Run example datasets and compare outputs (BRM block counts, peak positions, fisher p-values) vs v1 outputs as regression checks.

7) Performance & practicalities
   - Run per-SNP Fisher in parallel across variants; write streamed TSV to avoid memory overload.
   - Pre-filter low-depth/low-MAF SNPs before running Fisher.
   - Provide command-line flags: --brm-alpha, --fisher-alpha, --plot (bool), --plot-out (dir), --plot-static (bool).

Files to add/modify

- stats/fisher.go       (FisherExact2x2 + per-SNP wrapper)
- stats/brm.go          (calculateBRMBlocks and one-bulk variant)
- stats/plotting.go     (chart creation helpers ported from v1)
- run/run.go            (call the new functions and expose CLI flag for plotting)
- brm_and_plot.md       (this plan; added to repo root)
- tests/ (unit tests copied/adapted)

Milestones / Acceptance criteria

- M1 (dev): Add FisherExact2x2 with tests and per-SNP Fisher TSV output.
- M2 (dev): Add BRM block detection with tests; write BRM TSVs and validate against v1 for sample data.
- M3 (dev): Port plotting helpers; produce interactive HTML plots for a test dataset; plots correctly show thresholds and BRM shading.
- M4 (integration): Wire into bsaseq flow; add flags for enabling/disabling Fisher/BRM/plots; run full pipeline on example dataset and verify outputs match v1 behavior.

Notes / future work

- Add permutation bootstrap to compute empirical thresholds for all stats (expensive but more robust). Provide an option to run N permutations in parallel and combine results with BRM.
- Consider exposing multiple BRM alpha values in a single run and writing per-alpha BRM TSVs (v1 had Alphas array usage).
- Consider adding memory/disk-backed streaming for extremely large VCFs (>10M SNPs).

