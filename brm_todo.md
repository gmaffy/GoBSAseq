# BRM to-do

Goal: biological and code correctness and accuracy at 9/10 or better, with clean, fast code.
Current rating: biology 5.5, accuracy 5, code 6, cleanliness 5.5 (audit of 2026-11-21).

This plan comes from auditing `stats/brm.go` against the original BRM (Block Regression Mapping; Huang et al. 2020, Bioinformatics; github.com/huanglikun/BRM). The original's notes come from reading its source, not running it.

## 1. Biology and accuracy

### 1.1 Genome-wide error control (main gap) — DONE
- **Problem:** the threshold is u(1−α/2)·σ at every marker with α = 0.05 and no multiple-testing correction. In calibration, 22–24 of 100 no-QTL genomes get a block.
- **Original BRM:** Bonferroni over (chromosomes + genetic length / a cM unit that depends on generation). This is an analytic approximation that ignores smoothing and depth.
- **Change:** test a standardised statistic, z = signal / σ_smoothed (1.2), and take the threshold from the genome-wide null (`stats/hmmnull.go`): the (1 − α) quantile of max |z| over the simulated genomes, using `extremeQuantile`. `--brm-alpha` then becomes a genome-wide error rate.
- **Acceptance:** in `TestHMMLODCalibration`, about α·100 of 100 null genomes have a block (assert ≤ 2·α·100 + 3 at α = 0.05), and the QTL is still found.

### 1.2 Null SD of the smoothed signal, not of a single site — DONE
- **Problem:** σ is the variance of one site at that marker's own depth, but the value tested is averaged over the kernel. Smoothing removes most read noise but not pool noise, which is shared along the chromosome. σ is therefore wrong and jumps with each marker's depth, which breaks blocks into fragments.
- **Change:** compute σ²_smoothed = pool variance + Σ_m w_m²·read_m / (Σ_m w_m)². Here w_m are the kernel weights the smoother already computes, and read_m is the beta-binomial read variance, p(1−p)(1 + (d−1)ρ)/d for each bulk, with the estimated ρ. Accumulate Σw²·read in the same kernel pass (one extra field), so it costs almost nothing.
- **Alternative (simpler, exact for the model):** use the per-marker SD of the smoothed signal across the null genomes. This needs the per-marker outputs of the null, which are currently subsampled (`every`). Prefer the analytic form; use the null SD to check it.
- **Acceptance:** the null z has SD 1 ± 0.05 at markers in sparse, dense, shallow and deep regions; no fragmentation in a test with a sharp depth change.

### 1.3 Null frequency under H₀ — DONE
- **Problem:** signed modes use af = the mean of the observed smoothed SI. Under the null it should be p₀ (from the population). In one-bulk modes it should be the bias-corrected P0 for each polarity (`OneBulkP0`).
- **Change:** use p₀ (two-bulk) and the site's `OneBulkP0` (one-bulk). Remove the [0.05, 0.95] clamp, which is unnecessary with p₀.

### 1.4 Effective pool size
- **Problem:** both codebases assume every plant contributes equal DNA. Unequal contributions make the real pool smaller, so the null pool variance, and every threshold, is too low.
- **Change:** estimate a pool-variance inflation factor from the genome-wide spread of smoothed ΔSI (or AFDev) against the pool model, using a robust (MAD-based) estimate so QTLs don't inflate it. Use it in the null simulation (it then affects every statistic, not only BRM) and in σ_smoothed. Report it in the log beside ρ.
- **Acceptance:** with simulated unequal contributions (Dirichlet weights), calibration stays near α; with equal contributions the factor is ≈ 1.

### 1.5 Peaks and intervals instead of raw runs
- **Problem:** blocks are runs of markers above the threshold, so a broad QTL with a dip becomes several blocks, and no position uncertainty is reported.
- **Change** (from the original BRM, corrected):
  - Call peaks as local maxima of |z| above the threshold.
  - Suppress any peak inside a stronger peak's interval.
  - Report each peak's position and a support interval, the region where the signal stays within a fixed number of null SDs of the peak value.
  - Use a two-sided 95% bound (1.96), not the original's one-sided 1.65 per side, and state it.
- Keep the block extent (the region above the threshold) as a separate column.
- **Acceptance:** on simulated QTLs, the 95% interval contains the true position in about 95% of replicates (test over at least 100 replicates).

### 1.6 Polarity-free mode (`2b` without `--ref-parent`) — DONE
- **Keep:** √max(smoothed bias-corrected ΔSI², 0), as now.
- **Change:** its null is not normal (it is |ΔAF|-like and one-sided), so standardise it and calibrate its threshold from the same genome-wide null as 1.1, upper tail only.

### 1.7 Population variance factors
- **Keep** the exact population variance scale (`utils.Population.VarScale`: F2 2, F3 4/3, F4 8/7, BCn, BCnF2). Don't adopt the original's 2^T, which treats F3 like F2.
- **Change:** unknown populations must be an error, not a silent 1.0. `ExpectedAF` already rejects them; make the function return an error too.

## 2. Code correctness

- **Faulty one-bulk fallback:** `if s.SmAFDev != 0 { … } else { SmHighSI − expectedSI }` is wrong in low-bulk-only modes, where SmHighSI is 0 and the deviation becomes −p₀, always significant. Always use `SmAFDev`.
- **Block order changes between runs:** blocks are collected by iterating over a map, so the TSV order differs from run to run. Iterate chromosomes in sorted order.
- **Write errors ignored:** TSV `Fprintf` and `Close` errors are not checked. Use a `bufio.Writer` and return the flush and close errors, as `writeThresholdsTSV` does.
- **Merging:** `mergeBRMBlocks` merges across the merge gap by sign only. With peak calling (1.5) this becomes peak suppression; remove the ad-hoc merge or limit it to intervals of the same sign that overlap.
- **Interval ends:** use the same convention as the other detectors: interpolate where the curve crosses the threshold, as `FindPeakIntersections` does, rather than taking midpoints between markers.
- **Inputs:** with a missing bulk size or α outside (0, 1), return an error instead of silently returning nil.

## 3. Cleanliness and speed

- **One scanner:** replace `calculateBRMBlocksTwoBulk` and `calculateBRMBlocksOneBulk` (near-duplicates, each with its own `emitBlock`) with a single scanner over z, given a signal function and a σ function.
- **Reuse:** use the shared peak/interval code with the other detectors (`stats/detect.go`), so BRM, the individual statistics and CompositeZ call regions the same way.
- **Use gonum:** delete `inverseNormalCDF` (a hand-rolled Acklam approximation) and use `distuv.UnitNormal.Quantile`.
- ~~**Move functions:** move `PopulationVarianceScale` to `utils/population.go`.~~ Done: `utils.ParsePopulation`.
- **Speed and memory:**
  - Don't copy SmoothedStats for each chromosome. Group indices with `markerOrder` (already built), scan in place, and fill z in a reused buffer.
  - σ_smoothed comes from the existing kernel pass (one extra accumulated field), and the null threshold comes from replicates that are already simulated (one more extreme per replicate), so BRM adds almost no run time.
- **Naming:** BRM = Block Regression Mapping. Fix the README's "Bayesian regression model".

## 4. Tests

- **Calibration** (`stats/calibration_test.go`): assert the genome-wide BRM error rate (1.1) for signed `2p2b` (ρ = 0 and 0.01), one-bulk and polarity-free `2b`, and that the QTL is found.
- **σ_smoothed vs the null SD** (1.2): in regions of different density and depth.
- **Interval coverage** (1.5).
- **Unequal pool contributions** (1.4).
- **Regression tests:**
  - the low-bulk-only fallback is gone;
  - block order is stable across runs;
  - a write error returns an error (write to `/dev/full`);
  - unknown populations are rejected.
- **Unit tests:** scanner edge cases (signal above the threshold at the first or last marker, a single-marker peak, NaN threshold).

## 5. Docs

- **README:**
  - BRM = Block Regression Mapping (Huang et al. 2020);
  - what is tested (standardised smoothed ΔSI / AFDev; the polarity-free form);
  - `--brm-alpha` is a genome-wide error rate;
  - the peak and interval output columns;
  - how it differs from the original BRM (exact population variances, read noise and overdispersion, calibration from the null rather than Bonferroni).
- **`smoothing_algorithms.md`:** add a BRM section with the new rating after the changes.

## Order of work

1. Code fixes in section 2: small and independent.
2. σ_smoothed (1.2) and null frequency (1.3), then genome-wide calibration (1.1) and the polarity-free threshold (1.6).
3. One scanner with peaks and intervals (1.5, section 3).
4. Effective pool size (1.4): it changes every threshold, so recalibrate all statistics afterwards.
5. Tests and docs.

## Progress

- **2026-11-21 code fixes 1, 2, 3, 5, 6 (section 2):** done.
  - The one-bulk signal is always `SmAFDev`.
  - Blocks come out in sorted, stable order.
  - Write errors are returned.
  - Block ends are interpolated threshold crossings.
  - Bad inputs return errors.
  - Both scanners share `scanBRMBlocks`.
- **1.2 and 1.3:** done.
  - The Gaussian kernel pass accumulates Σ w²·r·dⁱ, so `ReadVar0 + ρ·ReadVar1` is the exact read-sampling variance of the smoothed ΔSI / AFDev. It accounts for bin aggregation and the local linear blend at chromosome ends.
  - BRM's σ² = pool variance at p₀ (one-bulk: the site's P0) + that read variance, with ρ from the null (`Thresholds.Rho`).
  - Empirical/predicted variance: 1.003 in the interior, 1.03 at chromosome ends, 1.04 in sparse stretches (`TestReadVarMatchesSmoothedVariance`).
  - The per-marker rate of |z| ≥ 1.96 is now 3.9–5.1% (target 5%); the pool term is a single position's variance, a slight upper bound.
  - **Consequence:** because α is still applied per marker, null genomes with a block rose from about 23% to about 50%. 1.1 (genome-wide calibration) is next.
- **1.1 and 1.6:** done.
  - BRM tests z = |signal|/σ (`brmModel`), with u taken from the genome-wide null: the (1 − α) quantile of each null genome's maximum |z| (`extremeQuantile`), stored as `Thresholds.BRMU`. `--brm-alpha` is now a genome-wide error rate, validated in `CalculateThresholds`.
  - The polarity-free signal is calibrated by the same null, one-sided.
  - Calibration at α = 0.05: 2, 4 and 2 of 100 null genomes with a block (2p2b at ρ = 0 and 0.01, 2b without `--ref-parent`); u ≈ 2.9–3.0; the QTL is found in all cases. Asserted in `TestHMMLODCalibration`.
  - The two scanners are now one model (code item 9), and `inverseNormalCDF` is deleted (item 8).
- **Merging (section 2, `mergeBRMBlocks`): kept as is, by decision (2026-11-21).** Blocks with peaks of the same sign within `--merge-distance` (default window/2) are joined without checking the signal between them. Known trade-off: two nearby QTLs can be reported as one block.
- **Population parsing (code item 7):** done.
  - `utils.ParsePopulation` parses the name once into P0, VarScale and MapExpansion, and errors on unknown names. It replaces three separate parsers (`ExpectedAF`, `PopulationVarianceScale` with its silent 1.0 fallback, and `mapExpansion`).
  - `--population` is validated at startup.
  - Results are unchanged for every valid population (`utils/population_test.go`).
- **Per-chromosome copies (code item 10):** done.
  - `byChromosome` (`stats/chroms.go`) uses each chromosome's markers in place when they are one contiguous, position-sorted run of the smoothed output (the usual case); it copies and sorts only otherwise.
  - Used by BRM and by both detection functions, which also gives them a stable chromosome order (detection iterated a map before).
  - Output is byte-identical.
- **Code nits (2026-11-21):** done.
  - Removed the unused `oneSided` field.
  - Removed the hidden single-site variance fallback; `RunBRM` now errors if a marker lacks the smoothed read variance.
  - `mergeBRMBlocks` no longer re-sorts the scanner's ordered blocks.
  - The TSV has a `PEAK_Z` column (signal/σ at the peak).
