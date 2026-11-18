# GoBSAseq Biological, Statistical, and Code Review

## Summary

`go test ./...` and `go vet ./...` pass. However, several correctness issues could materially change significance calls and QTL intervals.

## High severity

### 1. Low-bulk BRM mode can produce false block calls

Location: `stats/brm.go`, `calculateBRMBlocksOneBulk`.

When `SmAFDev == 0`, the one-bulk BRM function substitutes `SmHighSI - expectedSI`. Low-bulk analyses do not populate `SmHighSI`, so this becomes `0 - expectedSI`, regardless of the true low-bulk deviation.

**Impact:** Low-bulk-only BRM detection can classify blocks from an incorrect deviation.

**Recommended fix:** Use `SmAFDev` directly for one-bulk deviation detection. Ensure low-bulk raw statistics populate `SmAFDev` consistently.

### 2. Unequal bulk sequencing depths are not preserved

Locations: `stats/stats.go`, `RawStats`; `stats/thresholds.go`, `CalculateThresholds`.

`BSAstats.Depth` stores the smaller of the two bulk depths. Two-bulk simulations then set both `highDepth` and `lowDepth` equal to that value. The `simParams` type already has separate depth fields, but callers do not assign the distinct high and low depths.

**Impact:** Two-bulk experiments with unequal high and low sequencing depths receive incorrect empirical thresholds.

**Recommended fix:** Preserve separate high and low bulk depths throughout `BSAstats`, threshold assignment, cache keys, and simulation code. Avoid using one minimum depth as a proxy for both.

### 3. Composite-Z construction is inconsistent with its Monte Carlo model

Locations: `stats/smoothing.go`, `consolidate`; `stats/thresholds.go`, `simulateNullCompositeZTwoBulk` and `simulateNullCompositeZOneBulk`.

The real pipeline computes `|ZDeltaSI| + ZGstat + ZLOD` or `|ZAFDev| + ZOneBulkG + ZOneBulkLOD`. The comments imply all components are non-negative, but `ZGstat` and `ZLOD` are not guaranteed to be non-negative, and Monte Carlo columns do not implement the same complete transformation.

**Impact:** Negative robust Z-scores can cancel true positive signal, and the empirical Composite-Z threshold may not match the real statistic distribution.

**Recommended fix:** Define and implement one explicit, consistent rule for each component: for example, use absolute values or upper-tail only. Mirror that exact rule in both the real pipeline and Monte Carlo simulation.

## Medium severity

### 4. One-bulk Monte Carlo cache omits population structure

Location: `stats/thresholds.go`, `oneBulkCacheKey`.

The cache key does not include `population`. F2 and RIL can share `p0=0.5` while having different Stage-1 genotype distributions.

**Impact:** Thresholds generated for one population structure may be reused for another.

**Recommended fix:** Include the normalized population string in the one-bulk cache key. Also review whether the two-bulk cache key is complete under all supported modes.

### 5. Backcross population models are incomplete or inconsistent

Locations: `stats/stats.go`, `ExpectedAF`; `stats/thresholds.go`, `samplePoolAF`; `stats/brm.go`, `PopulationVarianceScale`.

`BCnF2H/L` is accepted by `ExpectedAF`, but `samplePoolAF` does not implement a distinct BCnF2 genotype model. Backcross forms are handled with assumptions that differ from the variance-scale implementation.

**Impact:** Backcross null distributions and BRM thresholds may not describe the intended genetics.

**Recommended fix:** Implement and document distinct genotype distributions for BCnH, BCnL, BCnF2H, and BCnF2L, including exact expected allele frequencies and segregation variances. Add independently calculated tests for each mode.

## Lower severity

- Formula tests are not sufficiently independent. Tests should assert output values from known count examples rather than only structural properties.
- Existing tests cover pool-sampling basics and QTL merging, but not formula correctness, cache isolation, unequal bulk depths, or low-bulk BRM behavior.
- Malformed and missing VCF fields need broader tests, especially missing or malformed DP/AD, alternate-allele polarity, and multiallelic records.
- Some comments describe Fisher/Z/Bayes methods more confidently than the current implementations demonstrate. Methods should be documented as approximate or validated against simulated or external reference data before publication.

## Validation performed

- `go test ./...` passes.
- `go vet ./...` passes.
- No repository files were changed during the review.
