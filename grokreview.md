# GoBSAseq review

Review of biological, statistical, and code accuracy. The core population model is largely right, but the primary QTL call (CompositeZ) is not calibrated to the statistic that is actually tested. Several other tests are internally consistent only because the same formula is used on both sides, while the number that is reported is not the estimator used in the BSA literature.

## What holds up

The diploid pool sampler for F2, F3, F4, and RIL is the right biological model. An F2 individual is AA, Aa, or aa at 1:2:1, an RIL is fixed, and the pool allele frequency is (allele copies) / (2N). Expected null frequencies are also right: 0.5 for selfing populations, 0.75 / 0.25 for BC1 to the high / low parent, and selfing after a backcross does not change the expected frequency.

`PopulationVarianceScale` matches those genotype frequencies. For a pool of n diploids the allele-frequency variance is p(1−p)/(scale·n), and the implemented scales (F2 = 2, F3 = 4/3, F4 = 8/7, RIL = 1, BC1 = 3, BC1F2 = 1.5) agree with the dosage variances. The two-bulk BRM variance also applies the law of total variance correctly: pool sampling, attenuated by depth, plus binomial read sampling, summed across independent bulks.

Parent handling is sound in the two-parent modes: sites are kept only when the parents confidently disagree (≥85% allele purity), bulk genotype quality is not used as a filter (it is allele-frequency biased in a pool), and a depth cap for collapsed repeats is available.

## Statistical errors that change calls

### 1. The CompositeZ threshold is the null of a different statistic

The scanned score is `|Z(ΔSI)| + Z(G) + Z(LOD)`, divided by √3:

```go
// stats/smoothing.go
czs = append(czs, math.Abs(s.ZDeltaSI), s.ZGstat, s.ZLOD)
```

The Monte Carlo null instead robust-Z’s the raw absolute difference:

```go
// stats/thresholds.go
// col 0: |DeltaSI|
mat[i*k+0] = math.Abs(hSI - lSI)
mat[i*k+1] = GStatistic(...)
mat[i*k+2] = lod(...)
```

`|Z(ΔSI)|` and `Z(|ΔSI|)` are different random variables. For a signed Δ centered at 0, `|Z|` is large in both tails and near 0 at the null. `Z(|Δ|)` is negative when the difference is smaller than the typical absolute difference, and large only in one tail. G and LOD are also not always ≥ 0 after robust Z: about half the genome sits below the median and gets a negative Z. CompositeZ itself is therefore not one-sided, despite the comments in `consolidate` and `calculateEmpiricalCompositeZThresholds`.

The null is also built from unsmoothed single-SNP draws at the genome-wide median depth. The score that is thresholded has already been Gaussian-smoothed, then robust-Z’d on the genome-wide smoothed distribution. Smoothing changes skewness of G and LOD, and one median depth does not represent low- and high-coverage sites. `validateCompositeZStatCount` only checks that both sides use k = 3, so this mismatch is invisible to the test.

### 2. Raw single-SNP thresholds are applied to smoothed statistics

`simulateTwoBulk` draws one SNP. `FindPeakIntersections` compares those quantiles with `SmGstat`, `SmDeltaSI`, `SmLOD`, and `SmED`. A sliding kernel with σ = window/2 (default σ = 1 Mb, cutoff 3σ) shrinks the null variance, especially the read-sampling part. The raw 99.5% threshold is then wider than the null of the smoothed track, so individual-statistic and BRM scans lose power wherever sequencing noise, rather than pool sampling, dominates. At very high depth the pool term is shared by linked SNPs and smoothing does not remove it, so the same threshold is closer to calibrated — which means calibration depends on coverage.

G, LOD, and ED⁴ should be recomputed from smoothed allele frequencies. Smoothing those statistics directly is not the same operation: ED⁴ = (Δ)⁴ is convex, so averaging it keeps noise that would cancel in smoothed ΔSI (Jensen). The usual ED pipeline raises a smoothed distance to the fourth power; it does not smooth the fourth power.

### 3. “P99” is not the same tail probability for every statistic

```go
// stats/thresholds.go
DeltaSIP99: q(dsiArr, 0.995), DeltaSIP95: q(dsiArr, 0.95),
DeltaSIMp99: q(dsiArr, 0.005), DeltaSIMp95: q(dsiArr, 0.05),
GstatP99: q(gsArr, 0.995), GstatP95: q(gsArr, 0.95),
```

For two-sided ΔSI, 0.995 and 0.005 are a 99% two-sided test, but the “P95” pair is 5% in each tail (90% two-sided), not 95%. G, LOD, ED⁴, and the Bayes factor are one-sided, so the column labeled P99 is the 99.5th percentile. CompositeZ uses the 99th. With the default 1,000 replicates, a 99.5% quantile is determined by about five tail draws, and that noisy estimate is reused for every SNP at that depth.

### 4. Pointwise thresholds are treated as QTL significance

Nothing controls the genome-wide chance of a smoothed bump. Neighboring SNPs in one pool share the same individuals, so a null excursion is a region, not a point. A chromosome-scale simulation (draw one pool genotype configuration, drop reads, smooth, record the maximum CompositeZ) is the threshold that matches the scan. Bulk-label permutation is the non-parametric version of the same idea.

### 5. Reported Bayes factors and G statistics are read-level, not pool-level

`betaBinomialLogBF` marginalizes a Beta(0.5, 0.5)–binomial over reads. It does not marginalize the finite-pool genotype distribution. At 200× and n = 20, a very small frequency difference produces a huge log₁₀ BF. The Monte Carlo threshold can still rank sites, but a value of 10 is not “strong evidence” in the Kass–Raftery sense, and G is not a χ² test on reads. The README’s “LOD > 3 = significant” is the human-pedigree rule of thumb; this program never uses it.

## Biological model issues

**Backcross-F2 pools are simulated as backcrosses.** `ExpectedAF("BC1F2H")` correctly returns 0.75. `samplePoolAF` parses the same string, then draws only recurrent homozygotes and heterozygotes with heterozygosity ½ⁿ:

```go
// stats/thresholds.go
if strings.HasPrefix(p, "BC") {
    suffix := p[2:]
    isHigh := strings.HasSuffix(suffix, "H")
    suffix = strings.TrimSuffix(strings.TrimSuffix(suffix, "H"), "L")
    suffix = strings.TrimSuffix(suffix, "F2")
    n, err := strconv.Atoi(suffix)
    if err == nil && n >= 1 {
        het := int(distuv.Binomial{N: float64(bulkSize), P: math.Pow(0.5, float64(n)), Src: rng}.Rand())
        copies := het
        if isHigh {
            copies = 2*bulkSize - het
        }
        // ...
    }
}
```

A BC1 individual is AA or Aa. A BC1F2 individual is AA, Aa, or aa at 0.625 / 0.25 / 0.125. The variance of allele frequency for one BC1F2 plant is twice the BC1 variance (0.125 vs 0.0625). Thresholds for every `BC{n}F2*` population are too tight, so false intervals are more likely. `PopulationVarianceScale` already has the F2-after-backcross factor; the sampler should use that same three-genotype distribution.

**ED⁴ is (ΔSI)⁴, not Euclidean distance to the fourth.** For a biallelic site the frequency-vector distance is √2 |Δ|, so ED⁴ = 4 (Δ)⁴. Calling and thresholding use the same function, so the scan is self-consistent, but the column cannot be compared with published ED⁴ cutoffs.

**Unbalanced bulks are collapsed to the shallower depth.** `BSAstats.Depth` becomes the minimum of the two allele totals, and both bulks are then simulated at that depth. The deeper bulk’s true sampling variance is smaller than the threshold assumes.

**A non-target ALT can be accepted as the target allele.** In `parentAllele`, any non-reference fraction ≥ 85% returns `targetAlt`, even when those reads support a different ALT. The site can pass as a confident parental contrast for an allele the parent does not carry. The check needs the target-ALT fraction itself, not “everything that is not REF”.

**One-bulk BRM treats a smoothed deviation of exactly 0 as missing.**

```go
// stats/brm.go
if s.SmAFDev != 0 {
    deviation = s.SmAFDev
} else {
    deviation = s.SmHighSI - expectedSI
}
```

In a low-bulk-only run `SmHighSI` is never set. An exact zero then becomes `−p₀` and can open a false block. Use a presence flag, or always use `SmAFDev`.

**Gene-space polarity is hard-coded as resistance.** High parent / high bulk are “resistant” and low parent / low bulk are “susceptible”. For `lplb`, the low parent is labeled resistant and the low bulk susceptible. That mapping is only meaningful for a resistance bulk.

**Equal DNA contribution is assumed.** There is no overdispersion for unequal tissue input, related individuals, or library bias. Real pools are wider than this null; thresholds will be somewhat anti-conservative on noisy libraries.

The README’s account of the two-stage null is reversed. Drawing reads straight from p₀ has variance p(1−p)/depth. Adding a finite pool increases variance. The two-stage model is the wider, more accurate null; it does not reduce false positives relative to a reads-only null. It removes false positives relative to interpreting read-level G or BF with χ² or Kass–Raftery cutoffs.

## Code correctness

- **One-bulk cache key omits the population.** `oneBulkCacheKey` hashes depth, p₀, replicate count, and bulk sizes. F2 and RIL both have p₀ = 0.5 and different variances. A single CLI run has one population, so this is latent, but the two-bulk key already includes the population and the one-bulk key should too.
- **Thresholds are not reproducible.** Seeds come from a shared counter, and depth simulations run on all CPUs. The same VCF can yield different quantiles on different runs. A scientific pipeline needs a user seed and a per-depth RNG that does not depend on scheduling.
- **Merged lower-tail peaks keep the weaker summit.** `mergePeakIntersections` keeps the larger `PeakValue`. For ΔSI⁻ the significant peak is the more negative one. Consolidated QTL rows also keep the left-hand peak when two intervals merge, rather than the stronger one.
- **QTL chromosome order is map iteration order** in `RunBRM` and `DetectIndividualQTLs`, so output order changes between runs.
- **If bulk size is 0, the two simulations disagree.** Per-statistic simulation substitutes depth; CompositeZ substitutes 100. The CLI default (`20,20`) hides this. The fallback should be one function, and a missing size should be a hard error: it dominates the null at high depth.
- **Mode detection is substring logic** (`strings.Contains(bsaType, "hb")`). Current names happen to work. A switch on the mode, shared by filtering, stats, and TSV writing, would keep the next mode from silently dropping a column.
- **Tests do not lock the science.** They check peak merging and that an F2 dosage is {0, 0.5, 1}. There is no fixture for G, LOD, ED⁴, expected AF, variance scale, `|Z|` versus `Z(|·|)`, or BC1F2 genotype frequencies.

## Documentation

The Magwene citation is the wrong 2011 paper. The G statistic used here is Magwene, Willis, and Kelly, *PLoS Computational Biology* 7:e1002255 (2011), “The statistics of bulk segregant analysis using next generation sequencing,” not the G3 paper listed in the README. Morton (1955) is a pedigree linkage LOD; the LOD here is a two-binomial likelihood ratio on read counts. The BRM formula in the README also omits the depth term that `calculateBRMBlocksTwoBulk` actually uses.

## Improvements, in order

1. Make the null identical to the scan: one pool draw, binomial reads, the same Gaussian smooth, the same robust Z, the same CompositeZ formula (`|Z(ΔSI)|`, not `Z(|ΔSI|)`). Take the threshold from the maximum of that null scan, or from bulk permutations. Use at least 10⁴ replicates and a fixed seed.
2. Smooth allele frequencies (depth-weighted), then calculate ΔSI, G, LOD, ED⁴, and the Bayes factor from the smoothed frequencies. Threshold those derived values.
3. Give `BC{n}F2` three genotypes in `samplePoolAF`, sharing the variance already coded in `PopulationVarianceScale`. Simulate each bulk at its own depth.
4. Fix `parentAllele` so support is the target ALT fraction. Treat exact-zero `SmAFDev` as a real zero. Put the population into the one-bulk cache key. Require `--bulk-sizes` when depth is much larger than the default.
5. Stop interpreting read-level BF and G as classical evidence. Either marginalize the pool size in the likelihood, or label the columns as scores whose only threshold is the simulation.
6. Add table tests: known 2×2 G and LOD, ED⁴ = 4(Δ)⁴ if you adopt the literature definition, AF and variance for F2 / RIL / BC1 / BC1F2, and a small chromosome where the simulated CompositeZ quantile matches the code path that emits `CompositeZ`.
