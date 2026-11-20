# Smoothing algorithms by statistic

As implemented in `stats/smoothing.go`.

| Statistic | Smoothing algorithm name | Description |
|---|---|---|
| HighSI | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | High-bulk high-parent allele fraction, averaged over markers within ±3σ (σ = window/2) with weight = marker precision × Gaussian distance. In `2b` without `--ref-parent` it is the REF fraction: plotted for QC only, with no thresholds. |
| LowSI | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | As HighSI, for the low bulk. |
| ΔSI | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | Per-site HighSI − LowSI, smoothed. In `2b` without `--ref-parent`, \|ΔSI\| is smoothed instead (`Sm_AbsDeltaSI`). |
| G | Tricube kernel (Magwene's G′), × marker weight | Per-site 2×2 G-test, tricube-weighted average over ±h with h = `--window-size` (QTLseqr's windowSize). |
| ED⁴ | Loess: local quadratic regression, tricube × marker weight, ±h | Per-site 4·ΔSI⁴; the fitted value at each marker (as MMAPPR's loess, but with a fixed ±window bandwidth instead of an AICc-chosen span). No edge bias; clamped at 0. |
| HMM LOD | None (pool hidden Markov model) | Computed in the raw-statistics stage along each chromosome by forward–backward over the hidden pool allele count (linkage and pool sampling built in), so it is not smoothed. Carried into the smoothed TSV in LOD units: no Z, not in CompositeZ or MaxAbsZ; it has its own genome-wide thresholds. Two-bulk: bulks differ (shared frequency free); one-bulk: bulk vs P0. Not computed in `2b` without `--ref-parent`. Replaces the per-SNP LOD and Bayes factor. |
| AFDev | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | One-bulk: per-site SI − P0 (P0 corrected for mapping bias by polarity), smoothed. |
| G1 | Tricube kernel (G′-style), × marker weight | One-bulk: per-site one-sample G-test against P0, tricube-weighted average over ±h. |
| ΔSI² (bias-corrected; `2b` without `--ref-parent` only) | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends (chromosomes ≥ 6σ) | (SI_H − SI_L)² minus each bulk's SI(1−SI)/(d−1), smoothed (`Sm_DeltaSI2c`); its square root estimates \|ΔAF\| without polarity and is what BRM tests in that mode. Not normalised; not in CompositeZ. |
| P0 | None | The site's own value is carried through unsmoothed. |
| Depth, HighDepth, LowDepth | None | The site's own values. They feed the marker weight and the threshold lookup. |
| NEff | Derived from the kernel | (ΣK)²/ΣK², the effective number of markers behind the smoothed value. Markers with NEff < 5 are dropped. |
| Z_* (each smoothed track) | None (normalisation) | Genome-wide robust Z: (value − median) / (1.4826·MAD), falling back to the ordinary mean/SD Z when MAD = 0. For the non-negative, right-skewed tracks (G′, G1, ED⁴) it is computed on log(x + c), with c = 10% of the median; on simulated no-QTL genomes this cuts the skew of G′/ED⁴ from +0.36/+0.49 to +0.14/+0.17. These log-based Z-scores are for display; CompositeZ uses the untransformed robust Z (see below). |
| CompositeZ | None (combination) | (\|Z_ΔSI\| + Z_G)/√2, or (\|Z_AFDev\| + Z_G1)/√2 in one-bulk modes. In `2b` without `--ref-parent`: (Z(\|ΔSI\|) + Z_G)/√2, with \|ΔSI\| not folded again. The HMM LOD is not included. G′ and G1 enter as their untransformed robust Z: with the log-based Z, the pointwise false-positive rate on overdispersed data rose from 1.0% to 1.4%. Thresholds are genome-wide: the 99%/95% quantiles of each simulated no-QTL genome's maximum CompositeZ (generalized extreme-value fit), with reads simulated at the read overdispersion estimated from the data. |
| MaxAbsZ | None (combination) | Largest deviation over the normalised tracks: \|Z\| for signed tracks (ΔSI, AFDev, one-bulk SI), the positive part of Z for one-sided tracks (G′, G1, ED⁴, \|ΔSI\|); the HMM LOD is not included. |

The marker weight is 1 / Σ_bulks [1/(c·n) + 1/depth], with c the population variance scale and n the bulk size. Gaussian kernel term: weight_m × exp(−½((pos_m − pos_k)/σ)²), σ = window/2, over ±3σ. Tricube term (G′, loess): weight_m × (1 − |Δ/h|³)³, h = window, over |Δ| < h.

Further details (all smoothed tracks):
- **Duplicate positions:** markers sharing a position (the views of a split multi-allelic site) are different allele contrasts from the same reads. Only the deepest view (most reads on its two alleles; ties go to the first) represents the position, in the kernels and in the HMM LOD, so the position counts once and contrasts aren't mixed. Every view still gets the position's output values.
- **Grid evaluation:** on chromosomes with more markers than grid points, the kernels are evaluated on a grid (`--step-size`, at most σ/20) and linearly interpolated to markers. Markers in sparse grid intervals (≤ 4 markers) are evaluated exactly.
- **Bin aggregation:** in grid mode, dense stretches enter the kernel as bins σ/100 wide (10 kb by default). Each bin sits at its weighted centroid and carries Σw, Σw² and the weighted-mean values, so the first-order position error cancels. On pure-noise input, the worst case, results are within 1% of each track's range of exact evaluation; Z-scores move by at most 0.09 and CompositeZ by at most 0.05.
- **Short contigs:** chromosomes spanning less than 6σ use only the local average; the edge local linear fit is not applied.

# Literature and best practice by statistic

Kernel shape matters little: Gaussian, tricube and Epanechnikov kernels give nearly the same results in kernel smoothing. Three things matter more:

- **Local linear vs local average.** The current smoother is a weighted local average, which flattens peaks and is biased at chromosome ends. Local linear regression (loess) corrects both.
- **Bandwidth in genetic distance, not bp.** In plant chromosomes recombination is far lower around centromeres, so a fixed bp window covers very different amounts of linkage in different places.
- **Whether to smooth the statistic at all.** For LOD and the Bayes factor, averaging is not the right operation.

| Statistic | What the literature does | Is the Gaussian kernel fine? | Best practice |
|---|---|---|---|
| SI / ΔSI | QTL-seq (Takagi 2013): simple sliding-window mean. QTLseqr (Mansfeld & Grumet 2018): tricube-weighted local average. | Yes; better than a plain sliding window. | A tricube or local linear kernel with the bandwidth in cM. For likelihood-based inference with confidence intervals, an HMM that models linkage directly (MULTIPOOL, Edwards & Gifford 2012). |
| G | G′ (Magwene et al. 2011): tricube-weighted local average of per-site G. | Yes for the smoothing itself. | Tricube G′, because Magwene's null procedure (log-normal fit after Hampel trimming, then false-discovery-rate control) was developed for it. It is the most widely used and validated BSA statistic. |
| ED⁴ | MMAPPR (Hill et al. 2013): loess regression of ED⁴. | Acceptable. | Loess, i.e. local regression, which also avoids the edge bias. |
| LOD | Smoothing a per-site LOD is not standard practice. | No: the operation is the issue, not the kernel. | Either drop it (it is G divided by a constant), or compute a model-based LOD (HMM or window likelihood) that accounts for linkage. |
| BBLogBF | No published BSA practice. | No: averaging log Bayes factors has no clear meaning, and summing them assumes independent markers, which linked markers are not. | Drop it, or replace it with HMM posterior probabilities. |
| AFDev (one-bulk) | MutMap (Abe 2012): sliding-window SNP-index for a single bulk. | Yes. | Same as SI. Edge correction matters for causal mutations near chromosome ends. |
| G1, LOD1, BBLogBF1 | No standard practice. | G1 yes; LOD1 and BBLogBF1 no, for the same reasons as LOD and BBLogBF. | Keep G1 (G′-style), drop LOD1 and BBLogBF1. |

## Recommendations for GoBSAseq

- **Keep the weighted kernel**, but switch to local linear for SI, ΔSI, AFDev and ED⁴ to remove the edge bias and peak flattening.
- **Allow the bandwidth in genetic distance**, either from a genetic map or a recombination rate in cM per Mb. Without a map, at least document that bp windows assume uniform recombination.
- ~~**Retire LOD/LOD1 and BBLogBF/BBLogBF1** as smoothed tracks.~~ Done: replaced by the HMM LOD.
- **Make G′ follow Magwene's method fully.** Tricube smoothing: done. Still to do: the trimmed log-normal null and FDR, as an option alongside the simulation-based thresholds.
- ~~**An HMM** (MULTIPOOL-style) for a calibrated LOD.~~ Done (`stats/hmm.go`, genome-wide null in `stats/hmmnull.go`). Confidence intervals for QTL position from the HMM posterior are still to do.

# Audit (2026-11-20): raw stats to consolidation

Ratings out of 10, for biological and code correctness and accuracy.

| Stage | Biology | Code |
|---|---|---|
| Raw stats | 8.5 | 9 |
| Smoothing | 8 | 8.5 |
| Normalising | 8 | 8 |
| Consolidating | 7 | 9 |

**Raw stats.** SI, ΔSI, the G-test, ED⁴, the effective depth and the one-bulk G-test are correct. The mapping-bias correction for one-bulk P0 is sound: it uses the interquartile mean, so QTL regions barely affect it. The HMM LOD is exact, with transitions per pool unit.
- Per-site statistics assume binomial read counts, so overdispersed data inflate G and ΔSI noise. Only the HMM thresholds take any account of this.
- `2b` without `--ref-parent` has no HMM LOD, which loses linkage-aware evidence in the mode most likely to be used without parents.

**Smoothing.** The kernel choices match the literature (Gaussian for SI/ΔSI, tricube G′, loess for ED⁴). Chromosome-end bias is corrected and split multi-allelic positions count once.
- The bandwidth is in bp, which assumes uniform recombination. `--cm-per-mb` only affects the HMM.
- G′ multiplies the tricube weight by the marker's precision weight, but G already grows with depth, so deep sites count twice. That departs from Magwene's G′.
- `smoothing.go` (about 1,400 lines: grid evaluation, bin aggregation, exact evaluation of sparse stretches) is fast and tested, but it is the hardest part of the app to maintain.

**Normalising.** Robust Z (median/MAD) is appropriate.
- It assumes most of the genome is null. Large or numerous QTLs widen the MAD and make the Z-scores conservative.
- Log-scale tracks have two Z-scores: the column shows the log-scale Z, while CompositeZ uses the untransformed one. This is correct but easy to misread.

**Consolidating.** The code is clean and driven by the track table.
- ~~CompositeZ thresholds come from a binomial null simulation, so overdispersed data (ρ = 0.01) give 1.43% of null sites over p99 instead of about 1%.~~ Fixed: the null uses the estimated overdispersion and the thresholds are genome-wide (1–3 of 100 null genomes over p99 in calibration).
- ΔSI and G′ come from the same reads and are highly correlated, so Stouffer's √2 is nominal and CompositeZ adds little beyond G′ alone. The empirical thresholds keep it valid.
- \|Z_ΔSI\| is always non-negative but Z_G′ can be negative, so the two components are not symmetric.

**Suggested fixes, in order**
1. ~~Add overdispersion, estimated from the data, to the null simulation.~~ Done.
2. Drop the marker weight from G′ (plain tricube, as Magwene).
3. Allow bandwidths in cM (a genetic map or `--cm-per-mb`) for the kernels.

# Thresholds (2026-11-20)

All thresholds now come from the genome-wide null simulation (`stats/hmmnull.go`); the single-SNP per-depth simulation is gone.
- **What:** for every scanned statistic (smoothed SI, ΔSI or |ΔSI|, G′, ED⁴, AFDev, G1; the HMM LOD; CompositeZ), the 99%/95% quantile of its maximum over a simulated no-QTL genome. Two-tailed tracks (SI, ΔSI, AFDev) split the level: 99.5%/97.5% of the maxima and of the minima. A p99 threshold has a 1% chance of being crossed anywhere in a genome without QTLs.
- **Why:** detection compares smoothed tracks, but the old thresholds were quantiles of a single SNP at that marker's depth, which is badly over-conservative for G′ and ignores linkage. The old ones were also inconsistent (99.5% for one-sided G/ED⁴), not reproducible (seeding depended on goroutine scheduling) and stored about 230 B per marker.
- **Null model:** real markers, depths and bulk sizes; linked pool frequencies from the HMM pool model; beta-binomial reads at the read overdispersion ρ estimated from the data (from tightly linked neighbours, so pool variation and QTLs do not inflate it); then the same scoring code as the data. Seeded, with per-chromosome streams: reproducible.
- **Quantiles:** a generalized extreme-value fit by probability-weighted moments (Hosking 1985) to the 200 genome extremes. A Gumbel fit understated the CompositeZ tail (q99 16.7 against 20.4 empirical); the GEV fit gives 22.3, against 22.1 from an independent plant-level simulator.
- **Calibration** (100 plant-simulated null genomes per case: 2p2b at ρ = 0 and 0.01, and 2b without `--ref-parent`): any detection track over its p99 in 0–2, over p95 in 2–5; CompositeZ 0–3 / 3–7; HMM LOD 0–1 / 1–4.
- `--rep` is deprecated and ignored.
