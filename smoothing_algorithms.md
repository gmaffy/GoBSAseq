# Smoothing algorithms by statistic

As implemented in `stats/smoothing.go`.

| Statistic | Smoothing algorithm name | Description |
|---|---|---|
| HighSI | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | High-bulk high-parent allele fraction, averaged over markers within ±3σ (σ = window/2) with weight = marker precision × Gaussian distance. In `2b` without `--ref-parent` it is the REF fraction: plotted for QC only, with no thresholds. |
| LowSI | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | As HighSI, for the low bulk. |
| ΔSI | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | Per-site HighSI − LowSI, smoothed. In `2b` without `--ref-parent`, \|ΔSI\| is smoothed instead (`Sm_AbsDeltaSI`). |
| G | Tricube kernel (Magwene's G′), × marker weight | Per-site 2×2 G-test, tricube-weighted average over ±h with h = `--window-size` (QTLseqr's windowSize). |
| ED⁴ | Loess: local quadratic regression, tricube × marker weight, ±h | Per-site 4·ΔSI⁴; the fitted value at each marker (as MMAPPR's loess, but with a fixed ±window bandwidth instead of an AICc-chosen span). No edge bias; clamped at 0. |
| HMM LOD | None (pool hidden Markov model) | Computed along each chromosome by forward–backward over the hidden pool allele count (linkage and pool sampling built in), so it is not smoothed. Two-bulk: bulks differ (shared frequency free); one-bulk: bulk vs P0. Replaces the per-SNP LOD and Bayes factor. |
| AFDev | Inverse-variance-weighted Gaussian kernel; local linear within 2σ of chromosome ends | One-bulk: per-site SI − P0 (P0 corrected for mapping bias by polarity), smoothed. |
| G1 | Tricube kernel (G′-style), × marker weight | One-bulk: per-site one-sample G-test against P0, tricube-weighted average over ±h. |
| P0 | None | The site's own value is carried through unsmoothed. |
| Depth, HighDepth, LowDepth | None | The site's own values. They feed the marker weight and the threshold lookup. |
| NEff | Derived from the kernel | (ΣK)²/ΣK², the effective number of markers behind the smoothed value. Markers with NEff < 5 are dropped. |
| Z_* (each smoothed track) | None (normalisation) | Genome-wide robust Z: (Sm − median) / (1.4826·MAD), falling back to the ordinary mean/SD Z when MAD = 0. |
| CompositeZ | None (combination) | (\|Z_ΔSI\| + Z_G + Z_HMMLOD)/√3, or (\|Z_AFDev\| + Z_G1 + Z_HMMLOD)/√3 in one-bulk modes. In `2b` without `--ref-parent`: (Z(\|ΔSI\|) + Z_G)/√2. Thresholds from the genome-wide null simulation. |
| MaxAbsZ | None (combination) | max \|Z\| over all tracks. SI tracks are excluded in `2b` without `--ref-parent`. |

The marker weight is 1 / Σ_bulks [1/(c·n) + 1/depth], with c the population variance scale and n the bulk size. Gaussian kernel term: weight_m × exp(−½((pos_m − pos_k)/σ)²), σ = window/2, over ±3σ. Tricube term (G′, loess): weight_m × (1 − |Δ/h|³)³, h = window, over |Δ| < h.

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
