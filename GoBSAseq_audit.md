# GoBSAseq audit

**Repository:** [gmaffy/GoBSAseq](https://github.com/gmaffy/GoBSAseq)
**Scope:** biological soundness, code correctness, accuracy and speed of VCF filtering, raw statistics, smoothing, normalising, consolidating, QTL detection and plotting; comparison with current flagship BSA-seq tools; review of `nomalising_ideas.md`.
**Date:** 2026-11-21

This is a source audit. `go test` could not be executed in the review environment (no Go toolchain). Ratings come from the current source, the unit tests as written, and the formulas, not from a fresh runtime bake-off on a shared VCF. Calibration numbers cited below are the repository's own (`TestHMMLODCalibration`, opt-in via `GOBSASEQ_CALIBRATION=1`) and the figures in `nomalising_ideas.md`.

## Verdict

GoBSAseq is statistically ahead of the usual BSA-seq packages on error control, and behind them on one thing that currently gates the final call: robust-Z CompositeZ.

Trust **HMM LOD + BRM** for calls. Treat CompositeZ as a display score until it is restandardised. Do not compare a "Z" from this tool to 1.96 or 3.

| Stage | Biology | Code | Speed | Notes |
|---|---|---|---|---|
| VCF filtering | 8 | 8.5 | Excellent | Right BSA gates; GATK filters are diploid defaults |
| Raw statistics | 9 | 9 | Excellent | SI, G, ED⁴, effective depth, one-bulk p₀ are correct |
| Smoothing | 8 | 8.5 | Excellent | Right kernels; G′ is not Magwene's G′; bandwidth is bp |
| Normalising | 6 | 8.5 | Fine | Robust Z is implemented correctly and is the wrong scale |
| Consolidating | 5.5 | 9 | Fine | CompositeZ is calibrated and weak |
| BRM | 8.5 | 9 | Cheap | Better genome-wide threshold than Huang et al. 2020 |
| HMM LOD | 8.5 | 8.5 | The bottleneck | MULTIPOOL-class, plus a real threshold |
| QTL calling | 8 | 9 | Fine | Careful rules; consensus is limited by CompositeZ |
| Plotting | 6.5 | 7 | Heavy HTML | Real thresholds, distorted axis |

## What is solid

Population moments in `utils/population.go` check out. For a bulk of *n* diploids, pool-frequency variance is p₀(1−p₀)/(c·n) with:

| Population | p₀ (high-parent allele) | c (variance scale) |
|---|---|---|
| F2 | 0.5 | 2 |
| F3 | 0.5 | 4/3 |
| F4 | 0.5 | 8/7 |
| RIL | 0.5 | 1 |
| BC1 to high | 0.75 | 3 |
| BC1F2 to high | 0.75 | 1.5 |

These match direct dosage calculations (F2, F3, F4, RIL, BC1, BC1F2). Selfing does not change allele frequency, so BCnF2 correctly shares BCn's p₀. The general BCn and BCnF2 variance formulas match the BC1 / BC1F2 cases.

Raw statistics in `stats/stats.go`:

- **SI** is the high-parent allele count over REF+ALT depth. Polarity is fixed once in the filter.
- **G** is the 2×2 likelihood-ratio statistic with no pseudo-counts, symmetric in allele order (Magwene et al. 2011).
- **ED⁴** = 4·ΔSI⁴, which is (√2·|ΔSI|)⁴, the biallelic Hill et al. 2013 definition.
- **Two-bulk depth** used for weights is the harmonic mean 2·dH·dL/(dH+dL): the equal depth that gives ΔSI the same read-sampling variance. That is the right reduction.
- **One-bulk p₀** estimates reference-mapping bias as half the gap between the interquartile means of sites where the high allele is REF versus ALT, using only sites with at least 2 reads of each allele. A QTL moves both classes together and cancels. Falling back to the theoretical p₀ below 200 sites per class is the right guard.

The genome-wide null is the right threshold design: real markers, depths and bulk sizes; linked pool frequencies; beta-binomial reads at ρ estimated from neighbours at most 0.1 cM apart; the same scoring code as the data; a generalised extreme-value fit to genome maxima. The old pointwise "p99" is gone. On the repository's own plant-level calibration, type I error for the HMM LOD, BRM and CompositeZ is in good shape. That calibration was not re-run for this review.

## VCF filtering

**Biology 8 / code 8.5 / speed excellent.**

`filter/fastfilter.go` is a streaming parser with GATK VariantFiltration semantics (a present annotation beyond its threshold fails; a missing annotation passes). BSA-specific gates are the right ones:

- Bulk depth is REF+ALT allele depth, not FORMAT/DP. Reads on other ALTs may be at most 5% of allele depth.
- No GQ gate on bulks. In a pool, GQ tracks allele frequency and would bias the SNP-index. Parent GQ is applied.
- Parents must be homozygous (85% of allele depth on one allele; other-ALT reads count against both) and, in two-parent modes, must differ.
- A segregation floor (default 0.4 × the expected minor allele frequency) applies only in two-bulk modes without two informative parents (`2b`, `hp2b`, `lp2b`), where nothing else proves the site segregates in the cross.
- Multi-allelic records are split, and the other-ALT depth is kept. Duplicate positions contribute only their deepest view to smoothing and the HMM.
- `--light-filtering` drops MQRankSum and ReadPosRankSum, which are unreliable near allele fixation in pools.

Two limits, neither of them a formula bug:

1. GATK QD/FS/SOR/MQ were built for diploid germline calls. In a cohort VCF they describe every sample, not the bulks. The README says so. Rank-sum tests are correctly optional.
2. Parent homozygosity at 85% still accepts a mildly heterozygous or contaminated parent. At 100× that is 15 alternate reads called homozygous. A binomial homozygosity test, or a 95% floor, would be tighter. True heterozygotes (~50%) are already rejected.

The filter architecture (kept-sample columns only, parallel workers, exact BGZF virtual offsets, own tabix writer) is fast. The repository's own figure, 300k sites × 60 samples with 4 used, is 0.58 s and ~0.5 GB RSS versus 13 s and 3.6 GB for the old vcfgo path. That was not re-benchmarked here. The design supports the claim.

## Raw statistics

**Biology 9 / code 9 / speed excellent.**

SI, ΔSI, G, ED⁴, effective depth and the one-bulk G-test are correct, stored at full precision, and rounded only on write. ED⁴ uses significant digits so small values are not printed as 0. A site with the wrong sample count is an error, not a silent skip. Workers format TSV chunks without `fmt.Sprintf` per site.

Remaining biological limits, all downstream of the per-site formulas:

- Per-site G and ΔSI assume binomial reads. Overdispersion is handled in ρ and in the null, not inside G. That is acceptable because thresholds are simulated, but a raw G is not a χ² p-value on overdispersed pools.
- `2b` without `--ref-parent` stores |ΔSI|. The sign would cancel under smoothing if REF is sometimes the high parent and sometimes the low parent. G and ED⁴ do not need polarity. The HMM does, and it is skipped. That is the mode most likely to be run without parents, and it loses the best statistic.
- BCnF2 map expansion is `n+1` (`utils/population.go`). The extra selfing meiosis only recombines heterozygotes (frequency 2⁻ⁿ), so the expansion should be nearer `n + 2⁻ⁿ` (BC1F2: 1.5, not 2). F3 (1.5), F4 (1.75) and RIL (2) already use the heterozygote-weighted series. This affects the HMM, not SI or G.

## Smoothing

**Biology 8 / code 8.5 / speed excellent.**

Kernel choices match the literature:

| Track | Kernel | Literature |
|---|---|---|
| SI, ΔSI, AFDev | Inverse-variance Gaussian, σ = window/2; local linear within 2σ of ends on chromosomes spanning at least 6σ | Better edge behaviour than QTLseqr's tricube |
| G, G1 | Tricube over ±window | Magwene G′, QTLseqr |
| ED⁴ | Local quadratic, tricube weights, clamped at 0 | MMAPPR loess, but a fixed bandwidth rather than an AICc span |
| HMM LOD | Not smoothed | Correct: it is already linkage-aware |

Markers are weighted by the precision of the allele-frequency estimate, 1 / Σ_bulks [1/(c·n) + 1/depth], so once depth exceeds the pool size, deeper markers gain little. That is the right weight for SI and ΔSI. Split multi-allelic views at one position count once (deepest view). Grid evaluation plus σ/100 bins makes cost follow chromosome length rather than marker density. The README's figure is 2M markers, smoothing from 114 s to about 1 s. The approximation claim (within 1% of each track's range versus exact evaluation on pure noise; CompositeZ moves by at most 0.05) is tested, not re-run here.

Three biological gaps:

1. **G′ is not Magwene's G′.** The tricube weight is multiplied by the allele-frequency precision weight (`smoothing.go`, `tw := nwt * tricube`). Per-site G already grows with depth, so deep sites are counted twice. Thresholds stay valid because the null uses the same weights. The track will not match QTLseqr, which is the G′ reviewers expect. Drop the marker weight on the tricube and loess kernels.
2. **Bandwidth is base pairs.** `--cm-per-mb` affects only the HMM. Under the HMM's default (1 Morgan per chromosome), a 2 Mb window is about 2 cM on a long chromosome and about 10 cM on a short one. BRM, G′ and the HMM are not on the same genetic scale, and the final call requires them to overlap.
3. **Interior smoothing is a local mean**, so peaks are flattened. Local linear is only an edge patch. Loess already does the right thing for ED⁴; SI and ΔSI should have the option.

NEff, (Σw)²/Σw², is computed from the Gaussian kernel and then applied to every track. A marker can pass that gate and still have a thin tricube window. Minor.

## Normalising

**Biology 6 / code 8.5.**

The implementation of the robust Z is correct and careful:

- Centre is the median, scale is 1.4826 × MAD (Rousseeuw & Croux consistency factor), with a mean/SD fallback when MAD is 0.
- For non-negative tracks (G′, G1, ED⁴) the **displayed** Z is of log(x + c), c = 10% of the median, which cuts skew.
- CompositeZ uses the **untransformed** robust Z of G′ and G1. The code comments say the log-based version raised the pointwise false-positive rate on overdispersed data (1.4% vs 1.0%). That split is intentional and easy to misread off the TSV.

The scale is the wrong one for inference. Z = (x − median) / (1.4826·MAD) is taken over the whole genome. Smoothed tracks are a handful of linked blocks, so the MAD is a small, noisy yardstick. A null genome's maximum robust Z lands around 15–25, not 3. The same QTL then scores differently in a quiet genome and a messy one (`nomalising_ideas.md`: CompositeZ 18.8 against a threshold of 18.4 on one simulated set, and 395 against 20 on another).

What the Z columns hide: individual-statistic thresholds and the BRM threshold are quantiles of the smoothed values (or of |signal|/σ) under the null. They do not inherit the robust-Z scale. A G′ or ΔSI region can be well calibrated while CompositeZ is not.

## Consolidating

**Biology 5.5 / code 9.**

CompositeZ is Stouffer on two components:

- Two-bulk: (|Z_ΔSI| + Z_G′) / √2
- One-bulk: (|Z_AFDev| + Z_G1) / √2
- `2b` without polarity: (Z_|ΔSI| + Z_G′) / √2, without folding again

The HMM LOD is correctly excluded. The code does exactly what it says, and the null simulation uses the same combination, so the false-positive rate is controlled. The combination itself is weak:

- G′ and ΔSI are two views of the same reads. Stouffer's √2 assumes independent tests. It does not buy a second piece of evidence.
- |Z_ΔSI| is folded (half-normal). Z_G′ is right-skewed and can be negative. The sum has an awkward null, which the simulation absorbs by raising the threshold.
- Because the scale depends on genome composition, CompositeZ cannot be compared across studies.

ED⁴ is a monotone function of ΔSI at the site level (4·Δ⁴) and is left out of CompositeZ. That is fine; putting it in would be a third copy of the same contrast.

## HMM LOD

**Biology 8.5 / code 8.5 / this is the runtime bottleneck.**

`stats/hmm.go` is a MULTIPOOL-style pool HMM:

- Hidden state is the high-allele count *j* in U = round(c·n) exchangeable units. The prior is Binomial(U, p₀), which matches the mean and the pool variance, and matches the distribution shape for F2 and RIL.
- A unit switches high→low with probability 2(1−p₀)r and low→high with 2·p₀·r, r the Haldane fraction times map expansion. p₀ stays stationary. The autocorrelation 1−2r is exact for F2.
- One-bulk LOD tests the bulk against its prior around p₀ (with the mapping-bias shift).
- Two-bulk LOD tests "bulks differ" against "both bulks drawn around one shared frequency". Segregation distortion that shifts both bulks is not scored as a QTL. That alternative is better than a naive per-bulk MULTIPOOL.
- Markers within 0.01 cM share one HMM step. Likelihood ratios are normalised so the per-bin constants cancel. The shared-frequency maximum is a 41-point grid plus a parabolic refinement. Real data and the null use the same function.

Limits:

- Read error is hardcoded at 1% and is not a flag. Modern Illumina after filtering is lower; the null uses the same value, so the test stays internally calibrated, but the model is slightly misspecified if the truth is ~0.1%.
- The likelihood is binomial. ρ inflates the threshold, not the LOD. A LOD of 5 is not a classical LOD of 5 once overdispersion is present. Compare LOD to the genome-wide line, never to 3.
- BCnF2 map expansion is high (see raw statistics).
- No HMM in `2b` without `--ref-parent`. A phase-free joint model was tried and dropped (about 20× slower; kept under `oldcode/`).
- Forward/backward state is stored as float32. Fine at usual bulk sizes; worth watching if U is in the hundreds and depths are very high.

## BRM

**Biology 8.5 / code 9 / cheap once smoothing has run.**

This is block regression mapping (Huang et al. 2020), not a Bayesian regression. The README's opening sentence still says "Bayesian regression modeling". The pipeline section names it correctly.

The signal is smoothed signed ΔSI, or √max(smoothed bias-corrected ΔSI², 0) without polarity, or smoothed AFDev in one-bulk modes. σ² is the single-locus pool variance plus the kernel's read variance (ReadVar0 + ρ·ReadVar1) accumulated in the same pass as the smooth. The threshold u is the (1−α) quantile of the genome-wide maximum |z|, so `--brm-alpha` is a genome-wide error rate. Blocks with fewer than `--min-qtl-markers` markers are dropped. Same-sign blocks within the merge gap are joined, keeping the stronger peak.

That threshold is a real improvement on the original paper's Bonferroni-on-map-length factor, which ignores smoothing and depth. The analytic σ matches the smoother, including the local-linear blend near ends.

Still open, and already written up in `brm_todo.md`:

- Pool variance is the single-locus variance. Exact when the window is one linked block; conservative when the window spans many centimorgans.
- No position confidence interval. The original BRM returns one. The called interval is the run of markers above the line, which grows with effect size and is not a 95% interval for the QTL position. The HMM 1.5-LOD support interval on the final QTL is a partial substitute.
- Unequal DNA contribution shrinks the effective bulk. Every threshold that uses *n* is then a bit low. ρ from tightly linked neighbours does not catch it, because both neighbours share the same pool.

## QTL detection

**Biology 8 / code 9.**

The three lists do what a careful caller wants:

- Each smoothed statistic is scanned at its genome-wide p99. p95 is used only when that statistic has **no** p99 region anywhere, and those regions are labelled suggestive. The fallback is not per chromosome.
- Opposite tails of ΔSI or AFDev are never joined. Unsigned statistics join the signed QTL they overlap most, or form a `none` QTL.
- The strongest peak is kept, not the first or the last.
- A final locus is the overlap (or merge-gap contact) of HMM LOD regions, BRM blocks and CompositeZ regions. It is a consensus QTL only when every available method is present. `2b` without polarity has no HMM, so the requirement is BRM + CompositeZ.
- Core (intersection) and span (union) are both reported. Gene space uses the core unless `--gene-space-span`.
- Partial-support loci stay in the file with a `Support` column.
- Lead position is the midpoint of the HMM LOD's 1.5-LOD support interval, which is the right choice when the LOD is a flat plateau.

The cost is the consensus rule. The calibration QTL is complete homozygote selection (high bulk AA only, low bulk aa only) in 20-plant bulks. That is the strongest contrast a bulk can show. CompositeZ only just clears its line (18.8 vs 18.4) while the HMM LOD is about 24 vs about 5 and BRM finds it every time. `TestHMMLODCalibration` asserts that BRM finds that QTL. Whether the consensus finds it is only logged. Requiring CompositeZ gives the weakest method a veto. Partial QTLs will miss.

## Plotting

**Accuracy 6.5 / code 7.**

`plots/plots.go` draws the real genome-wide p99 and p95 lines, not a hardcoded Z = 3. SI in `2b` without polarity is labelled as a REF fraction and drawn without thresholds. BRM blocks are shaded. That part is right.

Three accuracy problems:

1. **The x-axis is a category axis of position strings.** A SNP desert and a SNP pile have the same width per point. QTL width on the figure is not base pairs. QTLseqr plots on a linear physical axis. This is the main plotting defect.
2. **`Smooth: true` splines a series that is already kernel-smoothed**, so the line the user sees is not the series that was tested.
3. **The tooltip marks p99 only when `value > upper line`.** A lower-tail ΔSI QTL is never starred, even though the valley lines are drawn. MaxAbsZ is judged against the CompositeZ threshold, which is a different statistic.

`zSig = 3` and `zSugg = 2` are unused. On this scale, 3 is not a genome-wide 1% line. Per-chromosome HTML with one point per marker will be large; there is no downsampling and no genome-wide Manhattan.

## CLI bugs that will corrupt or block runs

These are not statistical, but they sit in front of every analysis:

1. **`run.Run` prompts on stdin for gene-space files whenever snpEff, GFF, PRG or descriptions are absent**, which is the default. A batch run blocks or fails before any statistic is computed. Gene space should be opt-in.
2. **`-P None` is documented as "skip this parent"** and is then looked up as a sample name. `None` is not in the VCF, so the run drops into an interactive prompt. Non-interactive `2b` and one-parent runs do not work as documented.
3. **Read input (`--hp-reads` and the other read flags) returns success and does nothing.**
4. **README line 1 calls BRM Bayesian.** The code is Huang's block regression. Reviewers will look for the wrong paper.

`BulkSizes` decides which bulk a mode uses with `strings.Contains(bsaType, "2b")` and similar. It happens to be correct for the current mode names (`2p2b`, `2phb`, `hplb`, …). It is a footgun if a mode name is added.

## Speed, relative to the flagships

Not a wall-clock bake-off on one VCF.

Filtering and smoothing are already in flagship territory. The bill is `--null-reps` (default 200) full HMM passes plus a smooth of each replicate. Dense SNPs stay cheap: a 0.01 cM bin has almost no crossovers, so each forward step is O(U) with U = c·n, not O(U²). Large bulks with large inter-marker gaps are the expensive case, not marker density. QTLseqr and the original BRM stay faster because they never resimulate an HMM. That extra time is buying the only threshold in this set that is both genome-wide and run through the same code as the data.

## Benchmark against current tools

| Tool | What it tests | Threshold | Where GoBSAseq lands |
|---|---|---|---|
| [QTLseqr](https://github.com/bmansfeld/QTLseqr) v0.7+ | Tricube ΔSNP-index and Magwene G′ | Takagi-style simulation, or a trimmed log-normal G′ then FDR | Same two statistics. GoBSAseq's null is the better error rate (chromosome-scale linkage, extreme value, ρ). QTLseqr's G′ is the one reviewers recognise; the depth-weighted G′ here will not match it. This is still what most plant BSA papers report. |
| MULTIPOOL (Edwards & Gifford 2012) | Pool HMM LOD | No genome-wide threshold. That was the published limitation. | GoBSAseq is the practical successor: F2/F3/F4/RIL/BCn, a two-bulk shared-frequency alternative, and a threshold. Expect similar power and fewer false intervals. It is not a line-by-line port. |
| BRM (Huang et al. 2020, [huanglikun/BRM](https://github.com/huanglikun/BRM)) | Smoothed allele-frequency difference | Analytic genome-wide factor from map length. Very fast. Also returns a position CI. | Same signal, a stricter and more honest threshold, no CI yet, and 200 HMM replicates instead of a closed form. Use GoBSAseq when false positives matter; use the original when you want their CI and a short run. |
| MMAPPR / MMAPPR2 | ED⁴ loess, span chosen by AICc | "Above the fit's noise" | The ED⁴ formula and the loess match the paper. The span does not adapt. The significance line is better. |
| [PyBSASeq](https://github.com/dblhlx/PyBSASeq) | Windowed SNP-index, G, ED, Fisher exact | Largely descriptive | Weaker model. GoBSAseq is ahead on every inferential point. |
| EXPLoRA (Pulido-Tamayo et al. 2016) | HMM with explicit multiple testing | Has a significance rule | Same idea as an HMM plus a genome-wide null. Poorly maintained. GoBSAseq covers this role if the HMM stays. |
| [bsaseq](https://github.com/rcac-bioinformatics/bsaseq) (Purdue, 2025) | Tricube window, Z and G | Hard Z = 3 | Fast (tens of seconds on a sorghum genome) and, on this null, anti-conservative. A genome-wide robust-Z maximum here is about 15–25. Z = 3 will call noise. |
| SHOREmap, NGM, SIMPLE | Allele frequency and candidate intervals for mutants | Heuristic | A different job (mapping a single induced mutation). One-bulk GoBSAseq competes with MutMap on AF deviation and has a real null. Those tools stay simpler for EMS screens. |

Shen et al. (and the Li & Xu 2022 review in *The Plant Journal*) already showed that ED, |ΔSNP-index| and G′ agree on major QTLs and diverge on the margins. GoBSAseq's contribution is not a new score. It is a null that those scores have usually lacked, plus an HMM those scores do not replace. The consensus rule currently throws away part of that advantage by letting CompositeZ veto the HMM and BRM.

## Review of `nomalising_ideas.md`

The filename is misspelled (`nomalising_ideas.md`). The diagnosis is right, and it matches `normalise` and `consolidate`.

What the note gets right:

- The scale comes from the genome's own median and MAD, not from a null model. Linked smooth tracks have few independent blocks, so null maxima in MAD units are large and the genome-wide line sits near 15–25.
- The scale depends on how much of the genome is null, so the same QTL is not comparable across datasets.
- G′ and ΔSI are near-copies. Stouffer gains nothing.
- |Z_ΔSI| and Z_G′ do not share a shape. The simulation hides this by raising the threshold.
- BRM's z is the textbook standardisation: null variance near 1, genome-wide line near 2.9.
- Do not turn either component into a normal-table p-value. Keep simulated genome-wide maxima.

One correction of emphasis. Individual G′, ED⁴ and ΔSI thresholds are already null quantiles of the smoothed tracks. They do not use the robust Z. The damage is specifically CompositeZ and the `Z_*` columns, not every significance line in the output.

On the three proposals:

### 1. Null-SD standardisation with a correlation-aware Stouffer

This is the right repair **if a single ΔSI+G′ score must stay**.

- z_ΔSI = smoothed ΔSI / σ_null, which is the BRM z already computed (pool variance + ReadVar0 + ρ·ReadVar1).
- z_G′ = (G′ − E₀[G′]) / SD₀[G′] from the null genomes already being simulated. Two hundred replicates put a per-marker SD within about 5%. Do not convert it with a normal table; G′ stays skewed.
- ConsensusZ = (z_ΔSI* + z_G′) / √(2+2r), with r estimated on the null genomes.

Expect r to be very high. ConsensusZ will then sit almost on top of |z_ΔSI|. If r > ~0.9 the second term is not evidence, and the log should say so. The checks in the note are the right gates: null rate near 1/100 with and without ρ, in `2b` and in one-bulk modes; power against the current CompositeZ, BRM and the HMM. Add one more: rank correlation of ConsensusZ with BRM z. If the ranks match, a second method has not been gained.

### 2. Use BRM's z and drop CompositeZ from the call

This is what should ship.

G′, ED⁴ and ΔSI are monotone functions of the same allele counts. Another Z of them does not confirm the QTL. The only partly independent contrast in the pipeline is the HMM (linkage along the chromosome) against the kernel smooth (a local average). Final QTLs should require the HMM LOD and BRM, same direction, core = overlap, span = union. Keep G′, ED⁴ and ΔSI as their own tracks with the genome-wide thresholds they already have. Those are already the honest version of a QTLseqr plot.

In `2b` without `--ref-parent` there is no HMM. The call is then BRM alone. CompositeZ should not be promoted to fill that gap.

### 3. Keep CompositeZ, stop requiring it for the final QTL

Smallest diff, and it removes the veto the calibration already shows. It leaves a dataset-dependent number in the main table, which will be quoted. Not recommended as the end state.

### Naming

If the scale is replaced, call the new columns null-standardised and do not keep writing `Z_*` for a MAD Z. If both exist, use a distinct prefix (`NZ_DeltaSI`, `NZ_Gstat`). "ConsensusZ" is only honest when r is not essentially 1.

### Before adopting it

The note's calibration QTL is a best case (allele frequency 1 versus 0), and CompositeZ barely clears. Partial QTLs will miss. Fix the gate before adding statistics. After the change, do not require three highly correlated tests and expect the power of the best one.

## Priority order

1. Stop requiring CompositeZ for a final QTL. Call on HMM LOD + BRM. In `2b` without polarity, call on BRM.
2. Make gene space opt-in, and treat the sample name `None` as "not provided" so batch runs do not block on stdin.
3. Drop the marker-precision weight from G′ (plain tricube) so the track is Magwene's and comparable to QTLseqr.
4. If a composite score is kept, standardise by the null SD and correct for r, under a new column name. Report r.
5. Plot on a linear base-pair axis; stop spline-smoothing tested series; star the lower tail; do not score MaxAbsZ against the CompositeZ line.
6. Later: kernel bandwidth in cM, BCnF2 map expansion, a BRM position interval, and a pool-size inflation factor for unequal DNA contribution.
