# Smoothing and normalising (`stats.SmoothAndNormalise`): todo

Scope: `stats/smoothing.go`. This covers `gaussianSmooth` (depth-weighted Gaussian kernel per chromosome), `normalise` (genome-wide robust Z-scores), `writeSmoothedTSV`, and the CompositeZ and MaxAbsZ fields that `consolidate` fills.
Sources: the code, `claude_review.md` (B1, C6), `grokreview.md` (#3, #4), `gemini_review.md` (2.4, 3.6, Priority 1), and measurements made in this review.
CompositeZ itself (LOD counted twice, mixing skewed Z-scores, its null) is **deferred** to a separate CompositeZ review. It is listed at the end only for completeness.

How it works now:
- **Grouping:** sites are grouped by chromosome.
- **Smoothed value:** at every marker k, each statistic is replaced by Σ wₘ·xₘ / Σ wₘ over the markers m within ±3σ, with wₘ = Depthₘ · exp(−½((posₘ − posₖ)/σ)²).
- **Kernel width:** σ = WindowSize/2 (default 1 Mb), so each smoothed value draws on 6 Mb.
- **Sparse markers:** a marker with fewer than 5 markers within ±3σ is dropped.
- **Z-scores:** each smoothed track becomes Z = (x − median)/(1.4826·MAD) over the whole genome.

Order: do Part 0 (golden test and benchmark) first. Then do Part 1 (results unchanged), then Part 2 (results change), one item at a time.

---

## Part 0: test harness first

- [x] **No tests exist** for `gaussianSmooth`, `normalise` or `robustZScores`. Add known-value tests:
  - a constant track stays constant;
  - a step function is smoothed symmetrically;
  - depth weights favour high-depth markers;
  - the sparse-marker drop works;
  - chromosomes are smoothed independently;
  - robust Z on a known vector, including the MAD = 0 fallback.
- [x] **Golden test:** fixed `[]BSAstats` input in all four layouts (two-bulk signed, two-bulk |ΔSI|, high-bulk only, low-bulk only). Compare `[]SmoothedStats` and the TSV byte for byte, so Part 1 can be shown not to change results.
- [ ] **Benchmark:** measured in this review on one 30 Mb chromosome with σ = 1 Mb. 25k markers took 6.7 s and 50k markers took 27.4 s, so the cost is quadratic in marker density.
  - A 400 Mb genome at one SNP per 600 bp: about 6 min.
  - 2M SNPs at one per 150 bp: about an hour.
  - Smoothing is now by far the slowest step after the filter.

## Part 1: code, speed and memory (results unchanged)

**Done (1.1–1.4):** one 30 Mb chromosome, σ = 1 Mb, 20 cores — 25k markers 6.7 s → 0.21 s; 50k markers 27.4 s → 0.69 s (40×); 200k markers (150 bp spacing) 10.2 s. `SmoothedStats` 560 → 280 bytes. Output identical to the golden files except CompositeZ, which moved by ≤ 1e-6 because it is now summed from unrounded Z-scores. `BulkFlags` now matches modes exactly (same result for all 10 modes). Still quadratic in marker density: the grid/binning change is open, and it alters results. Per-chromosome streaming is also still open.

### 1.1 Kernel cost (HIGH for dense marker sets)
- [x] **Compute weights once per query marker, not once per statistic.** `smooth1` runs a binary search and an `exp` for every (marker, neighbour) pair *for each of up to 7 statistics*. Computing the weights once and accumulating all statistics in the same inner loop makes it about 7× faster with identical results.
- [x] **Use two moving pointers instead of two `sort.Search` calls per marker.** Positions are sorted, so the window edges only move forward.
- [x] **Run chromosomes (or chunks of query markers) in parallel.** It currently runs on a single thread.
- [x] **Possibly change the algorithm** — DONE: grid evaluation (--step-size, ≤ σ/20) with interpolation, exact in sparse intervals; 2M markers 114 s → 4 s. (this would change results slightly, so it belongs in Part 2 if adopted):
  - bin markers at about σ/20 (50 kb), then convolve the bins (O(n + bins·kernel)); or
  - evaluate only on a regular grid. `--step-size` exists (default 100 kb) but is **never used** (`StepSize` is set and never read).
  - Either removes the quadratic cost.

### 1.2 Memory
- [x] **Remove the unused field:** `SmoothedStats` carries an unexported `thresholds Thresholds` field (280 B) that nothing reads. It makes up half of the 560-byte struct: 2M sites × 560 B = 1.1 GB.
- [x] **Avoid copying structs:** `extractFloat` passes `BSAstats` (224 B) by value to a closure once per statistic per site. Read the fields directly, or let the single pass from 1.1 read each site once.
- [x] (partly) Smoothing/normalising now process one chromosome at a time and normalise uses one buffer per field (≤ 3 at once): 2M-marker peak heap 1.5 GB → 1.18 GB. Full raw→output streaming is not possible while Z is genome-wide and later steps need all sites. Was: stream one chromosome at a time from raw stats to smoothing to writing (gemini 3.6). This goes with the in-memory pipeline that the genome-wide null (B1) needs.

### 1.3 TSV writer (`writeSmoothedTSV`)
- [x] **Return errors:** write, flush and close errors are ignored (deferred `Flush`/`Close`, unchecked `Fprintln`), so a full disk leaves a truncated file with exit status 0. This is the same fix as in raw stats.
- [x] **Faster writer:** `Sprintf`/`Join` per field. Use the reusable-buffer approach from `appendRawRow`.

### 1.4 Code health
- [x] **`BulkFlags` matches substrings** (`strings.Contains(bsaType, "hb")`). It is used throughout `stats/` and `plots/`. Derive the layout once from the sample roles, as raw stats now does.
- [x] **`validateCompositeZStatCount` panics at runtime** to check two constants against each other. Replace it with a unit test. The comment about an "import cycle" is wrong (both are in package `stats`), and "Fix 6" refers to nothing.
- [x] **Stale comments:**
  - `OneBulkP0` "taken from last contributor": it is the site's own value, which now differs by polarity.
  - "Mark as invalid by leaving POS = 0": the code uses POS = −1.
  - The `CompositeZ` field comment "always ≥ 0" is false, because Z_G and Z_LOD are negative below the median.
- [x] **Unnecessary rounding:** Z-scores, CompositeZ and MaxAbsZ are rounded to 6 dp in memory. As in raw stats, keep full precision and round only in the TSV.
- [x] **Report dropped markers per chromosome**, not just the total. A whole contig can disappear without notice.

## Part 2: biology and statistics (results change)

### 2.1 Depth weighting over-weights high-depth markers (MEDIUM) — DONE
The weight is proportional to depth, which is right only if read sampling is the only noise. SI variance is the pool term plus the read term, about p(1−p)·[1/(c·n) + 1/d], where c is the population variance scale and n the bulk size. Once depth exceeds about 2n, pool sampling dominates and markers deserve roughly equal weight.
- **Example (F2, 20 plants):** a 200x marker gets 5× the weight of a 40x marker, but its true precision is only 1.7× higher.
- **Why it matters:** high-depth markers are where collapsed repeats and paralogue mismapping concentrate, which are the least trustworthy markers.
- [x] Use inverse-variance weights: w = 1 / Σ_bulks [1/(c·n_b) + 1/d_b], with `PopulationVarianceScale` and the bulk sizes already in the config. For one-bulk modes, use the single bulk term.

### 2.2 Depth-dependent baselines of G, LOD and the Bayes factor (MEDIUM)
With no QTL, G grows with depth: about 1 at low depth, rising as pool variance becomes visible to deep sequencing (the mean G in the 60x one-bulk simulation was about 2.5). The smoothed G therefore follows local coverage, and regions with higher coverage look like weak signal. The threshold is looked up at the **central** marker's depth, while the smoothed value averages its neighbours' depths.
- [ ] At minimum, look up thresholds at the kernel-weighted mean depth (or the depth pair) of the window rather than the central marker's.
- [ ] Better: compare smoothed statistics with a null that is smoothed the same way (claude_review B1). Or use Magwene's G′ with a null fitted to the smoothed values.

### 2.3 Thresholds come from single SNPs but are applied to smoothed tracks (HIGH; claude_review B1, grok #3)
This is not a smoothing bug, but smoothing is where the mismatch arises.
- **Read noise:** smoothing averages away read-sampling noise across about 6 Mb of markers.
- **Pool noise:** it barely reduces pool-sampling noise, because linked markers share the same plants.
- **Result:** the single-SNP null is too wide wherever read noise dominates (low depth), so the scan loses power there. Calibration also changes with depth.
- [ ] Fix with the chromosome-scale null run through this same smoothing code (B1). Refactor `gaussianSmooth` so it is a pure function on in-memory data that the null can call.

### 2.4 The minimum-marker rule doesn't catch sparse regions (MEDIUM) — DONE
"At least 5 markers within ±3σ" counts markers up to 3 Mb away whose Gaussian weight is about 0.01. Near-zero-weight markers still count, so in a sparse region (around centromeres, introgressions, poorly called regions) the smoothed value is in effect the central marker alone. That gives noisy spikes that can cross thresholds.
- [x] Use the effective number of markers, n_eff = (Σw)² / Σw², with a threshold (e.g. ≥ 10), and write n_eff to the TSV.
- [ ] Optionally flag or drop markers whose nearest neighbour is more than σ away.

### 2.5 Chromosome ends (LOW–MEDIUM) — DONE (local linear near ends for SI/ΔSI/AFDev; loess for ED⁴)
Near a chromosome end the kernel is one-sided: there are fewer markers, the variance is higher, and the estimate is biased toward the interior. A QTL near a telomere is both flattened and noisier.
- [ ] Report n_eff (from 2.4) so this is visible. Consider a local-linear kernel (the Nadaraya–Watson edge correction) for SI and ΔSI.

### 2.6 Genome-wide robust Z (MEDIUM; gemini 2.4)
- **Skewed statistics:** G, LOD, ED⁴ and the Bayes factor are right-skewed and bounded at 0. (median, MAD) assumes symmetry, so their Z-scores are not normal deviates and their tails are uneven: the upper tail is heavier than the Z suggests.
- **Large or distorted regions shift the reference:** a large QTL region, or several, moves the median and MAD and so every other Z. In crosses with **segregation distortion** (common in plant crosses), the distorted regions shift AFDev in one-bulk modes exactly as a QTL would. They also shift HighSI/LowSI, though not ΔSI.
- **Autocorrelation:** the smoothed values are strongly autocorrelated, so the genome holds far fewer independent values than markers. The median and MAD are dominated by a few long blocks.
- [ ] Use the Z columns for plots only, and call QTLs on the smoothed values in their natural units against the null (claude_review B1 recommends this). If Z-scores stay, transform skewed statistics first (e.g. log G, as in Magwene's G′), and estimate centre and spread with QTL regions excluded (Hampel-trimmed, as in G′).

### 2.7 Meaningless SI tracks in `2b` without `--ref-parent` (LOW–MEDIUM) — DONE
SI there is the REF fraction, so smoothed HighSI and LowSI hover around 0.5 even at a QTL (the polarity cancels). They are still plotted with thresholds and enter MaxAbsZ.
- [x] In that mode, leave Sm/Z HighSI and LowSI out of MaxAbsZ, and label or drop the plots. (Labelled "REF fraction (polarity unknown: not a QTL track)", no thresholds, removed from the Z overlay.)

### 2.8 Split multi-allelic sites are counted twice (LOW) — DONE in smoothing (one kernel neighbour per position: mean weight, weighted-mean value, NEff counts it once). Still open in the HMM, whose emission multiplies the reads of every view at a position.
After splitting, two ALTs at one position become two markers at the same position, and both enter the kernel with full weight. They come from the same reads and plants.
- [ ] Keep one view per position for smoothing (e.g. the one with the higher depth), or halve the weights.

### 2.9 What gets smoothed (LOW; grok #4) — DECIDED
**Decided:** per-site statistics are smoothed, as in the literature. G and G1 use a tricube kernel (G′); ED⁴ uses loess; SI, ΔSI and AFDev keep the inverse-variance-weighted Gaussian kernel. See `smoothing_algorithms.md`.
Smoothing per-site G, ED⁴ and LOD directly is what the published methods do: G′ is a kernel-smoothed G (Magwene 2011), and MMAPPR fits a loess to ED⁴ (Hill 2013). This is acceptable as it stands. The alternative of computing the statistics from smoothed allele counts is less noisy but is a different method.
- [ ] Decide which to follow and document it. No change is needed if the published approach is kept.

### 2.10 Kernel width and the `--window-size` name (DOCS)
σ = WindowSize/2 and the cutoff is ±3σ, so the default "2 Mb window" uses markers up to 3 Mb either side (6 Mb in all). QTLseqr uses a 1 Mb tricube window.
- [ ] Document what `--window-size` means, or rename it (e.g. `--kernel-sigma`). Remove or wire up `--step-size`.
- [ ] Give guidance: the width should match the expected size of the linkage block, roughly the genome length divided by the expected number of crossovers per bulk. Too wide a kernel merges neighbouring QTLs and lowers peak height.

---

## Deferred (CompositeZ review)
- ~~LOD is G divided by a constant, so CompositeZ counts the same evidence twice.~~ Fixed: LOD replaced by the pool-HMM LOD.
- CompositeZ sums Z-scores of skewed statistics, which are negative below the median, so it is not a Stouffer Z, and the `√k` scaling assumes independent inputs.
- ~~The CompositeZ null uses single SNPs, one median depth and no smoothing.~~ Fixed: thresholds now come from simulated no-QTL genomes run through the same scoring code (pointwise p99/p95; calibration checked in `TestHMMLODCalibration`).

## Suggested order
1. Part 0 (known-value tests, golden test, benchmark).
2. 1.1 (single pass, two pointers, parallel) and 1.2 (drop the `thresholds` field), then 1.3 and 1.4.
3. 2.4 (n_eff) and 2.7 (SI tracks in `2b`): small and local.
4. 2.1 (inverse-variance weights) and 2.2 (threshold depth from the window).
5. 2.3 and 2.6 together with the genome-wide null (B1) and the CompositeZ review.

## Done since the review
- Short contigs (< 6σ span) keep the local average; the edge local linear fit needs two full kernel widths.
- Null-simulation speed: HMM log read probabilities cached per shift (byte-identical); grid-mode kernel neighbours aggregated into σ/100 bins at weighted centroids. 2M markers: 2.54 → 0.98 s per null replicate.
