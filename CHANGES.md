# Changes: review fixes (2026-11-17)

Fixes for issues found in a code review of GoBSAseq, each verified on simulated
data (see `Data/`). Build, `go vet` and all tests pass.

## Summary

| # | Problem | Fix |
|---|---------|-----|
| 1 | Tabix index (`.tbi`) of the hard-filtered VCF was unusable | New index writer, built from the finished file |
| 2 | Threshold level fell back from P99 to P95 per chromosome | QTLs called at P99 on every chromosome |
| 3 | Thresholds came from a null that didn't match the tested statistic | Genome-wide thresholds from simulated null genomes |
| 4 | Thresholds changed between identical runs | Deterministic seeding, new `--seed` flag |
| 5 | Run stopped at a `[y/N]` prompt (fails in batch/HPC jobs) | Gene-space check is non-interactive |
| 6 | Hard filtering slow and memory-hungry on multisample VCFs | Unused sample columns dropped before parsing; compact records; parallel BGZF decompression |
| 7 | Smoothing slow (single core, weights recomputed per statistic) | One shared weight pass on all cores; identical output |
| 8 | Plotting slow, memory-hungry, and HTML too large to open | Peak-preserving downsampling, 2-point threshold lines, numeric axis |
| 9 | A strong QTL was reported as the whole chromosome | LOD-drop support interval around the peak, reported first in every QTL table |
| 10 | Null-genome simulation too slow for plant-sized genomes | Smoothing per 10 kb bin instead of per marker; tabulated read sampling and log-gamma |
| 11 | BRM called QTLs in 13–40% of genomes with none, more with deeper sequencing | Drift-only σ and a genome-wide multiplier from the null simulation; `--brm-alpha` is now genome-wide (default 0.01) |
| 12 | CompositeZ found only 4 of 10 weak QTLs (G/LOD found 8) | Z-scores use a fixed median and MAD from simulated null genomes instead of each genome's own |
| 13 | `--min-qtl-width` / `--min-qtl-markers` did nothing; BRM one-bulk fallback bug | Both flags now filter every method's regions; fallback removed |
| 14 | Support-interval drop (1.5 LOD) uncalibrated; too wide at low depth | Drop of 6.5 null LOD scales, calibrated to 95% coverage on 216 simulated datasets |
| 15 | Peak memory ~2.5 GB at 900k markers (duplicated per-marker data) | One genome-wide threshold value, unused field removed, chromosomes scanned in place: −40% |
| 16 | `--min-SOR-SNP` / `--min-FS-SNP` were upper limits | Renamed `--max-SOR-SNP` / `--max-FS-SNP`; old names kept as deprecated aliases |
| 17 | README out of date | Rewritten to match the current app |
| 18 | genome-whisperer update broke the build | Adapted to the new `genespace.GeneSpace(Options)` and `variants.RunPipeline(Options)`; gene descriptions optional |

---

## 1. Tabix index

**Problem.** Region queries on `GoBSAseq.<type>.hardfiltered.vcf.gz` returned
almost nothing: `bcftools view -r chr3:10000000-12000000` gave 0 records where a
full scan found 166. There were two causes:

- The chunk offsets were taken from the compressed byte count while the BGZF
  writer was still compressing in the background, so they lagged behind the
  real positions.
- The index library itself is broken. In `github.com/biogo/hts` v1.4.5 (the
  latest release), `tabix.Index.Add` never records new contig names, so every
  record became its own contig. A 5-chromosome file got an index with 12,315
  contigs.

**Fix.** New `filter/tabix.go` writes the `.tbi` itself. It reads the finished,
closed file back to get exact offsets, then writes standard UCSC bins, merged
chunks and a 16 kb linear index. The biogo tabix package is no longer used.

**Verified.** 60 random region queries and 10 whole-chromosome queries return
exactly the records a full scan finds. New tests: `filter/tabix_test.go`.

**Files.** `filter/tabix.go` (new), `filter/tabix_test.go` (new), `filter/filter.go`.

## 2. One threshold level for the whole genome

**Problem.** For each statistic, each chromosome reported only its P99 regions
if it had any, and otherwise its P95 regions. Whether a region was called
therefore depended on what else was on its chromosome. In a weak-QTL simulation
this let two false QTLs through at P95 (chr3 and chr4).

**Fix.** Regions are called at P99 on every chromosome (`findPeaksAtP99` in
`stats/detect.go`). P95 values are still written to the thresholds TSV.

**Files.** `stats/detect.go`.

## 3. Genome-wide thresholds from simulated null genomes

**Problem.** QTLs are called on Gaussian-smoothed, genome-normalised
statistics, but thresholds came from simulating single, unlinked SNPs. That
null leaves out the main source of variation in smoothed data: genetic drift.
A bulk holds a limited number of plants, so long parental segments are over- or
under-represented by chance, and linked markers move together. The CompositeZ
P99 threshold swung between 10.3 and 16.9 depending on the dataset, and null
chromosomes could pass it.

**Fix.** New `stats/nullsim.go`. For each of `--rep` simulated genomes with no
QTL:

1. Every bulk plant is simulated through the cross (`-m`: F2/F3/F4/RIL by
   selfing, BCn by backcrossing), chromosome by chromosome, with Haldane
   crossovers.
2. At every real marker, reads are drawn from the simulated pool frequency
   using that marker's observed depth. Real marker positions and depths are
   kept.
3. The pipeline's own statistics, smoothing, robust-Z normalisation and
   CompositeZ formula are applied. The CompositeZ formula is now a single
   shared function (`compositeZ` in `stats/smoothing.go`).
4. The genome-wide maximum (and minimum, for signed tracks) of every track is
   recorded.

P99 is the 99th percentile of those maxima: only 1% of genomes with no QTL
produce a call anywhere, so it controls the **genome-wide** false-positive
rate. Signed tracks (SI, ΔSI, AFDev) split it across both tails.

**Assumption.** Every chromosome is taken to be 100 cM long, about one
crossover per meiosis, which is typical for plants. It is the constant
`morgansPerChromosome` in `stats/nullsim.go` and is printed in the run log.
There is no flag for it yet.

**Verified.**
- Calibration: across 20 simulated genomes with no QTL, each run's observed
  maximum CompositeZ fell between the 6th and 91st percentile of its simulated
  null (mean rank 0.43). The simulated null matches reality, if anything
  slightly conservatively.
- Strong QTL: called on every track (CompositeZ 160 vs threshold 27.6).
- Weak QTL: no longer called. In bulks of 50 plants, chance drift alone
  reached CompositeZ 18.2 in one genome with no QTL, higher than the weak
  QTL's 17.6. Weak QTLs need larger bulks.
- Runtime: 150k markers with 1,000 null genomes took 104 s on 16 cores
  (1.4 GB peak memory).
- New tests in `stats/regression_test.go`: recombination follows Haldane's
  formula, each cross type gives its expected allele frequency, and thresholds
  hold their error rate on independent null genomes.

**Behaviour changes.**
- `CalculateThresholds` now takes the raw statistics as well:
  `CalculateThresholds(cfg, bsaType, raw, smoothed)`.
- Thresholds are the same for every marker (genome-wide); the format of
  `*.thresholds.tsv` is unchanged.
- New output: `stats/GoBSAseq.<type>.null_genome_extremes.tsv`, the per-genome
  maxima behind the thresholds.
- `--rep` now means the number of simulated null genomes (default 1000).
- `--bulk-sizes` now directly sets how much drift the null includes, so it must
  be accurate.
- The old per-depth single-SNP simulation (~650 lines) was removed, which also
  dropped the `golang.org/x/sync` dependency (`go mod tidy`).

**Files.** `stats/nullsim.go` (new), `stats/thresholds.go`,
`stats/smoothing.go`, `stats/regression_test.go`, `run/run.go`, `go.mod`,
`go.sum`.

## 4. Reproducible thresholds

**Problem.** Two runs on identical input gave different thresholds (Gstat P99
11.79 vs 10.54), because random seeds were handed out in whatever order worker
goroutines happened to run.

**Fix.** Each simulated genome's seed is derived from a base seed plus the
genome's index, independent of scheduling. New flag `--seed` (default 1). The
same seed and input give byte-identical thresholds; a different seed gives a
different, equally valid set.

**Files.** `stats/thresholds.go`, `stats/nullsim.go`, `utils/config.go`,
`cmd/root.go`.

## 5. No gene-space prompt

**Problem.** If gene-space parameters were missing, the run asked
`Continue without gene space analysis? [y/N]` and exited with an error when
there was no keyboard input (batch, HPC and CI jobs).

**Fix.** The check is now non-interactive:
- None of `--snpEffDB`, `--gene-descriptions`, `--prg`, `--gff` given: a note
  is logged and gene-space analysis is skipped.
- Only some given: the run stops immediately and lists which are missing, since
  that is almost certainly a configuration mistake.

**Unchanged:** choosing the high/low parent and bulk interactively. If `-P`/`-B`
are left out, or a name isn't found in the VCF, the sample list is still shown
and you pick by number. Batch runs should pass the names on the command line.

**Files.** `run/run.go`.

## 6. Hard-filtering speed and memory

**Problem.** Filtering time and memory grew with the number of samples in the
VCF, even though at most four are used. The VCF reader (vcfgo) fully parsed
every FORMAT field of every sample into a per-sample map, on one thread, and
the filter then threw all but four away. On top of that, every passing record
was kept in memory as a full parsed `vcfgo.Variant` for the rest of the run.

**Fixes.**
1. **Unused samples dropped before parsing.** The filter re-reads the VCF
   through a lightweight reader that keeps only the selected sample columns
   (`filter/vcfreader.go`), so vcfgo only ever parses those. Sample role
   indices are remapped to the subset automatically. This also makes the
   filter tolerate Windows (CRLF) line endings.
2. **Compact records.** Passing sites are stored as `FilteredVariant` with
   only what the statistics need (position, alleles, and per sample: GT, REF
   depth and target-ALT depth), instead of the full parsed record.
   `RawStats` reads these directly.
3. **Parallel decompression.** BGZF input (the standard output of bcftools,
   GATK and DeepVariant) is decompressed on all cores. Benchmarked on
   realistic data at 1,020 MB/s against 202 MB/s for Go's standard gzip
   reader. Plain gzip and uncompressed VCFs fall back automatically.

**Verified.** Every output (filtered and rejected VCFs, index, raw, smoothed,
thresholds, BRM) is byte-identical to before on four inputs: 200 samples
BGZF, 200 samples plain gzip, 4 samples uncompressed, and
`Data/sim_F2_strongQTL.vcf.gz`. New test `TestSubsetFastPathMatchesFullParse`
compares the fast path with the original full parse on a multisample VCF.

Filter step only:

| Input | Before | After |
|---|---|---|
| 150k records, 200 samples (BGZF) | 10.3 s, 1,281 MB | 1.2 s, 146 MB |
| 20k records, 200 samples, realistic data (BGZF) | 2.7 s | 0.3 s |
| 150k records, 4 samples (uncompressed) | 0.7 s, 430 MB | 0.6 s, 125 MB |

**API change.** `filter.FilteredVariant` no longer holds `*vcfgo.Variant`;
build one with `filter.NewFilteredVariant(v, targetAlt)`.

**Files.** `filter/vcfreader.go` (new), `filter/filter.go`,
`filter/filter_test.go`, `stats/stats.go`.

---

## 7. Smoothing speed

**Problem.** `gaussianSmooth` ran one pass per statistic (nine for two-bulk
runs). Each pass recomputed the same Gaussian × depth weight, an `exp()` call,
for every marker pair within 3σ, all on one core: about 8 billion `exp()`
calls for 150k markers, almost 90% of them repeats.

**Fix.** New `smoothAll` (`stats/smoothing.go`) computes each weight once and
applies it to all statistics in a single pass, with marker positions split
across all cores. Each statistic still accumulates in the same order, so
output is **bit-for-bit identical**.

**Verified.** `smoothed_and_normalised.tsv` and every downstream output is
byte-identical to before.

**Files.** `stats/smoothing.go`.

## 8. Plotting speed, memory and file size

**Problem.** Every chart drew every marker in every series, about 25 series
per chromosome, on a categorical x-axis with a text label per marker. Threshold
lines repeated a constant value at every marker, and each chromosome's data
was copied before plotting. On 150k markers the plotting step peaked at 1.3 GB
and wrote 144 MB of HTML (`individual_plots.html` alone was 107 MB), which a
browser struggles to open.

**Fixes** (`plots/plots.go`).
- **Peak-preserving downsampling.** Each series keeps only the minimum and
  maximum point of every bin of 1% of the smoothing σ (10 kb at the default
  2 Mb window). Curves are smoothed over about 1 Mb, so no visible detail is
  lost, and every peak and dip is kept exactly.
- **Threshold and reference lines** are 2-point series instead of one point
  per marker.
- **Numeric position axis** (`type: value`) instead of a text label per
  marker. BRM shading now uses exact block coordinates.
- **No per-chromosome copies** of the smoothed and threshold records; charts
  index into the existing slices.
- Tooltips read `[position, value]` points. The individual-statistic tooltip
  now has the p99/p95 values embedded, since the threshold lines are no longer
  sampled at every position.

**Verified.**
- All three pages rendered in a real browser (headless Edge). Axes, threshold
  lines, BRM shading and hover tooltips (including ★ p99 at the QTL peak) all
  work.
- Across all 80 plotted series (every statistic on every chromosome, all
  three pages), the drawn maximum and minimum match the full data exactly,
  with 24% of the points drawn on this 1 kb-spaced data. Sparser marker sets
  are reduced less; denser ones more.
- TSV and Excel outputs are unchanged and still contain every marker.

**Files.** `plots/plots.go`.

## Speed and memory summary (steps 6–8)

150k markers, 4 samples, 16 cores (`--rep 10`, so threshold simulation is
excluded):

| Step | Before | After |
|---|---|---|
| Filtering | 0.7 s, 436 MB | 0.6 s, 130 MB |
| Smoothing | 30.3 s, 636 MB | 1.6 s, 256 MB |
| Plotting | 7.1 s, 1,279 MB | 0.9 s, 306 MB |
| Whole run | 39.9 s | 4.9 s |
| HTML output | 144 MB | 16.8 MB |

With 200 samples in the VCF, filtering went from 10.3 s and 1,281 MB to 1.2 s
and 146 MB (section 6).

---

## 9. QTL support intervals (localisation)

**Problem.** A QTL was reported as the stretch where the statistic is above the
threshold. With a strong QTL in an F2 that is the whole chromosome: after one
round of meiosis even markers 50 cM away are still strongly linked to the QTL
(expected ΔSI = e^(−2d) for complete selection, d in Morgans, so 0.37 at the
ends of a 1 M chromosome with the QTL in the middle). On
`Data/sim_F2_strongQTL.vcf.gz` every output reported chr2:1,000–29,995,075.
The original per-marker thresholds behaved the same way (all 2,472 chr2
markers above P99), so this was not caused by the switch to genome-wide
thresholds (section 3). The threshold answers *whether* a chromosome carries a
QTL; it cannot say *where*.

**Fix.** The threshold still decides which regions are QTLs. Each region now
also gets a **support interval**: starting from the highest smoothed LOD in the
region, it extends in both directions while LOD stays within 1.5 of the peak
(the classical 1.5-LOD drop; `SupportLODDrop` in `stats/support.go`). One-bulk
designs use the one-bulk LOD.

- **QTL tables** (`final`, `compositez`, `individual` `.qtl.xlsx`): the first
  columns are now `Support_Start`, `Support_End`, `Support_PeakPos`,
  `Support_PeakLOD`; the above-threshold region follows as
  `Region_Start`/`Region_End` (previously `Start`/`End`). When regions are
  merged, support intervals are combined (span of both, higher peak).
- **BRM blocks TSV**: `SUPPORT_START`, `SUPPORT_END`, `SUPPORT_PEAK_POS`,
  `SUPPORT_PEAK_LOD` appended, so existing column positions are unchanged.
- **Plots**: support intervals are shaded red and labelled on every chart,
  on top of the amber BRM blocks.
- **Gene-space analysis** now searches the support interval instead of the
  whole above-threshold region.
- The run prints each QTL's support interval at step 8.

**Result** on `sim_F2_strongQTL` (true QTL at chr2:15.0 Mb):

| Interval | chr2 | Width |
|---|---|---|
| Above P99 threshold (before) | 0.00–30.00 Mb | 30 Mb |
| 1.5-LOD support interval (now) | 14.60–16.29 Mb, peak 15.44 Mb | 1.7 Mb |

The interval is the same whether bulk sizes are given correctly (`-S 50,50`) or
left at the default (`20,20`); only the thresholds change. Weak-QTL and null
datasets still produce no calls. Statistics and thresholds outputs are
byte-identical to before.

**Not yet calibrated.** 1.5 comes from classical interval mapping, where each
marker is tested once. Here LOD is smoothed across neighbouring markers, so
the right drop for ~95% coverage of the true position is unknown. Next step:
simulate a few hundred datasets with `Data/simulate_bsa.py` (random QTL
position and strength) and choose the drop that gives ~95% coverage.

**Files.** `stats/support.go` (new), `stats/support_test.go` (new),
`stats/detect.go`, `stats/brm.go`, `stats/genespace.go`, `plots/plots.go`,
`run/run.go`.

---

## 10. Null-genome simulation speed

**Problem.** The threshold simulation (section 3) smooths every simulated
genome at every marker, using a 601-bin kernel when σ = 1 Mb (the default
2 Mb window). Cost grew linearly with marker count: 0.6 s per simulated genome
at 150k markers, 3.8 s for a tomato-sized genome (900k markers) and 12.5 s at
3M markers, times 1,000 genomes. The full tomato-sized run took 7.3 minutes.

**Fix** (`stats/nullsim.go`):

1. **Smooth each 10 kb bin once.** About 10 markers share each bin and got
   bit-identical smoothed values, so each occupied bin is now computed once.
   Robust-Z normalisation weights each bin by its marker count
   (`robustZScratch`), giving exactly the same median, MAD and Z as before.
   Per-worker buffers are now sized by bins, not markers.
2. **Tabulated beta-binomial Bayes factor.** It always takes log-gamma of
   whole or half-integer arguments, so these come from a table. The values are
   bit-identical.
3. **Tabulated read sampling.** Pool frequencies can only be c/2n, and depths
   are the observed ones, so each Binomial(depth, c/2n) is tabulated once
   (depth ≤ 500; deeper markers use the old sampler). Each draw is one uniform
   number and a binary search, replacing gonum's rejection sampler, which was
   half of the remaining time. This changes the random stream, so thresholds
   change within Monte Carlo noise.

**Checks.**

- Steps 1–2 alone: all 1,000 null-genome extremes byte-identical to before
  (150k markers).
- Step 3: 3,000 simulated genomes with each sampler. For every statistic the
  null distributions are indistinguishable (Kolmogorov–Smirnov D 0.018–0.029,
  5% critical value 0.035), and P99s agree within 2%.
- New tests: weighted robust Z matches the expanded version bit for bit,
  tabulated BF matches `betaBinomialLogBF`, and sampler mean and variance
  match the binomial.

**Result.** Per simulated genome on one core (`BenchmarkNullGenome`):

| Genome | Before | After |
|---|---|---|
| 150k markers | 609 ms | 76 ms |
| Cucumber-sized (7 × 50 Mb, 350k markers) | 1.47 s | 0.18 s |
| Tomato-sized (12 × 75 Mb, 900k markers) | 3.81 s | 0.49 s |
| Tomato-sized, 3M markers | 12.5 s | 0.85 s |

Full default runs (`--rep 1000`, 16 cores):

| Input | Original code | Before step 10 | Now |
|---|---|---|---|
| 150k markers | 46 s, 2.1 GB | 70 s, 0.56 GB | 20 s, 0.45 GB (steps 1–2 only) |
| Cucumber-sized, 350k markers | 122 s, 9.4 GB | 164 s, 1.26 GB | 32 s, 0.92 GB |
| Tomato-sized, 900k markers | not run (projected > 20 GB) | 435 s, 3.1 GB | 84 s, 2.7 GB |

Support intervals and QTL calls are unchanged on these files.

**Files.** `stats/nullsim.go`, `stats/regression_test.go`,
`stats/nullsim_bench_test.go` (new: benchmark and tests).

---

## 11. BRM accuracy

BRM is kept with its formula (a block is significant where |ΔSI| ≥ u·σ;
one bulk: |SI − p0| ≥ u·σ), but four problems made its error rate wrong and
unstable.

**Problems.** Measured on simulated F2 genomes with no QTL (bulks of 50, 15
genomes per setting):

1. **Per-position significance level.** u was `qnorm(1 − α/2)` with α = 0.05
   at each position, so the chance of a false block anywhere grew with genome
   size.
2. **Single-marker noise model for a smoothed value.** σ added read-sampling
   noise, af(1−af)/depth per bulk, as for one marker. The tested ΔSI is
   smoothed over hundreds of markers, which removes almost all read noise. The
   assumed SD (0.119 at depth 60) was about twice the real one (0.062).
3. **Error rate depended on depth.** The oversized read term shrank as depth
   rose, so false calls went from 2/15 genomes at 30× to 3/15 at 60× and
   6/15 at 150× (4/15 with 12 chromosomes at 60×).
4. **Threshold jumped between markers.** σ used each marker's raw depth, so
   the threshold zig-zagged (0.213–0.254 at 60×). False blocks were 0.1–0.2 Mb
   long, far narrower than the 1 Mb smoothing scale.

**Fix** (`stats/brm.go`, `stats/nullsim.go`, `stats/thresholds.go`):

- σ is pool sampling (drift) only:
  σ² = af(1−af)·[1/(V·n1) + 1/(V·n2)], or p0(1−p0)/(V·n) for one bulk, with V
  the population variance scale (unchanged). It no longer uses depth, so it is
  smooth along the chromosome.
- u is genome-wide. Each simulated null genome records its largest |ΔSI|/σ
  (same σ, same smoothing as the real data), and u is the (1 − α) quantile.
  Residual read noise at the real depths is included through the simulation,
  so the error rate no longer depends on depth or genome size.
- `--brm-alpha` now means "chance that a genome with no QTL gets any BRM
  block". The default changed from 0.05 (per position) to **0.01**, matching
  the P99 used for the other statistics. A note is printed if α × `--rep` < 10.
- New outputs: `BRM_U` column in the thresholds TSV, `BRM_maxAbsZ` in the
  null-extremes TSV. The unused `inverseNormalCDF` was removed.

**Result** (same simulated datasets as before):

| Test | Before | After |
|---|---|---|
| No QTL, 30× | 2 / 15 genomes with a call | 0 / 15 |
| No QTL, 60× | 3 / 15 | 0 / 15 |
| No QTL, 150× | 6 / 15 | 0 / 15 |
| No QTL, 12 chromosomes, 60× | 4 / 15 | 0 / 15 |
| Strong QTL: found / wrong-chromosome calls | 20 / 20, 2 | 20 / 20, 0 |
| Weak QTL: found / wrong-chromosome calls | 8 / 10, 1 | 8 / 10, 0 |
| Cucumber-sized test (350k markers): false chr4 block | yes | no (chr3 QTL still found) |

The BRM cutoff now adapts to the data. With bulks of 50 it is ΔSI ≈ 0.33 at
30× (more read noise left after smoothing), 0.25 at 60× and 150×, and 0.27
with 12 chromosomes (more chances for drift). Tests:
`TestNullThresholdCalibration` now also checks that BRM's simulated cutoff
gives the stated error rate on independent null genomes; `TestBRMSigma`
checks σ.

**Not changed:** BRM blocks still join the final table without a minimum
width or marker count (problem 5 of the BRM review), and the one-bulk
`SmAFDev == 0` fallback (problem 6) is still there.

---

## 12. CompositeZ: Z-scores against the simulated null

**Problem.** Each statistic was turned into a Z-score with the genome's *own*
median and MAD, in the real data and in every simulated genome. Two effects
cost power:

- A QTL inflates its own genome's spread (with 5 chromosomes, a whole-chromosome
  QTL is 20% of the markers), which shrinks its own Z-scores. In one weak-QTL
  run G was 2.2× its threshold but CompositeZ only 0.6× its own.
- The MAD varies from genome to genome, so the per-genome maximum of
  CompositeZ had a heavy tail (P99 3.3× the median), giving high, noisy
  thresholds.

**Fix** (`stats/nullsim.go`, `stats/smoothing.go`, `stats/thresholds.go`).
A separate batch of 50 simulated null genomes (own seeds, about 5% extra
run time) gives each statistic a fixed median and scale (1.4826 × MAD,
weighted by markers per bin; very long genomes are thinned to 2M values per
statistic). The same transform, z = (x − median) / scale, is applied to the
real data and to every threshold genome, so Z measures distance from what a
QTL-free genome produces. The values are written to
`GoBSAseq.<type>.null_scales.tsv`.

Because Z-scores now need the simulation, they are filled in the threshold
step, and the smoothed TSV is written there (same file name and columns).
The per-genome robust Z code (`normalise`, `robustZScores`, `medianOf`) was
removed. The robust-Z overlay plot now shows Z against the simulated null.

**Result** (same 60 simulated datasets as the earlier benchmark):

| Test | Before | After |
|---|---|---|
| Weak QTL found by CompositeZ | 4 / 10 | **8 / 10** (same as G/LOD) |
| No QTL: genomes with a CompositeZ call | 0 / 30 | 0 / 30 |
| Strong QTL found / wrong-chromosome calls | 20 / 20, 0 | 20 / 20, 0 |
| Null CompositeZ maximum, P99 ÷ median | 3.3 | 2.4 |
| `Data/` CompositeZ P99 (strong / weak / null) | 28.4 / 29.2 / 31.2 | 22.9 / 24.6 / 22.3 |

All 60 datasets now give 0 false calls from every method (CompositeZ,
individual statistics, BRM). Tests: `TestWeightedMedianMatchesExpanded`,
`TestEstimateScales` (ΔSI centred near 0, reproducible for a seed).

---

## 13. Minimum region size; BRM one-bulk fallback

**Problem.** `--min-qtl-width` and `--min-qtl-markers` were parsed but never
used, so a region of any width or marker count could become a QTL, including
BRM blocks entering the final table. Separately, one-bulk BRM fell back to
`SmHighSI − p0` whenever `SmAFDev` was exactly 0; in a low-bulk-only run
`SmHighSI` is 0, so that marker got a deviation of −p0 and was called
significant.

**Fix.**

- New `RegionFilter` (`stats/regionfilter.go`) applied to the final regions
  of all three methods after merging: individual statistics
  (`DetectIndividualQTLs`), CompositeZ (`DetectCompositeZQTLs`) and BRM
  (`RunBRM`). A region must be at least `--min-qtl-width` wide (default 0 =
  smoothing σ, i.e. `--window-size`/2; negative disables) and contain at
  least `--min-qtl-markers` markers (default 5; 0 disables). The run prints
  how many regions each method dropped. Thresholds are calibrated without the
  filter, so it can only make calls more conservative.
- One-bulk BRM uses `SmAFDev` directly (no fallback).
- BRM blocks are now sorted by chromosome and position (they followed Go's
  random map order).

**Result.** On the 60-dataset benchmark nothing was dropped and results are
unchanged (strong 20/20, weak 8/10, 0 false calls in 30 null genomes), so
the filter costs no power there; it is a safeguard against narrow artefacts.
Test: `TestRegionFilter`.

---

## 14. Calibrated support intervals

**Problem.** Support intervals (section 9) used the classical 1.5-LOD drop,
which was never checked for smoothed BSA-seq LOD. Smoothed LOD grows with read
depth, so a fixed drop in LOD units is relatively wide at low depth and
narrow at high depth.

**Calibration.** 216 simulated F2 datasets (`Data/simulate_bsa.py`, which
gained a `moderate` strength): depth 30, 60 and 150×; bulks of 30, 50 and 100;
strong, moderate and weak QTLs at random positions; 8 replicates each; default
2 Mb window. In 197 the QTL was detected. For each, I measured the smallest
drop from the LOD peak at which the interval reaches the true position, in
three candidate units:

| Unit | Spread of the 95% drop across conditions (max ÷ min) |
|---|---|
| Raw LOD | 3.5× |
| Fraction of peak LOD | 13× |
| **Null LOD scale (Z_LOD, section 12)** | **1.8×** |

Null-scaled units carry over best, so the drop is now **6.5 × the null LOD
scale** (`SupportZDrop` in `stats/support.go`; the value in LOD units, about
0.3 at 30×, 0.6 at 60× and 1.2 at 150×, is printed at the threshold step and
stored in `Thresholds.SupportDrop`).

| Drop | Coverage of true QTL | Median width | 90th pct width |
|---|---|---|---|
| 1.5 LOD (before) | 99% | 3.8 Mb | 26.8 Mb |
| 6.5 null scales (now) | 95% (92–98% per depth, bulk size, strength) | 2.3 Mb | 5.9 Mb |

Weak QTLs: median width 16.7 → 4.9 Mb. Strong: 1.6 → 1.0 Mb.

**Held-out check** (the 60-dataset benchmark, different seeds): intervals
contained the true QTL for 26 of 28 detected QTLs (93%); median width strong
1.69 → 1.05 Mb, weak 13.1 → 6.4 Mb; no change to calls (0 false in 30 null
genomes). `Data/sim_F2_strongQTL`: 14.93–15.94 Mb (truth 15.0 Mb), was
14.60–16.29 Mb.

**Limits.** Calibrated for F2 with the default window and 1 Morgan per
chromosome; other crosses and window sizes are untested. A 95% interval
misses the true QTL about 1 time in 20 by design.

**Found along the way.** The default `--bulks-depth 40,40` removes every
marker where a bulk has fewer than 40 reads. With bulks sequenced at ~30×,
only about 80 of 15,000 markers passed and the localisation was poor; those
calibration runs were redone with `-b 10,10`. See Known issues.

---

## 15. Memory

**Problem.** Peak memory was about 3 KB per marker (2.5 GB for 900k markers).
Three things duplicated per-marker data: the genome-wide thresholds were
stored as one 304-byte `Thresholds` value per marker; every `SmoothedStats`
carried an unused 304-byte `thresholds` field; and QTL detection copied each
chromosome's records and thresholds before scanning them.

**Fix.** `CalculateThresholds` returns a single `Thresholds` (the thresholds
TSV still lists it per marker); the unused field is gone (`SmoothedStats`
568 → 264 bytes); `chromosomes()` (`stats/detect.go`) splits the sorted
input into per-chromosome sub-slices without copying (it copies only if
chromosomes are interleaved). `FindPeakIntersections` now takes the threshold
as a function of marker index.

**Result** (tomato-sized, 900k markers, full default run):

| | Before | After |
|---|---|---|
| Peak memory (RSS) | 2.48 GB | **1.49 GB** |
| Max heap, threshold step | 1.96 GB | 1.37 GB |
| Max heap, QTL detection | 2.03 GB | < 0.7 GB |
| Run time | 79 s | 77 s |

Every output TSV and every QTL call is byte-identical. Individual-statistic
and CompositeZ QTLs are now listed in chromosome order (they followed Go's
random map order). Test: `TestChromosomes`.

## 16. SOR / FS flag names

`--min-SOR-SNP` and `--min-FS-SNP` set upper limits (variants above them
fail). They are now `--max-SOR-SNP` and `--max-FS-SNP`, matching
`--max-SOR-INDEL` / `--max-FS-INDEL`. The old names still work, print a
deprecation notice and are hidden from `--help`; giving both an old and a new
name is an error. Defaults (3 and 60) are unchanged.

## 17. README

Rewritten to match the app: current flags (filtering, depth caps, light
filtering, region filter, seed, variant calling, config file), output files
and columns, pipeline steps, CompositeZ formula, BRM as Block Regression
Mapping (Huang et al. 2020), support intervals, test data, troubleshooting
(including the `--bulks-depth` default), project structure, and corrected
references (the Magwene 2011 citation was wrong; QTL-seq and BRM added).

---

## 18. genome-whisperer update

`github.com/gmaffy/genome-whisperer` was updated from
`v0.0.0-20260706090856-3d8bbe76817b` to `v0.0.0-20261118091955-61bdb2472be2`
(`go get …@latest`, `go mod tidy`). Two API changes broke the build:

- **`genespace.GeneSpace`** now takes an `Options` struct and returns
  `(rows, reportPath, error)`. `stats/genespace.go` fills it from the config
  and prints the number of genes scored and the report path per interval. The
  new options are exposed as `--trait` (default `disease-resistance`),
  `--min-concordance` (default 0.9, validated to (0, 1]) and
  `--promoter-length` (default 2000), using genome-whisperer's own defaults.
- **`variants.VariantCalling` was removed**; the BAM path in `run/run.go` now
  calls `variants.RunPipeline(variants.Options{…})`. `Species` and `RefVer`
  are required but only label output files outside genome-whisperer's
  data-directory mode, so they are set to the snpEff DB name (or `bsaseq`) and
  the reference file name, as the old call effectively did. `MinGQ` is 0:
  genome-whisperer's GQ filter would apply to bulks, where GQ is correlated
  with allele frequency; GoBSAseq applies `--min-GQ` to parents itself.

**Gene descriptions are optional.** Gene-space now needs `--snpEffDB`,
`--prg` and `--gff` (all or none; a partial set is an error), and
`--gene-descriptions` may be left out. The PRG file is still required:
genome-whisperer's annotation step (`CreateSuperVcf`) opens it
unconditionally.

**Checks.** Build, `go vet` and tests pass; a normal run's outputs are
byte-identical; the parameter checks behave as described. Gene-space and BAM
variant calling were not run end to end (snpEff, GATK and DeepVariant are not
installed in the test environment).

**README.** The read (FASTQ) input flags were documented as working, but a
read-based run only prints a message and stops; the README now says so.

---

## Also added

- `Data/`: three simulated F2 datasets (strong QTL, weak QTL, no QTL), the
  simulator that made them, and a README with the expected results.
- README: documents the new thresholds, `--seed`, and the gene-space rule.

## Known issues not addressed

- **One support interval per region.** Two linked QTLs on one chromosome fall
  in one above-threshold region and only the higher peak gets an interval.
- **Default `--bulks-depth 40,40` discards most markers at low coverage.**
  With bulks at ~30× mean depth, ~99% of markers fail it (section 14). Users
  sequencing bulks below ~60× should lower it (e.g. `-b 10,10`); the
  threshold simulation uses the observed depths, so lower-depth markers are
  handled correctly.
- **Support intervals were calibrated for F2 with the default window only.**

- Chart tooltips mark only upper-tail exceedances (★ p99 / ● p95); a value
  below a lower (valley) threshold, such as LowSI at a QTL, is not flagged.

- **Bulks-only runs (`-P None,None`)** can't tell which allele came from which
  parent, so ΔSI's sign is effectively random per marker and smoothing cancels
  part of the true signal. The null simulation reproduces this, so thresholds
  stay valid, but power is lost. G and LOD are not affected.
- Write errors on output files are not checked everywhere.
- Many files have Windows (CRLF) line endings, which shows up as whole-file
  changes in git. A `.gitattributes` with `* text=auto eol=lf` would fix it.
