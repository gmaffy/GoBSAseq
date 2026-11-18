# GoBSAseq — Path to 10/10

Code review of accuracy and biological soundness across the pipeline. Items are
grouped by module and ordered within each group by impact (🔴 correctness / can
change QTL calls, 🟠 calibration / statistical soundness, 🟡 robustness & polish).

---

## 1. Filtering (`filter/filter.go`, `filter/split.go`) — **8.5 / 10**

Strong, biologically sound, and well-tested (12 filter tests). Caller-aware,
AF-safe (the GQ-on-bulks caveat is now fixed), decomposition-correct, and
depth-bounded. Not yet a 10: the coverage cap is manual not adaptive (1.b), there
is no indel-proximity masking (1.c), per-allele INFO is dropped on split (1.f),
and the parse path is still single-threaded (speed, §7). Adding an adaptive
max-depth cap (1.b) is the next-biggest quality lever.

**Resolved this pass:**
- [x] 🔴 ~~**Strict biallelic gate** rejected sites over a single stray read.~~
  Replaced with a fractional second-ALT tolerance (`secondAltMaxFraction = 5%`)
  on *contaminating* alleles only — never truncates the REF-vs-target-ALT
  distribution, so no SNP-index bias.
- [x] 🟠 ~~**GATK filters silently no-op'd on DeepVariant VCFs.**~~
  `DetectFilterProfile` inspects the header; `PassesHardFilter` applies each
  threshold only when the annotation is present; `reportActiveFilters` prints the
  caller, mode, and active filters.
- [x] 🟠 ~~**Light mode wrongly dropped FS/SOR.**~~ It now keeps FS/SOR/MQ (strand
  bias stays valid in a pool) and drops only MQRankSum/ReadPosRankSum (unreliable
  near fixation, i.e. at QTLs).
- [x] 🔴 ~~**No maximum-depth cap** (collapsed-repeat/CNV artifacts passed).~~
  Added `--bulks-max-depth` / `--parents-max-depth` (0 = off).
- [x] 🟠 ~~**No genotype-quality gate for DeepVariant.**~~ Added `--min-GQ`
  (applied only when GQ present). *See new caveat 1.a below.*
- [x] 🟠 ~~**Mixed SNP+indel multiallelic records classified into one bucket.**~~
  Resolved by decomposition (`filter/split.go`, on by default;
  `--no-split-multiallelic` to disable) — a pure-Go `bcftools norm -m-` that
  reindexes GT and FORMAT R/A/G (AD, PL, GL) by header Number; each decomposed
  allele is then classified on its own.

**New findings from this review:**
- [x] 🟠 ~~**1.a — GQ floor applied to bulks (AF-correlated bias).**~~ Fixed:
  `bsaSeqFilterAllele` now gates bulks with `bulkOK` (depth + allele signal, no
  GQ) and only parents with the GQ floor (`parentSigOK` / `parentDepthOK`). The
  config/flag docs and the run report state "parents only"; a regression test
  (`TestMinGQAppliesToParentsNotBulks`) locks in that a low bulk GQ passes while a
  low parent GQ is rejected.
- [ ] 🟠 **1.b — Max-depth cap is a manual absolute per role.** The ideal is a
  data-driven cap (e.g. drop depth > k × median, or > 99th percentile per bulk),
  which requires a cheap pre-pass over depths. Add an `--max-depth-factor` (auto)
  alongside the absolute caps.
- [ ] 🟠 **1.c — No indel-proximity / low-complexity masking.** SNPs within a few
  bp of an indel carry alignment-artifact allele ratios. Add a position-window
  mask (needs a sorted sliding-window pass; awkward in the current out-of-order
  worker design — likely a small serial pre-pass or a post-collection filter).
- [ ] 🟡 **1.d — Validate `max-depth ≥ min-depth`.** If a user sets a max below the
  min, every variant is silently rejected. Warn or error at config time.
- [ ] 🟡 **1.e — Decomposition drops other-ALT reads from the SNP-index
  denominator.** Standard `norm -m-` behaviour, and genuine balanced triallelic
  segregation is rare in a biparental cross, but a real multiallelic site emits
  two records each over-estimating its allele frequency. Note it, and consider
  dropping sites with ≥2 substantial ALTs even when splitting. (When split is on,
  the 5% second-ALT tolerance is moot — it only guards `--no-split-multiallelic`
  runs; document that the two mechanisms are complementary, not redundant.)
- [ ] 🟠 **1.f — Per-allele INFO (AC, AF, MLEAC…) is dropped on split**, not
  reindexed, because vcfgo's `Info_` is unexported. Fine for this pipeline (only
  Number=1 INFO is read) but a cleaner decomposed VCF for external tools needs a
  serialize-and-reparse path. Documented in `split.go`.

**Remaining from the original review:**
- [ ] 🟠 **Make `parentAllelePurity = 0.85` configurable** for heterozygous/noisy
  parents.
- [ ] 🟡 **`BsaSeqTargetAlt` takes the first passing ALT**, not the best-supported.
  Now largely moot with decomposition on (records are biallelic), but still
  applies under `--no-split-multiallelic`.
- [ ] 🟡 **Warn when neither DP nor AD/RO is present** so `effectiveDP`'s zero
  isn't silently treated as "below threshold".

## 2. BSA-seq statistics (`stats/stats.go`, `stats/brm.go`)

- [ ] 🟠 **G-statistic continuity correction is applied to every cell always**
  (`+0.5` on all four counts). This is defensible (Anscombe-style) but biases G
  downward at low counts and differs from the usual "add 0.5 only to zero cells."
  Document the choice explicitly, or switch to Woolf/Haldane on zero cells only,
  and cite it.
- [ ] 🟠 **`GStatistic` parameter names don't match call sites.** Signature says
  `(highBulkAlt, highBulkRef, lowBulkAlt, lowBulkRef)` but callers pass
  `(HighBulkH, HighBulkL, LowBulkH, LowBulkL)` (high-allele / low-allele counts).
  G is symmetric so results are correct, but rename params to avoid a future
  editing error.
- [ ] 🟡 **ED⁴ is `(ΔSI)⁴`, a monotone transform of |ΔSI|.** It carries no
  information independent of ΔSI (the code already excludes it from CompositeZ
  for this reason). Keep it only as a display statistic and label it as derived,
  or replace with the canonical two-allele Euclidean distance if independence is
  intended.
- [ ] 🟡 **Document `ExpectedAF` for backcrosses.** `p0 = 1 − 0.5^(n+1)` for BCnH
  needs a derivation comment / reference so reviewers can confirm the expected
  allele frequency per cross.
- [ ] 🟡 Per-variant `Depth` is set to `min(highTotal, lowTotal)` — conservative and
  reasonable, but note it feeds the threshold simulation depth; document it.

## 3. Smoothing (`stats/smoothing.go`)

- [ ] 🟠 **Smoothing–null mismatch (most important statistical item).** Real
  CompositeZ is computed on Gaussian-smoothed, spatially autocorrelated values,
  but the Monte Carlo null (`thresholds.go`) draws *single, unsmoothed* SNPs.
  Both are robust-Z scaled by their own MAD, which only partly compensates —
  smoothing changes the tail shape. Simulate a smoothed window in the null
  (draw a neighbourhood and apply the same kernel), or derive the effective
  variance reduction analytically. This directly affects false-positive rate.
- [ ] 🟠 **Global (genome-wide) MAD normalisation loses power on small genomes /
  large QTLs.** If a big linked region occupies much of the genome, it inflates
  the MAD and shrinks Z everywhere. Consider per-chromosome or trimmed-null
  normalisation, and warn when a large fraction of markers exceed the threshold.
- [ ] 🟡 **`minVariantsInKernel = 5` silently drops sparse positions.** Make it
  configurable and report how many/where were dropped so real sparse regions
  (low-recombination, low-coverage) aren't invisibly lost.
- [ ] 🟡 Confirm and document that averaging per-SNP G/LOD across the kernel (rather
  than recomputing on pooled counts) is the intended Magwene-style G′; add the
  citation.

## 4. Threshold calculation (`stats/thresholds.go`)

- [ ] 🔴 **Raise default `rep` above 1000 for tail estimation.** The P99/P99.5
  quantiles are estimated from the top 0.5–1% of draws — only ~5–10 points at
  rep=1000, so thresholds are noisy and run-to-run unstable. Default to ≥10,000
  for the tail quantiles (caching keeps it affordable), or interpolate the tail.
- [ ] 🟠 **CompositeZ threshold uses a single median `repDepth` for all variants**,
  while per-stat thresholds are per-variant depth. Low-depth regions therefore
  get an anticonservative CompositeZ threshold. Either make CompositeZ depth-
  aware or document that CompositeZ is only valid at ~median depth.
- [ ] 🟠 **Verify the BRM analytical variance formula.** In
  `calculateBRMBlocksTwoBulk`: `pool*(1-1/d) + af(1-af)/d` per bulk. The
  `(1-1/d)` attenuation and the twofold binomial term deserve a written
  derivation (law of total variance) with a reference; unit-test it against the
  Monte Carlo two-bulk null at a few depths to confirm they agree.
- [ ] 🟡 **Seeding is not reproducible.** `nextSeed()` is a process-global counter,
  so results aren't reproducible across runs and can't be pinned in tests. Add a
  `--seed` flag and thread it through.
- [ ] 🟡 The `empiricalZThresholds` fallback (2.576/1.960) silently replaces the MC
  thresholds on error — log loudly and record in output which path was used.

## 5. QTL detection (`stats/detect.go`, `run/run.go`)

- [ ] 🔴 **`MinQTLWidth` and `MinQTLMarkers` are parsed but never used.** No
  minimum-width or minimum-supporting-marker filter is applied, so single-marker
  and ultra-narrow spurious peaks are reported as QTLs. Wire both into
  `ConsolidateQTLs` / final QTL emission.
- [ ] 🟠 **`findPeaksWithFallback` is all-or-nothing per chromosome.** If *any* p99
  peak exists on a chromosome it returns only p99 peaks and never reports p95
  (suggestive) regions elsewhere on that chromosome. Collect p99 and p95
  independently and annotate each peak with its level, rather than suppressing
  p95 whenever a p99 exists.
- [ ] 🟡 **`StepSize` is dead.** Smoothing is per-SNP Gaussian, not stepped windows;
  remove the flag or document it as vestigial to avoid implying a windowing
  behaviour that doesn't exist.
- [ ] 🟡 Remove dead helpers `minI64`, `maxI64` (`detect.go`) and the large
  commented-out `FindPeakIntersections` block.

## 6. Plotting (`plots/plots.go`)

- [ ] 🔴 **X-axis is a category axis of position *labels*, not a value axis.** SNPs
  are plotted at equal index spacing regardless of genomic distance, so peaks
  render at visually wrong genomic positions and gaps/dense regions are
  distorted. Switch to a numeric value axis so bp spacing is faithful.
- [ ] 🟠 **Individual charts draw flat *average* threshold lines** while detection
  uses per-variant thresholds. A marker can sit below the drawn line yet be
  called significant (or vice versa). Draw the per-variant threshold as its own
  series so the plot matches the calls.
- [ ] 🟡 Remove dead `absMax` in `plots.go`.
- [ ] 🟡 BRM mark-area index mapping relies on the category axis; re-verify after
  moving to a value axis.

## 7. Flow, architecture & correctness (`run/`, `cmd/`, repo)

- [ ] 🔴 **`runType == "reads"` returns `nil` immediately** — the reads pipeline is
  unimplemented but exits success silently. Either implement it or return a clear
  "not implemented" error.
- [ ] 🟠 **Fully non-interactive mode.** `Run` uses `fmt.Scan` for sample selection
  and the gene-space prompt, blocking scripting/CI. Allow every choice to be
  supplied by flag/config and only prompt when a TTY is attached.
- [ ] 🟠 **Harden the config-file parser (`cmd/root.go`).** It positionally splits
  reads on spaces and silently ignores malformed lines / unknown keys. Move to a
  typed format (YAML/TOML) with validation and clear errors.
- [ ] 🟠 **Add end-to-end and statistical tests.** Current tests cover peak
  interpolation, merging, and pool sampling. Add: a golden-file run on a tiny
  synthetic VCF, a null-calibration test (no QTL → false-positive rate ≈ α), and
  a positive-control test (planted QTL is recovered at the right locus).
- [ ] 🟡 **Decompose the 150-line `bsaseq()` step function** into named stage
  functions returning typed results; makes each of the 10 steps independently
  testable.
- [ ] 🟡 Delete `oldcode/` (backups belong in git history), and the per-model report
  files if not needed in-tree.
- [ ] 🟡 Replace ad-hoc `color.*`/`fmt.Print` with a leveled logger (`--verbose`,
  `--quiet`) and write a machine-readable run manifest (parameters, seed,
  versions, counts) alongside outputs for reproducibility.
- [ ] 🟡 Remove unused `OneParent*` config fields or wire them in.

---

## Suggested order of attack

1. **Correctness that changes calls:** unused QTL width/marker filters (§5),
   reads-runType no-op (§7), plot x-axis (§6), biallelic gate (§1).
2. **Statistical calibration:** smoothing↔null mismatch (§3), rep count (§4),
   depth-aware CompositeZ (§4), p95/p99 fallback (§5).
3. **Reproducibility & testing:** seeding (§4), null/positive-control tests (§7).
4. **Robustness & cleanup:** non-interactive mode, config parser, dead code,
   logging, decomposition.
