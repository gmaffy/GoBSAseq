# GoBSAseq review: biology, statistics, code

Scope: `stats/`, `filter/`, `run/`, `cmd/`, the three root reviews (`copilotfusion_review.md`, `gemini_review.md`, `grokreview.md`) and `indep_bsa.py`.
`go build ./...` and `go test ./...` pass. Every finding below was checked against the source, and the calibration claims were checked by simulation (scripts are listed in the appendix).

## Bottom line

The genetics is mostly right: pool genotype sampling for F2/F3/F4/RIL/BCn, `ExpectedAF`, `PopulationVarianceScale` and the law-of-total-variance BRM formula all check out by hand. Significance is the weak part. None of the three detectors (per-statistic, CompositeZ, BRM) controls the chance of a false QTL somewhere in the genome. In a simulated null F2 genome (no QTL, 20+20 plants, 50x, 10 chromosomes x 30 Mb, σ = 1 Mb):

| Detector, default settings | Null genomes with ≥ 1 false call |
|---|---|
| CompositeZ at "p99" (13.2) | **30%** (40% at 150x) |
| CompositeZ at "p95" fallback (8.0) | **80%** |
| BRM, α = 0.05 | **45%** |

`final.qtl.xlsx` is the union of CompositeZ and BRM, so its false-positive rate is higher than either one alone.

This contradicts `gemini_review.md` §2.1, which says smoothed CompositeZ "virtually never exceeds 4" and that the threshold causes catastrophic false negatives. Linked markers share the same pool genotypes, so the smoothed tracks are heavy-tailed blocks rather than near-Gaussian averages. In the simulation the pointwise null p99 of the smoothed CompositeZ was about 8.7 and the genome maximum about 10.3. The code's single-SNP threshold sits between those two values. That makes it too strict pointwise and too loose genome-wide, so it is not calibrated for either use.

---

## A. Findings no other review reported

### A1. `OneBulkSize` is never set, so one-bulk nulls use depth as the number of plants (HIGH)
`utils.AnalysisConfig.OneBulkSize` is read in `stats/thresholds.go:285,618` and `stats/brm.go:316`, but `cmd/root.go` never assigns it.
- `simulateOneBulk` falls back to `bulkSize = depth`. At 60x with a 20-plant bulk, the simulation assumes 60 plants. The pool-sampling variance is then about 3x too small, the AFDev/G1/LOD1 thresholds are too tight, and false positives follow in every `*hb` and `*lb` mode.
- CompositeZ and BRM fall back to `HighBulkSize`, which is the wrong bulk for `2plb`, `hplb` and `lplb`.

**Fix:** derive the size from the mode, and make a missing size a hard error.
```go
// utils/config.go
func (c AnalysisConfig) OneBulkN(bsaType string) (int, error) {
    n := c.HighBulkSize
    if strings.HasSuffix(bsaType, "lb") { n = c.LowBulkSize }
    if n <= 0 { return 0, fmt.Errorf("bulk size required for %s (-S)", bsaType) }
    return n, nil
}
```
Use this in `CalculateThresholds`, `calculateEmpiricalCompositeZThresholds` and `RunBRM`. Delete the "bulkSize = depth" and "= 100" fallbacks. Each one silently picks a different null.

### A2. Splitting multiallelic records hides a parent that carries the other ALT (MEDIUM)
`SplitMultiallelic` reduces AD to `[ref, alt_i]` and remaps every other ALT in GT to `0` (`filter/split.go:127`). Two problems follow:
1. A parent that is homozygous for ALT2 has AD `[0,0]` in the ALT1 record. `parentAllele` (`filter/filter.go:215–246`) sees a total of 0, falls back to the remapped GT `0/0`, and reports it as a **confident REF** parent. The site passes as informative, and the bulk SI is then measured as REF vs ALT1, where REF comes from neither parent.
2. The 5% second-ALT tolerance in `hasUsableAlleleSignal` stops working, because after the split `altDepths` has one element and `other` is always 0.

**Fix:** in `splitSample`, keep the off-target depth, for example `fields["_OTHER_AD"] = strconv.Itoa(sum(other alts))`. Read it in `hasUsableAlleleSignal` and `parentAllele`. In `parentAllele`, fall back to GT only when AD is absent, not when it sums to zero.
```go
if errR == nil && errA == nil { // AD present
    if total == 0 { return 0, false } // no informative reads: not confident
    ...
}
```

### A3. Bulks-only (`2b`) mode smooths a signed ΔSI whose sign is arbitrary at each SNP (HIGH for `2b`)
With no parent, `determineHighAllele` returns REF, so ΔSI = SI(REF, high) − SI(REF, low). Unless the reference genome *is* one of the parents, the REF allele comes from parent A at some SNPs and from parent B at others. The sign of ΔSI therefore flips from SNP to SNP. `gaussianSmooth` averages signed ΔSI, so at a real QTL the positive and negative values cancel toward 0. That suppresses `ZDeltaSI`, the `DeltaSI±` peaks and **BRM**, which tests the signed `SmDeltaSI`. G, LOD and ED⁴ do not depend on the sign and are unaffected.
`indep_bsa.py` does this correctly: it smooths `|dsi|` and simulates random polarity in its null.

**Fix:** in `2b` mode, smooth `|ΔSI|`, threshold it against a null that also flips polarity at random, and drop the `DeltaSI-` track. Alternatively, add a `--ref-is-high-parent` / `--ref-is-low-parent` flag for when the reference genome is a parent.

### A4. Single-parent and bulks-only modes have no segregation filter (MEDIUM, HIGH for single-bulk modes)
In `hp2b`, `lp2b`, `2b` and every single-bulk mode, nothing checks that a site actually segregates in the cross. Many sites are polymorphic only in other samples of a multi-sample VCF, or are shared differences from the reference.
- In two-bulk modes those sites add ΔSI ≈ 0 noise and dilute the signal.
- In `hphb`, `hplb`, `lphb` and `lplb`, a site where both parents share the allele gives AFDev = ±0.5 and a very large G1. That produces genome-wide false signal.

The lone parent is also not required to be homozygous. `parentSigOK` checks depth and second-ALT reads only, so a 55/45 parent still sets the polarity.

**Fix:** apply `parentAllele(...).confident` (the 85% purity rule) to the lone parent. For `2b`, `hp2b` and `lp2b`, require a pooled bulk allele frequency in [0.2, 0.8], as `indep_bsa.py` does. For single-bulk modes, require at least 2 reads of each allele in the bulk, and print a warning that one-parent, one-bulk designs cannot tell a fixed QTL apart from a non-segregating site.

### A5. In one-bulk mode, G1 and LOD1 are the same statistic, so CompositeZ counts it twice (MEDIUM)
`oneBulkGStatistic` = 2·ln(10)·`oneBulkLOD` exactly: both are the binomial likelihood ratio without pseudo-counts. Robust Z is affine-invariant, so `ZOneBulkG ≡ ZOneBulkLOD`, and one-bulk CompositeZ = (|Z_AFDev| + 2·Z_G)/√3.
In two-bulk mode the same holds up to the +0.5 pseudo-count in `GStatistic`.
The comment in `consolidate()` that describes "low-redundancy" statistics is therefore wrong. **Fix:** drop LOD from CompositeZ (k = 2), or drop CompositeZ altogether (see B1).

---

## B. Findings shared with the other reviews, with my verdict

### B1. Thresholds do not match the scanned statistic and are not genome-wide (HIGH; agree with grok #1, #2, #4)
- **Per-statistic thresholds** (`simulateTwoBulk`) come from single, unsmoothed SNPs but are applied to `Sm*` tracks smoothed with σ = 1 Mb. Smoothing removes read noise but not pool noise, because linked SNPs share the same plants. Calibration therefore changes with depth.
- **The CompositeZ null** computes `Z(|ΔSI|)`, while the real pipeline computes `|Z(ΔSI)|`. These are different random variables. The null also uses one median depth, no smoothing and no linkage.
- **Pointwise thresholds are used as QTL calls.** See the table at the top.

**Fix (recommended), following the `null_maxima` approach in `indep_bsa.py`:**
1. For each null replicate and each chromosome, simulate 2N gametes per bulk with crossovers (Haldane, user-set genetic length, default 1 M per chromosome, or a cM/Mb rate). Assign pool frequencies at the **real marker positions**.
2. Draw reads at each marker's **real per-bulk depths**.
3. Pass the replicate through the **same** `RawStats → gaussianSmooth → normalise → consolidate` code. Refactor so these functions take `[]BSAstats` in memory and write no files.
4. Record the genome-wide maximum of each scanned track. The 95% and 99% quantiles of those maxima become the genome-wide thresholds. Use ≥ 1000 replicates. The simulation can thin markers, but smoothing on a thinned set overstates read noise, so thinning should stay modest.

This one change replaces `simulateTwoBulk`, `simulateOneBulk`, `simulateNullCompositeZ*` and the per-depth caches, and fixes B1, B2, B4 and B6 together.

**Cheaper alternative:** Magwene's G′. Smooth G, fit a log-normal null to the Hampel-trimmed log G′, then apply Benjamini–Hochberg. This is what `indep_bsa.py` does, and QTLseqr uses the same method.

Also: compare smoothed statistics in their natural units against the null maxima. The robust Z is relative to the observed genome, so a large QTL inflates the MAD and changes every other Z. Keep the Z columns for plots only.

### B2. BCnF2 is simulated as BCn (HIGH for BC*F2; agree with all three)
`samplePoolAF` (`stats/thresholds.go:129`) strips `"F2"` and draws only recurrent-homozygous or heterozygous plants. BC1F2 should be 5/8 : 2/8 : 1/8, which gives twice the per-plant variance of BC1. `PopulationVarianceScale` is already correct (I checked BC1F2 = 1.5 and BC2F2 = 1.4 by hand).
```go
q := math.Pow(0.5, float64(n))          // heterozygote frequency in BCn
pHet, pDonor := q, 0.0
if isF2 { pHet, pDonor = q/2, q/4 }     // one generation of selfing
rc := 0                                 // recurrent-parent allele copies
for range bulkSize {
    u := rng.Float64()
    switch { case u < pDonor: rc += 0; case u < pDonor+pHet: rc += 1; default: rc += 2 }
}
f := float64(rc) / float64(2*bulkSize)
if !isHigh { f = 1 - f }                // SI is the high-parent allele
return f
```
Add a test asserting that the sampled variance equals p(1−p)/(scale·n) for each population, using `PopulationVarianceScale`.

### B3. Unequal bulk depths are collapsed to the minimum (MEDIUM; agree)
`BSAstats.Depth = min(high, low)` (`stats/stats.go:319`), and both simulated bulks use it (`thresholds.go:758,886`). **Fix:** store `HighDepth` and `LowDepth`. B1's fix uses real per-marker depths and solves this directly. If you keep per-SNP thresholds, cache by `(hd, ld)` pairs with depths binned (for example to 5x).

### B4. Tail labels (LOW–MEDIUM; agree with grok #3)
For two-sided ΔSI and AFDev, the "P95" pair uses 0.95 and 0.05, which is a 90% two-sided interval. Use 0.975 and 0.025. One-sided "P99" statistics use 0.995. Choose one convention and label it correctly. With 1000 replicates, the 0.995 quantile rests on about 5 draws.

### B5. One-bulk BRM treats `SmAFDev == 0` as missing (LOW; agree)
`stats/brm.go:250`. In low-bulk modes `SmHighSI` is 0, so the deviation becomes −p₀. **Fix:** `deviation := s.SmAFDev`.

### B6. BRM α is pointwise (MEDIUM)
`uAlpha = Φ⁻¹(1−α/2)` is applied to every SNP and compared with the single-SNP variance. In the null simulation this gave a 45% chance of at least one false block per genome. **Fix:** use the B1 null maxima for `|SmΔSI|`, or Bonferroni over the effective number of independent tests (about the genome's genetic length divided by the smoothing span).

### B7. Allele polarity via `parentAllele` counts all non-REF reads (LOW; agree with grok)
This is covered by A2. Count the target ALT only.

### B8. Reproducibility (MEDIUM; agree)
Seeds come from a global counter that goroutines take in scheduling order (`thresholds.go:91`), so thresholds change from run to run. **Fix:** add `--seed`, and make each simulation's seed depend only on `(seed, depth/chrom, replicate)`.

---

## C. Code-level issues (verified)

| # | Location | Problem | Fix |
|---|---|---|---|
| C1 | `detect.go:363` `findPeaksWithFallback` | Once a chromosome has any p99 peak, every p95 peak on it is dropped | Collect both levels, then merge. The p99 label wins on overlap. |
| C2 | `detect.go` `mergePeakIntersections` | Keeps the larger `PeakValue`, which is the wrong summit for lower-tail (`DeltaSI-`, `AFDev-`) tracks | Compare by tail: `|v|` or the sign of the tail |
| C3 | `detect.go` `addEvidence` | `DeltaSI+` and `DeltaSI-` (and AFDev ±) both write `q.DeltaSI`; the last one wins regardless of strength | Keep the peak with the larger \|value\|, or store separate +/− fields |
| C4 | `detect.go` `mergeConsolidatedQTL` | Keeps the first peak, not the strongest | Compare values |
| C5 | `RunBRM`, `DetectIndividualQTLs`, `DetectCompositeZQTLs` | Iterate over a Go map, so chromosome order in the output is random | Sort the keys |
| C6 | `smoothing.go:311` | Output keeps input order, and peak finding assumes POS is sorted within each chromosome. An unsorted VCF silently breaks intervals. | Sort `valid` by (CHROM, POS), or reject unsorted input |
| C7 | `cmd/root.go:420–422` | `--min-qtl-width`, `--min-qtl-markers` and `--step-size` are parsed and never used | Apply them in `ConsolidateQTLs` and `DetectCompositeZQTLs`, or remove the flags |
| C8 | `run/run.go:475` | The reads path prints a message and returns nil (exit 0) | Return an error |
| C9 | `run/run.go:52,446` | `fmt.Scan` prompts block batch jobs | Add `--non-interactive`, or check `term.IsTerminal(os.Stdin)` |
| C10 | `cmd/root.go` | `--min-FS-SNP` and `--min-SOR-SNP` are actually maximums | Rename to `--max-*`, keeping the old names as deprecated aliases |
| C11 | `smoothing.go` | Unplaced scaffolds, `Chr00` and organelle contigs are included in the genome-wide MAD and in calling | Add `--exclude-chroms` (regex) |
| C12 | `detect.go` | `minI64` and `maxI64` are dead code | Remove |

## D. Documentation
- Reference 2 in the README is the wrong Magwene paper. The G-statistic paper is Magwene, Willis & Kelly (2011), *PLoS Comput Biol* 7:e1002255.
- "LOD > 3 = significant" is a rule of thumb for human pedigree linkage and is not what the code uses. Remove it.
- Describe BBLogBF and LOD as read-level scores: they ignore pool size, so Kass–Raftery and χ² cut-offs do not apply.
- ED⁴ in the code is (ΔSI)⁴. The literature value for a biallelic site is 4·(ΔSI)⁴. The difference is cosmetic, but cut-offs from papers will not transfer.

## E. Notes on `indep_bsa.py`
This is a sound cross-check, and its null (linked F2 haplotypes, random polarity, smoothed |ΔSNP|, genome-wide maximum) is the model GoBSAseq should adopt. Minor points:
- **Null replicates:** 200 replicates leave the 99% quantile resting on about 2 maxima. Use ≥ 1000.
- **Marker thinning:** smoothing on a thinned marker set overstates read noise, so the threshold is somewhat conservative.
- **Genetic length:** the null assumes 1 M per chromosome. Real plant chromosomes are often longer. More crossovers mean more independent blocks and a higher maximum, so the threshold is slightly anti-conservative.
- **G′ null fit:** the code uses `np.std` of all retained log G′ values. Magwene and QTLseqr estimate spread from the left of the mode, so the QTL tail inflates the SD here, which is conservative.
- **Hard-coded inputs:** `NH, NL = 20, 16` and `Cmo_Chr00` should become arguments. The line `g = ... if False else None` is dead code.

## Suggested order of work
1. A1 (bulk size), B2 (BCnF2), B5. Small and safe.
2. A2, A4 (filter correctness), A3 (`2b` polarity).
3. B1 (chromosome-scale null run through the real pipeline), which also covers B3, B4, B6 and B8. Remove LOD from CompositeZ (A5).
4. C1–C12 and the documentation.
5. Tests: G and LOD on known 2×2 tables; AF and variance per population (including BC1F2); one-bulk size wiring; a `parentAllele` case with ALT2 after a split; a null-genome test asserting that the FWER at the 95% threshold is ≤ 0.05 ± MC error.

## Appendix: simulation
`czsim.py` and `brmsim.py` are in `/tmp/claude-1000/-mnt-d-GitHub-GoBSAseq/60710c43-b713-4657-8781-15acef220dd5/scratchpad/`.
Setup: F2, 20+20 plants, 50x (and 150x), SNPs every 5 kb, 10 chromosomes x 30 Mb, 1 M per chromosome, depth-weighted Gaussian with σ = 1 Mb and a 3σ cutoff. The robust Z, the CompositeZ formula and the code's single-SNP null replicate the Go functions.
