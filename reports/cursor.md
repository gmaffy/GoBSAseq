# GoBSAseq Comprehensive Code Review

**Date:** July 10, 2026  
**Reviewer role:** Experienced computational biologist / Go pipeline developer  
**Repository:** `e:\GitHub\GoBSAseq`  
**Scope:** Full pipeline review — mode selection, filtering, statistics, smoothing, thresholds, QTL detection, plotting, orchestration  
**Verification:** Independent code read of all core packages; `go test ./...` passed (filter package only has tests)

---

## Executive Summary

GoBSAseq is a **serious, well-engineered BSA-seq pipeline** that goes well beyond a simple SNP-index calculator. It implements a complete workflow from multisample VCF input through hard filtering, multi-metric statistics, Gaussian smoothing, Monte Carlo null modelling, QTL interval calling, and interactive HTML plots. The statistical backbone — especially the **two-stage Monte Carlo null model** for deep sequencing and **population-structure-aware variance scaling** — reflects genuine expertise in bulk segregant analysis.

The pipeline is **suitable for exploratory and many production BSA-seq analyses** on diploid F₂/backcross populations with well-structured VCF input. It is **not yet publication-grade without caveats**: several README claims diverge from the implementation, peak calling uses simplified heuristics that can miss or merge QTLs, test coverage is thin, and scalability limits remain for very large genomes.

### Overall Score: **8.6 / 10**

| Dimension | Score | Notes |
|-----------|------:|-------|
| Implementation correctness (Go) | 9.0 | Solid concurrency, edge-case handling, numerically stable stats |
| Biological / statistical accuracy | 8.3 | Strong core methods; simplifying assumptions in polarity, RIL/F₃ nulls, LD |
| Pipeline completeness | 8.0 | VCF path complete; reads path is a stub; doc/output mismatches |
| Usability & reproducibility | 8.0 | Good CLI; interactive prompts; no checkpointing |
| Test coverage | 7.0 | Two integration tests in `filter/` only |
| Performance / scalability | 8.2 | Fast for typical datasets; in-memory design limits huge genomes |

**Prior reports in this repo** (`copilot.md`, `geminipro.md`, `maicode_report.md`, `vibe_report.md`, `gemini35flash.md`) largely agree on strengths (two-stage MC, multi-metric stats, pooled-aware filtering) and weaknesses (memory, FDR, test coverage). This review adds **code-verified findings** on documentation drift, unused parameters, and QTL-calling logic that prior reports under-emphasised.

---

## Pipeline Architecture (As Implemented)

```
User input (VCF + sample roles)
        │
        ▼
┌───────────────────┐
│ Mode detection    │  run/bsaseqType() — 10 modes from parent/bulk presence
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 1: Filter    │  filter/HardFilterVcf — GATK hard filter + BSA segregation filter
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 2: Raw stats │  stats/RawStats — SI, ΔSI, G, ED⁴, LOD, Bayes factor
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 3: Smooth    │  stats/SmoothAndNormalise — depth-weighted Gaussian + MAD Z-scores
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 4: Thresholds│  stats/CalculateThresholds — two-stage MC per depth + CompositeZ null
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 5: BRM       │  stats/RunBRM — analytical block detection on smoothed ΔSI / AFDev
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 6–8: QTLs    │  stats/detect — per-stat MC peaks, CompositeZ peaks, merge w/ BRM
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 9: Plots     │  plots/GeneratePlots — interactive HTML (go-echarts)
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ STEP 10: Genes    │  stats/GeneSpaceFromMerged (optional)
└───────────────────┘
```

The orchestration in `run/bsaseq()` is clear, sequential, and well-logged. Error propagation is generally good.

---

## 1. Analysis Mode Selection (Two-Bulk vs One-Bulk)

**Implementation:** `run/bsaseqType()` infers mode from which of `-P` (parents) and `-B` (bulks) are provided (use `None` to skip a role).

| Mode | Parents | Bulks | Statistics computed |
|------|---------|-------|---------------------|
| `2p2b` | both | both | ΔSI, G-stat, ED⁴, LOD, BB log-BF |
| `2phb`, `2plb`, `hp2b`, `lp2b`, `hphb`, … | partial | partial | One-bulk: AFDev, G, LOD, BB log-BF |
| `2b` | none | both | Two-bulk stats (no parental anchoring) |

**Strengths**
- Ten modes cover realistic experimental designs (Takagi-style two-bulk, single-bulk contrasts, bulks-only).
- Mode drives both filtering logic (`filter/bsaSeqFilterAllele`) and which statistics are written.

**Weaknesses**
- **No explicit `--mode` flag** — mode is inferred silently. Mis-labelled samples produce wrong statistics without a hard error.
- **`2b` (bulks-only):** `determineHighAllele()` defaults to REF (allele index 0) when no parent is present (`stats/stats.go:179`). ΔSI then measures REF-vs-ALT frequency shift with **no biological anchor** for which allele is trait-associated. Bulks-only BSA without parents is only valid if the user accepts arbitrary polarity or both bulks were crossed from known homozygous parents whose alleles happen to align with REF/ALT ordering — a fragile assumption.
- **Interactive sample selection** (`run/Run()`) blocks unattended/batch runs unless all four sample names are supplied on the CLI upfront.

**Recommendation:** Add `--mode` override with validation; for `2b`, require explicit `--high-allele REF|ALT` or warn prominently.

---

## 2. Variant Filtering (`filter/filter.go`)

**Score: 9.0 / 10**

### Biological accuracy

This is one of the strongest modules.

1. **Hard filtering** implements GATK best-practice SNP/INDEL thresholds with a **`LightFilter` mode** that correctly skips FS, SOR, MQRankSum, and ReadPosRankSum — filters that penalise the allele-frequency skew **expected** in pooled BSA samples.

2. **BSA segregation filter** enforces:
   - Minimum depth per role (`--parents-depth`, `--bulks-depth`)
   - Parent allele purity ≥ 85% (`parentAllelePurity = 0.85`)
   - ≤ 20% reads from non-target alleles in bulks (`maxOtherAlleleFrac = 0.20`)
   - Parents must carry **different** predominant alleles at informative loci

3. **Multi-allelic handling:** `BsaSeqTargetAlt()` tests each real ALT independently rather than discarding multi-allelic sites — important for complex loci.

4. **MNP classification fix:** Same-length substitutions with `len(REF) > 1` are treated as SNPs for filtering purposes, not silently dropped.

### Engineering

- Producer–consumer parallelism with **sequence-ID ordering** preserved in output VCF and tabix index — verified by `filter/filter_test.go`.
- Safe AD/RO/AO parsing avoids `vcfgo` panics on malformed FORMAT fields.
- Graceful recovery from `vcfgo` numeric parse errors on `.` missing values.

### Gaps

- No optional **MAF / polymorphism** filter between parents (sequencing errors can pass if depth thresholds are met).
- **Diploid-only** — no polyploid support (wheat 6×, etc.).
- No warning when required INFO fields (`QD`, `MQ`, …) are absent from the VCF header (variants may pass/fail unpredictably).

---

## 3. Raw Statistics (`stats/stats.go`)

**Score: 8.8 / 10**

### Methods implemented (verified)

| Statistic | Formula / approach | Assessment |
|-----------|-------------------|------------|
| **SNP Index (SI)** | `H / (H + L)` where H/L are counts of high- vs low-phenotype-associated alleles | Correct; uses AD not GT |
| **ΔSI** | `SI_high_bulk − SI_low_bulk` | Standard BSA metric (Takagi et al. 2013) |
| **G-statistic** | `2 Σ O ln(O/E)` on 2×2 table | Correct; **Yates +0.5** on two-bulk (`GStatistic`) |
| **ED⁴** | `(ΔSI)⁴` | Standard (Hill et al. 2013) |
| **LOD** | `log₁₀(L_alt / L_null)` with binomial likelihoods | Correct; epsilon clamping |
| **Bayes factor** | Beta-binomial with Beta(0.5, 0.5) prior | Correct; log-gamma stable |
| **One-bulk stats** | AF deviation from `ExpectedAF(population)` | Reasonable for single-contrast designs |

### Parent polarity (`determineHighAllele`)

Uses **allele depth** from the high (or low) parent to assign which VCF allele corresponds to the high phenotype — a sound choice when parents carry residual heterozygosity. Falls back to GT when depth is unavailable.

**Limitation:** Hard threshold (`refDepth >= altDepth` → REF is high allele). No confidence weighting; near 50/50 parents can flip polarity site-by-site.

### Issues found

1. **One-bulk G-statistic lacks Yates correction** (`oneBulkGStatistic`) while two-bulk uses +0.5 — minor inconsistency; can produce slightly anti-conservative one-bulk p-values at low counts.

2. **`ExpectedAF()` oversimplifies advanced populations:**
   - `F2`, `F3`, `F4`, `RIL` all return **p₀ = 0.5**
   - For **RIL**, most markers are fixed (0 or 1), not 0.5 — null expectations and one-bulk tests are miscalibrated at non-segregating loci (though filtering removes many of these).
   - Backcross parsing (`BC1H`, `BC2L`, …) is implemented and correct: p₀ = 1 − 0.5^(n+1) for high backcross.

3. **G-statistic parameter naming** passes phenotype-oriented counts (`HighBulkH`, `HighBulkL`, …) into parameters named `highBulkAlt/highBulkRef`. Logic is internally consistent but confusing for maintainers.

4. **No locus-level exclusion** when both bulks have identical allele counts (ΔSI = 0 everywhere) — variants are retained with zero signal.

---

## 4. Smoothing & Normalisation (`stats/smoothing.go`)

**Score: 9.2 / 10**

### Strengths

- **Depth-weighted Gaussian kernel:** `w = depth × exp(−0.5 × (d/σ)²)` — high-confidence variants dominate; biologically appropriate.
- **3σ cutoff + binary search** — O(n × m) per chromosome instead of O(n²).
- **Minimum 5 kernel contributors** (`minVariantsInKernel`) — suppresses spurious signal in marker-sparse regions (centromeres, low-recomb blocks).
- **Robust Z-scores:** MAD with 1.4826 consistency factor; mean/SD fallback when MAD = 0.
- **CompositeZ (Stouffer):** Uses minimal non-redundant set — `|ZΔSI| + ZGstat + ZLOD` (k = 3) — with explicit rationale for excluding correlated stats (ED⁴, BBLogBF). Runtime validation (`validateCompositeZStatCount`) prevents silent k mismatch with threshold simulation.

### Issues

1. **`StepSize` (`-s` / `--step-size`) is stored in config but never used.** Smoothing evaluates **every variant position** as a kernel centre. README and CLI describe a sliding window step — this is **documentation/implementation drift**. Either implement stepped evaluation or remove the flag.

2. **Fixed σ (= WindowSize / 2).** No adaptation to local marker density; telomeres (dense) and centromeres (sparse) get the same physical bandwidth.

3. **Global MAD normalisation** across all chromosomes — reasonable for autosomal BSA but can dilute signal on small chromosomes/contigs with few variants.

4. **No FDR / genome-wide multiple testing correction** on Z-scores.

---

## 5. Threshold Calculation (`stats/thresholds.go`)

**Score: 9.5 / 10 — best module in the pipeline**

### Two-stage Monte Carlo null model (verified correct)

```
Stage 1:  Alt alleles in pool ~ Binomial(n = bulk_size, p = p₀)
Stage 2:  Observed reads     ~ Binomial(n = depth,     p = realised_AF)
```

This correctly separates **finite-pool sampling variance** from **sequencing sampling variance**. At high coverage (>50×), single-stage models that draw reads directly from p₀ **underestimate null variance** and inflate false positives — a common failure mode in deep BSA-seq. GoBSAseq addresses this properly when `--bulk-sizes` (`-S`) is set.

### Additional rigour

- Per-depth threshold caching with `singleflight.Group` (no duplicate simulations).
- Parallel warm-up across all observed depths.
- **Full-pipeline CompositeZ null:** simulates raw stats → column-wise robust-Z → Stouffer combine, capturing correlation between ΔSI, G, and LOD under H₀.
- Empirical CompositeZ P99/P95 from simulated null (not fixed at 3.0).

### Issues

1. **Cache keyed on exact depth** — high depth variance (1×–500×) causes cache bloat. Depth binning (e.g. round to nearest 5 for depth > 50) would help.

2. **Default bulk size fallback = 100** when `-S` not provided (`calculateEmpiricalCompositeZThresholds`) — may miscalibrate thresholds if actual pools differ greatly.

3. **`--rep 1000`** may be insufficient for stable P99.5 tails; no convergence checking.

4. **Sequencing error rate not modelled** in null simulations.

---

## 6. QTL Detection (`stats/detect.go`)

**Score: 7.8 / 10 — functional but coarse**

### What works

- Per-chromosome peak finding with linear interpolation at threshold crossings.
- Separate tracks for each statistic (G-stat, ΔSI±, ED⁴, LOD, BBLogBF) using **MC-derived P99/P95** thresholds.
- CompositeZ peaks use **empirical** `CompositeZP99/P95` from full null simulation — not the hardcoded 3.0 from README.
- `ConsolidateQTLs()` merges overlapping peaks from different statistics into enriched intervals.

### Critical limitations

1. **`FindPeakIntersections()` uses a chromosome-wide average threshold** (`avgThresh`), not the per-variant depth-specific thresholds already computed in `thresholds[]`. A variant at 200× depth is compared against the same average as one at 20× — **wastes the per-depth MC work** and can mis-call peaks in mixed-depth datasets.

2. **`MergeCompositeBRM()` retains only the single best CompositeZ peak per chromosome** (`bestCompositePeakByChrom`) and one BRM summary per chromosome. **Multiple QTLs on the same chromosome are collapsed** into one merged interval — a major biological limitation for oligogenic traits.

3. **No minimum QTL width**, **no merge-distance parameter** (both commented out in `config.go`), **no FDR control**.

4. **README claims outputs `*.qtls.tsv` and `*.mc_qtls.tsv`** — **these files are not written.** Actual outputs are Excel workbooks:
   - `GoBSAseq.{mode}.individual.qtl.xlsx`
   - `GoBSAseq.{mode}.compositez.qtl.xlsx`
   - `GoBSAseq.{mode}.final.qtl.xlsx`

5. **README states primary threshold is |CompositeZ| ≥ 3.0** — code uses MC-derived empirical thresholds; 3.0 appears only as a **plot reference line** in `plots/plots.go` (`zSig = 3.0`).

---

## 7. BRM (`stats/brm.go`)

**Score: 8.5 / 10**

Despite the name, this is an **analytical threshold block caller**, not a full Bayesian regression model.

- Two-bulk threshold: `u_α × √[(n₁+n₂)/(V_scale·n₁·n₂) × p(1−p)]`
- One-bulk threshold: `u_α × √[p₀(1−p₀)/(V_scale·n)]`
- `PopulationVarianceScale()` correctly adjusts for F₂ (2.0), F₃ (4/3), F₄ (8/7), RIL (1.0), and backcross generations.

**Issues:** Sequential per-variant processing; local p estimated from average of smoothed SI values; no posterior credible intervals; name oversells the method.

---

## 8. Plotting (`plots/plots.go`)

**Score: 8.3 / 10**

- Three HTML pages per run: individual stats, Z-score overlay, composite signal.
- MC thresholds, BRM blocks, and z = ±2/±3 reference lines overlaid.
- Professional styling via go-echarts.

**Issues:**
- **No downsampling** — millions of variants produce multi-hundred-MB HTML files that hang browsers.
- No Manhattan-style genome-wide view, no PNG/PDF export, no BED output for IGV.
- Plot z = 3 line may **contradict** MC-derived CompositeZ thresholds shown elsewhere — confusing for users reading README.

---

## 9. Orchestration & Usability (`run/run.go`, `cmd/root.go`)

**Score: 8.0 / 10**

**Strengths:** 10-step pipeline; BAM → variant calling integration via `genome-whisperer`; optional gene-space annotation; basic config file support (`--config`).

**Weaknesses:**
- **Reads-based path is a stub** — `runType == "reads"` prints a message and returns without analysis.
- **No checkpoint/resume** — long MC runs restart from step 1 on failure.
- **Interactive prompts** when sample names missing — problematic for HPC/CI.
- Misleading log message at step 9: `"MC-based QTL detection completed"` even though MC QTL TSV output does not exist.

---

## 10. Test Coverage

**Score: 7.0 / 10**

`filter/filter_test.go` contains two valuable tests:
1. Multi-allelic ALT selection flows into raw stats correctly.
2. Hard-filter output preserves coordinate sort under parallelism.

**Missing tests (high value):**
- G-statistic, LOD, Bayes factor against known analytical values
- Two-stage MC null distribution properties (mean/variance at high depth)
- Smoothing kernel weights and CompositeZ consolidation
- Peak detection on synthetic QTL spike
- End-to-end mini-VCF regression with expected QTL interval

`go test ./...` passes; all other packages have **no test files**.

---

## Biological Accuracy Summary

### What GoBSAseq gets right (publication-quality concepts)

1. Multi-metric BSA (ΔSI, G, ED⁴, LOD, Bayes factor) — captures diverse genetic architectures.
2. Two-stage MC null for deep sequencing — **state-of-the-art** for pooled BSA.
3. Population-structure variance scaling in BRM and expected AF for backcrosses.
4. Pooled-sample-aware filtering (LightFilter).
5. Depth-based parent polarity instead of naive GT.
6. Stouffer consolidation with deliberate exclusion of redundant correlated statistics.

### Simplifications that limit biological fidelity

| Assumption | Impact |
|-----------|--------|
| Hardy-Weinberg / binomial pools | Moderate — ignores overdispersion from pooling bias |
| RIL/F₃ null p₀ = 0.5 | Moderate — wrong for fixed loci in near-homozygous lines |
| No LD modelling | Moderate — adjacent variants correlated; peaks appear broader |
| No sequencing error in null | Low–moderate |
| Bulks-only arbitrary polarity | High in `2b` mode |
| One QTL per chromosome in final merge | High for oligogenic traits |
| Uniform recombination / fixed Gaussian σ | Moderate — resolution varies by region |

---

## Documentation vs Implementation (Important)

Several README claims do **not** match the current code:

| README claim | Actual behaviour |
|-------------|-----------------|
| Primary QTL output: `*.qtls.tsv` with \|CompositeZ\| ≥ 3.0 | Excel files; threshold = MC empirical P99/P95 |
| Secondary: `*.mc_qtls.tsv` | Not implemented as separate output |
| `--step-size` controls sliding window step | **Parameter unused** in smoothing |
| MC QTL detection run by default | Individual-stat MC peaks yes; dedicated MC QTL TSV no |

These should be reconciled before citing the tool in a methods section.

---

## Roadmap to 10 / 10

### Priority 1 — Correctness & trust (≈ +0.5 points)

| # | Item | Effort | Impact |
|---|------|--------|--------|
| 1 | Use **per-variant thresholds** in `FindPeakIntersections`, not chromosome average | Low | High |
| 2 | Support **multiple QTLs per chromosome** in `MergeCompositeBRM` | Medium | High |
| 3 | Fix **README/output drift** (TSV exports, threshold description, StepSize) | Low | Medium |
| 4 | Implement or remove **`StepSize`** | Low | Medium |
| 5 | Add **`--high-allele`** for bulks-only mode | Low | Medium |

### Priority 2 — Statistical rigour (≈ +0.3 points)

| # | Item | Effort | Impact |
|---|------|--------|--------|
| 6 | **FDR control** (Benjamini–Hochberg) on genome-wide tests | Medium | High |
| 7 | **Probabilistic parent polarity** with confidence scores | Medium | Medium |
| 8 | Population-specific p₀ for **RIL** (skip or adjust fixed loci) | Medium | Medium |
| 9 | Yates correction on **one-bulk G-stat** | Trivial | Low |
| 10 | Model **sequencing error rate** in MC null (optional `--error-rate`) | Low | Medium |

### Priority 3 — Scale & usability (≈ +0.2 points)

| # | Item | Effort | Impact |
|---|------|--------|--------|
| 11 | **Chromosome-streaming** to cap memory on 50M+ variant genomes | High | High |
| 12 | **Plot downsampling** (retain all peaks, thin background) | Low | High |
| 13 | **Checkpoint/resume** after filter, smooth, or threshold steps | Medium | Medium |
| 14 | Complete **reads-based** analysis path | High | Medium |
| 15 | **Depth binning** for MC cache | Low | Low |

### Priority 4 — Validation (≈ +0.2 points)

| # | Item | Effort | Impact |
|---|------|--------|--------|
| 16 | Unit tests for all statistics against hand-calculated examples | Medium | High |
| 17 | Synthetic VCF integration test with planted QTL at known position | Medium | High |
| 18 | Benchmark against QTLseqR / BSA on published dataset | Medium | High |
| 19 | Bundled **toy dataset** with documented expected outputs | Low | Medium |

**Estimated effort to reach ~9.5/10:** 40–60 developer hours focused on items 1–6, 11–12, 16–17.

---

## Comparison to Existing BSA Tools

| Feature | GoBSAseq | Typical SNP-index scripts | QTLseqR (R) |
|---------|----------|--------------------------|-------------|
| Performance (Go) | Excellent | Moderate (Python/R) | Moderate |
| Two-stage MC null | Yes | Rarely | Partial |
| Multi-metric + Stouffer | Yes | Usually ΔSI only | Yes |
| Pooled-aware filtering | Yes | Variable | Yes |
| Interactive HTML plots | Yes | Rare | Yes |
| Maturity / validation | Early | Variable | High (published) |
| Polyploid support | No | Variable | Limited |

GoBSAseq's statistical engine is **competitive with or ahead of** most open-source BSA tools on null modelling; it lags on validation literature, edge-case handling, and output standardisation.

---

## Final Verdict

### Is GoBSAseq production-ready?

| Use case | Verdict |
|----------|---------|
| F₂ / backcross BSA-seq, standard genomes, well-annotated VCF | **Yes** — recommended with manual plot review |
| Deep sequencing (>50×) with accurate `-S` bulk sizes | **Yes** — two-stage MC is a major advantage |
| Oligogenic traits (multiple QTLs per chromosome) | **Caution** — merge logic will under-call |
| Bulks-only (`2b`) without parents | **Not recommended** without polarity workaround |
| Wheat-scale / 50M+ variants | **Caution** — memory limits |
| Publication methods text | **Reconcile README vs code first** |
| Polyploid crops | **Not supported** |

### Score justification

GoBSAseq earns **8.6/10** because it combines a **genuinely advanced statistical core** (two-stage MC, population variance scaling, thoughtful CompositeZ design) with **production-quality Go engineering** (concurrency, VCF robustness, modular pipeline). It loses points for **QTL-calling heuristics** (average thresholds, one-QTL-per-chromosome merge), **documentation drift**, **dead CLI parameters**, **thin tests**, and **biological simplifications** (RIL null, bulks-only polarity) that matter in real breeding programmes.

Implementing per-variant peak thresholds, multi-QTL chromosome support, FDR control, and chromosome streaming would raise this to **~9.3–9.5/10** — arguably best-in-class among open-source BSA-seq pipelines.

---

## Key Code References

Parent polarity via depth:

```133:179:stats/stats.go
func determineHighAllele(v *vcfgo.Variant, highParentIdx, lowParentIdx, altIdx int) int {
	// ... uses AD/RO predominant allele, falls back to GT
}
```

Two-bulk G-statistic with Yates correction:

```445:477:stats/stats.go
func GStatistic(highBulkAlt, highBulkRef, lowBulkAlt, lowBulkRef int) float64 {
	hba := float64(highBulkAlt) + 0.5
	// ...
}
```

Two-stage MC simulation:

```132:178:stats/thresholds.go
for i := 0; i < p.rep; i++ {
	// Stage 1: Binomial(bulk_size, p0) → realised AF
	// Stage 2: Binomial(depth, realised_af) → observed reads
}
```

Peak detection using average threshold (not per-variant):

```216:220:stats/detect.go
var threshSum float64
for _, t := range thresholds {
	threshSum += threshFn(t)
}
avgThresh := threshSum / float64(len(thresholds))
```

One QTL per chromosome merge:

```796:804:stats/detect.go
func bestCompositePeakByChrom(peaks []PeakIntersection) map[string]PeakIntersection {
	// keeps only highest CompositeZ peak per chromosome
}
```

---

**Reviewed by:** Claude (Computational Biology Agent)  
**Recommendation:** Suitable for deployment on standard BSA-seq workflows; address Priority 1 items before methods-section publication.
