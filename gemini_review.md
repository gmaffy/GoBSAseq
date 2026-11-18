# GoBSAseq Codebase Review: Biological, Statistical, and Code Accuracy

This report provides a comprehensive evaluation of the **GoBSAseq** codebase across biological accuracy, statistical correctness, and software implementation.

---

## Executive Summary & Scorecard

| Category | Score | Primary Strengths | Critical Vulnerabilities |
| :--- | :---: | :--- | :--- |
| **Biological Accuracy** | **7.5 / 10** | Multiallelic decomposition (`bcftools norm -m-` equivalent); caller-aware hard filtering (GATK vs. DeepVariant); finite-pool diploid dosage modeling in $F_2/F_3/\text{RIL}$. | $BC_n F_2$ genotypic sampling completely omits homozygous recessive donor alleles; category-based X-axis distorts physical chromosomal distances and recombination slopes; arbitrary REF polarity in bulks-only mode (`2b`). |
| **Statistical Correctness** | **6.0 / 10** | Accurate law-of-total-variance derivation for finite pooled sampling; robust Z-score standardization using Rousseeuw–Croux MAD. | **Severe Smoothing–Null Distribution Mismatch** (comparing 2 Mb Gaussian-smoothed data against single-marker $\chi^2_1$ simulation); Stouffer's Composite Z combines collinear statistics; absence of genome-wide FWER/FDR multiple testing correction; depth asymmetry in two-bulk MC simulation. |
| **Code Implementation** | **7.0 / 10** | Concurrent worker pools (`GOMAXPROCS`, atomic work distribution); robust VCF error handling preventing `vcfgo` panics; pure-Go tabix indexing. | Interactive `fmt.Scan` breaks headless/CI pipelines; `--hp-reads` exits 0 silently without running; `MinQTLWidth`, `MinQTLMarkers`, and `StepSize` flags are parsed but dead; $p95$ peaks are suppressed on chromosomes with $p99$ peaks. |
| **Overall Score** | **6.8 / 10** | *High potential pipeline with strong infrastructure, hindered by statistical calibration mismatches and CLI control flow bugs.* |

---

## 1. Biological Accuracy and Soundness

### 1.1. Cross Structure & Genotypic Sampling in $BC_n F_2$ Progeny
* **File:** `stats/thresholds.go:125-139`
* **Status:** 🔴 **Critical Biological Error**

In `samplePoolAF`:
```go
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
        return float64(copies) / float64(2*bulkSize)
    }
}
```

#### Biological Breakdown:
1. In plant breeding and genetics, a $BC_n F_2$ population is created by **selfing** a backcross generation ($BC_n$). The entire biological rationale for creating $BC_n F_2$ lines is to allow recessive alleles from the donor parent to segregate as **homozygotes** ($LL$) so they can be phenotyped in the bulks.
2. Under Mendelian segregation:
   - For $BC_1$: $50\% HH$ and $50\% HL$.
   - Selfing $BC_1$ to $BC_1 F_2$: $HH \to 100\% HH$ ($0.5 \times 1.0 = 5/8 = 62.5\%$). $HL \to 25\% HH, 50\% HL, 25\% LL$ ($0.5 \times 0.25 = 1/8 = 12.5\% LL$; $0.5 \times 0.50 = 2/8 = 25\% HL$).
3. In `samplePoolAF`, line 129 executes `suffix = strings.TrimSuffix(suffix, "F2")` and models every individual as **either $HH$ or $HL$**. It assumes **$0\%$ $LL$ homozygotes**.
4. **Internal Code Contradiction:** In `stats/brm.go:86-91`, the analytical variance scale factor for $BC_1 F_2$ was correctly derived as $V_{scale} = 1.5$:
   $$\text{Var}(X) = 0.50 \implies \text{Var}(p_{pool}) = \frac{0.50}{4N} = \frac{3/16}{1.5 N} = \frac{p(1-p)}{1.5 N}$$
   However, the Monte Carlo simulation in `thresholds.go` completely omits the selfing step, producing a null simulation that contradicts its own analytical formula.

---

### 1.2. Allele Polarity in Bulks-Only Analysis (`2b` Mode)
* **File:** `stats/stats.go:133-180`
* **Status:** 🟠 **Biological Caveat**

```go
func determineHighAllele(v *vcfgo.Variant, highParentIdx, lowParentIdx, altIdx int) int {
    ...
    return 0 // Defaults to ref (0) when no parent is present.
}
```
* In `2b` mode (two bulks sequenced without parental lines, common in QTL-seq when parents are unavailable), `highAllele` unconditionally defaults to the **Reference (REF)** allele.
* If the causative or enriched QTL allele in the High bulk is the **Alternate (ALT)** allele, $\Delta SI = SI_{high} - SI_{low}$ becomes negative ($\Delta SI < 0$).
* While `CompositeZ` takes $|\Delta SI|$ and individual QTLs track `DeltaSI-`, the raw TSVs and plots present $\Delta SI$ inverted relative to standard phenotype-enrichment conventions.
* **Recommendation:** Document that in `2b` mode, $\Delta SI > 0$ denotes REF enrichment in High bulk, or allow the user to designate reference polarity.

---

### 1.3. Distortion of Chromosomal Coordinates via Category X-Axis
* **File:** `plots/plots.go:365-373`
* **Status:** 🔴 **Visualization & Biological Interpretation Defect**

```go
line.SetXAxis(positionLabels(x))
```
* Genomic positions are passed as string categories (`positionLabels(x)`).
* In Echarts, category axes assign an **equal horizontal spacing to every array element** regardless of genomic distance.
* **Biological Consequence:** Recombination deserts (centromeres) with low marker density shrink into narrow slivers, while dense SNP clusters (100 SNPs across 5 kb) expand to consume vast visual canvas. The user cannot visually evaluate physical QTL interval sizes, peak shoulders, or recombination slopes.

---

### 1.4. Indel Proximity & Low-Complexity Artifacts
* **File:** `filter/filter.go:373-428`
* **Status:** 🟠 **Filter Sensitivity Gap**
* Short-read aligners (BWA, Bowtie2) suffer from high mismatch and gap-opening ambiguity near insertions and deletions. SNPs located within 5–10 bp of an indel exhibit heavily biased allele depths ($AD$).
* GoBSAseq applies GATK hard filters and second-ALT tolerances, but lacks an indel-proximity window mask. Alignment artifacts near indels can manifest as false-positive single-marker spikes.

---

## 2. Statistical Accuracy and Correctness

### 2.1. The Smoothing–Null Model Mismatch (Primary Statistical Flaw)
* **Files:** `stats/smoothing.go:156-227`, `stats/thresholds.go:160-240`, `stats/thresholds.go:454-517`
* **Status:** 🔴 **Critical Statistical Inconsistency**

GoBSAseq smooths raw variant statistics across chromosomes using a depth-weighted Gaussian kernel ($\sigma = \text{WindowSize} / 2$, spanning $3\sigma = 3 \text{ Mb}$ by default). 

#### Mathematical Proof of Mismatch:
1. **Real Data Smoothing:** Let $X_i$ be the raw test statistic at marker $i$ (e.g. $G$-statistic). The smoothed statistic is:
   $$\bar{G}_k = \frac{\sum_m w_m G_m}{\sum_m w_m}$$
   Under $H_0$, single-marker $G_m \sim \chi^2_1$.
   When averaging $M$ independent markers in a window, by the Central Limit Theorem:
   $$\bar{G}_k \xrightarrow{d} \mathcal{N}\left(1, \frac{2}{M_{eff}}\right)$$
   Standardizing by its genome-wide MAD yields a robust Z-score $Z_{\bar{G}} \sim \mathcal{N}(0, 1)$, whose 99th percentile is:
   $$Z_{p99} \approx 2.33$$
2. **Null Simulation in `thresholds.go`:**
   The Monte Carlo simulation draws **isolated, unsmoothed single markers**.
   A single marker has $G \sim \chi^2_1$, with median $\approx 0.455$, $\text{MAD} \approx 0.45 \times 1.4826 \approx 0.67$, and 99th percentile $\approx 6.63$ (or higher with continuity correction).
   The simulated robust Z-score for a single marker reaches:
   $$Z_{\text{single}, p99} = \frac{6.63 - 0.455}{0.67} \approx 9.25 - 11.1$$
3. **The Resulting CompositeZ Threshold:**
   When summing $|Z_{\Delta SI}| + Z_G + Z_{LOD}$ for single markers, the empirical 99th percentile threshold generated by `CalculateThresholds` is **$q99 \approx 9.7 - 15.2$**.
   In real data, the smoothed CompositeZ is an average of Gaussianized tracks, which under $H_0$ virtually never exceeds **$4.0$**.
4. **Impact:** Real biological QTLs with smoothed $CompositeZ \in [4.0, 9.0]$ are **rejected as statistically insignificant**. The pipeline generates catastrophic false negatives.

---

### 2.2. Collinearity & Stouffer's Composite Z-Score Formulation
* **File:** `stats/smoothing.go:451-540`
* **Status:** 🔴 **Statistical Formulation Flaw**

```go
czs = append(czs, math.Abs(s.ZDeltaSI), s.ZGstat, s.ZLOD)
s.CompositeZ = math.Round((sum/math.Sqrt(float64(len(czs))))*1e6) / 1e6
```

#### Theoretical Breakdown:
1. **Mathematical Equivalence of $G$ and $LOD$:**
   In a $2 \times 2$ contingency table of read counts:
   $$\text{LOD} = \log_{10}\left(\frac{L_1}{L_0}\right) = \frac{\ln(L_1 / L_0)}{\ln 10}$$
   $$G = 2 \ln\left(\frac{L_1}{L_0}\right) = 2 \ln(10) \times \text{LOD} \approx 4.605 \times \text{LOD}$$
   $G$ and $\text{LOD}$ are the exact same likelihood ratio test. Their standardized Z-scores $Z_{Gstat}$ and $Z_{LOD}$ are collinear with pairwise correlation $\rho > 0.999$.
2. **Violation of Stouffer's Assumption:**
   Stouffer's formula assumes $k$ independent variables:
   $$\text{Var}\left(\frac{\sum_{i=1}^k Z_i}{\sqrt{k}}\right) = \frac{1}{k}\left(\sum_{i=1}^k \text{Var}(Z_i) + \sum_{i \neq j} \text{Cov}(Z_i, Z_j)\right)$$
   For $k=3$ nearly identical variables ($\rho \approx 1$):
   $$\text{Var}\left(\frac{3 Z_1}{\sqrt{3}}\right) = \frac{9}{3} \text{Var}(Z_1) = 3$$
   The true standard deviation is $\sqrt{3} \approx 1.732$, not $1.0$.
3. **Half-Normal Skew:** Taking $|Z_{\Delta SI}| \ge 0$ creates a half-normal distribution with mean $\sqrt{2/\pi} \approx 0.798$. The composite sum is neither zero-mean nor symmetric.
4. **Catastrophic Fallback:** If simulation fails, line 358 falls back to $Z_{crit} = 1.960$. Evaluating a variable with mean $\approx 0.8$ and standard deviation $\approx 1.73$ against 1.960 yields a false positive rate exceeding $25\%$.

---

### 2.3. Absence of Multiple Testing Correction (FWER / FDR)
* **File:** `stats/detect.go:363-382`
* **Status:** 🔴 **Multiple Testing Inflation**
* Detection evaluates variants against **pointwise** quantiles ($p95$, $p99$).
* In an organism with 200 independent resolution windows (genome size / window bandwidth):
  $$E[\text{False Positives at } \alpha = 0.05] = 200 \times 0.05 = 10 \text{ spurious QTLs per genome}$$
* Standard quantitative genomics packages (e.g. *QTLseqr*, *R/qtl*, *TASSEL*) apply Benjamini–Hochberg False Discovery Rate (FDR) or genome-wide permutation testing (Churchill & Doerge 1994) to control the family-wise error rate.

---

### 2.4. Global MAD Normalization Distortion by Major QTLs
* **File:** `stats/smoothing.go:371-381`
* **Status:** 🟠 **Power Loss in Asymmetric Architectures**
* Robust Z-scores are computed globally across all chromosomes combined:
  $$Z = \frac{X - \text{median}(X)}{1.4826 \times \text{MAD}(X)}$$
* In crosses with a dominant major gene (e.g., a resistance locus spanning 15 Mb in a 150 Mb Arabidopsis genome, or 10% of total markers), the massive linked region inflates the genome-wide MAD.
* This inflated denominator artificially deflates Z-scores across all other chromosomes, masking real minor/modest QTLs.

---

### 2.5. Sequencing Depth Asymmetry in Two-Bulk Null Model
* **File:** `stats/stats.go:319-321`, `stats/thresholds.go:886-887`
* **Status:** 🟠 **Sampling Variance Overestimation**
* In `stats.go`:
  ```go
  s.Depth = lowTotal
  if highTotal < s.Depth { s.Depth = highTotal }
  ```
* In `thresholds.go`:
  ```go
  p.highDepth = sm.Depth
  p.lowDepth = sm.Depth
  t.TwoBulk = simulateTwoBulkCached(p)
  ```
* `simulateTwoBulk` sets **both** bulk depths to $\min(\text{depth}_{\text{high}}, \text{depth}_{\text{low}})$.
* In unbalanced libraries (e.g., High bulk = 80x, Low bulk = 25x), the true variance of $\Delta SI$ is $\text{Var}_{pool} + \frac{p(1-p)}{80} + \frac{p(1-p)}{25}$. The simulation models $(25x, 25x)$, overestimating sequencing noise and artificially raising thresholds.
* For CompositeZ, a single median depth (`repDepth`) is used genome-wide, under-correcting low-depth regions and over-correcting high-depth regions.

---

### 2.6. Continuity Correction in the G-Statistic
* **File:** `stats/stats.go:445-450`
* **Status:** 🟡 **Methodological Discrepancy**
* `GStatistic` adds $+0.5$ to all four cells unconditionally (Anscombe-style pseudo-counts).
* In canonical biometry (Woolf 1957, Sokal & Rohlf 1995, Magwene et al. 2011), $+0.5$ is added **only to cells with count 0** (Haldane correction) to avoid $\ln(0)$. Adding $0.5$ unconditionally biases $G$ downward at low-to-moderate counts.

---

## 3. Code Accuracy, Architecture, and Implementation

### 3.1. Interactive Prompts Block Headless / Automated Execution
* **File:** `run/run.go:442-456`, `run/run.go:564-640`
* **Status:** 🔴 **Operational Blocker**

1. **Gene Space Prompt:** If gene-space flags are omitted, line 446 executes `fmt.Scan(&answer)`. In non-interactive environments (Nextflow, Snakemake, Slurm, Docker), this hangs indefinitely or immediately errors with `unexpected newline / EOF`.
2. **Sample Selection with `-P None`:**
   In `cmd/root.go`, default flag values are `","`. If a user runs `-P None,None -B hb,lb`:
   Line 581 in `run.go` checks:
   ```go
   if !slices.Contains(sampleNames, cfg.HighParentName)
   ```
   Because `"None"` is not present in the VCF header sample names, it triggers `promptForSampleChoice`, prompting on stdin even though the user explicitly supplied `-P None,None`.

---

### 3.2. Unimplemented Reads Pipeline Returns Exit Code 0 Silently
* **File:** `run/run.go:475-478`
* **Status:** 🔴 **Silent Pipeline Failure**

```go
if runType == "reads" {
    fmt.Printf("Running BSAseq with read-based analysis\n")
    return nil
}
```
* When passing raw FASTQ reads (`--hp-reads`, `--hb-reads`, etc.), the program logs a message and terminates with exit code 0 without executing BWA alignment or variant calling.
* **Fix:** Either call the alignment module or return a clear `fmt.Errorf("read-based analysis pipeline is not yet implemented")`.

---

### 3.3. Dead Flags and Unused QTL Region Parameters
* **Files:** `stats/detect.go:591-666`, `cmd/root.go:214-221`
* **Status:** 🔴 **Functional Incompletion**

* `--min-qtl-width` (`MinQTLWidth`): Parsed into config, but never referenced in `DetectIndividualQTLs`, `DetectCompositeZQTLs`, or `MergeCompositeBRM`. Single-base transient crossings are emitted as QTLs.
* `--min-qtl-markers` (`MinQTLMarkers`): Parsed, but never checked during region consolidation.
* `--step-size` (`StepSize`): Never used; smoothing is marker-centered, not stepped.
* Dead helpers: `minI64` and `maxI64` in `detect.go`; `absMax` in `plots.go`.

---

### 3.4. Peak Calling Suppresses Suggestive Peaks on Strong Chromosomes
* **File:** `stats/detect.go:363-382`
* **Status:** 🟠 **Detection Logic Flaw**

```go
func findPeaksWithFallback(...) []PeakIntersection {
    ...
    if peaks := tryLevel(threshP99, "p99"); peaks != nil {
        return peaks
    }
    return tryLevel(threshP95, "p95")
}
```
* If a chromosome possesses a major locus exceeding $p99$, `findPeaksWithFallback` returns immediately.
* Any valid secondary $p95$ (suggestive) peaks elsewhere on that same chromosome are **completely discarded**.

---

### 3.5. Flag Naming Inversion for GATK Hard Filters
* **File:** `cmd/root.go:517-518`
* **Status:** 🟡 **CLI Usability Bug**

```go
rootCmd.Flags().Float64("min-SOR-SNP", 3.0, "StrandOddsRatio SNPs")
rootCmd.Flags().Float64("min-FS-SNP", 60.0, "FisherStrand SNPs")
```
* In GATK best practices, FS and SOR are **maximum allowed thresholds** ($FS < 60.0, SOR < 3.0$). Variants with $FS > 60$ or $SOR > 3$ are rejected.
* In GoBSAseq's code: `!atMost("FS", hfcfg.SNP_FS_Max)`.
* Calling the flag `--min-FS-SNP` when it enforces a maximum upper bound contradicts standard bioinformatics conventions (and contrasts with the correctly named `--max-FS-INDEL`).

---

### 3.6. Memory Scalability on Whole-Genome Datasets
* **Files:** `stats/stats.go:230`, `stats/smoothing.go:179`
* **Status:** 🟡 **Performance Bottleneck**
* Slices of `BSAstats` and `SmoothedStats` are held fully in memory across the entire genome.
* For crops with large genomes (e.g. 16 Gb wheat, 2.5 Gb maize) yielding 20M+ variants, keeping multiple deep structs in memory simultaneously produces high GC pause times and risks Out-Of-Memory (OOM) termination.
* **Fix:** Stream and write results chromosome-by-chromosome.

---

## 4. Prioritized Actionable Recommendations

### Priority 1: Correct the Smoothing–Null Model Calibration
* **Approach A (Empirical $G'$ Calibration - Magwene et al. 2011):**
  Instead of simulating single markers, fit the empirical null distribution directly to the smoothed $G'$ track across the genome (excluding peaks), or simulate a window of $M$ markers matching the local marker density.
* **Approach B (Analytic Smoothing Adjustment):**
  Compute the effective marker count $M_{eff} = \frac{(\sum w_i)^2}{\sum w_i^2}$. Scale the variance of smoothed statistics by $1 / M_{eff}$ before computing thresholds.

### Priority 2: Fix $BC_n F_2$ Pool AF Sampling
Replace the trimmed logic in `stats/thresholds.go:125-139` with the full trinomial genotypic distribution for selfed backcross generations:
```go
// For BC1F2: P(HH) = 5/8, P(HL) = 2/8, P(LL) = 1/8 (if High recurrent)
// Sample individual dosages (0, 1, or 2 copies) directly per individual in bulkSize
```

### Priority 3: Eliminate Interactive Blocking
* Support a `--non-interactive` / `--batch` flag (or detect whether `os.Stdin` is a TTY via `term.IsTerminal`).
* Recognize `"None"` as a valid sample specifier to skip parent selection without prompting.

### Priority 4: Enforce `MinQTLWidth` and `MinQTLMarkers`
In `stats/detect.go:591`, filter consolidated QTLs by minimum base-pair span (`cfg.MinQTLWidth`) and count of supporting markers (`cfg.MinQTLMarkers`).

### Priority 5: Switch Plot X-Axis to Value Axis
In `plots/plots.go:365-373`, configure the X-axis as `Type: "value"`, passing coordinates as `[pos, value]` 2D tuples rather than category strings.

---

## Conclusion
GoBSAseq possesses an impressive, high-performance Go architecture with multiallelic splitting and caller-aware filtering. However, its QTL detection pipeline currently suffers from a fundamental calibration mismatch between smoothed signals and unsmoothed null simulations, alongside multiple testing and CLI interaction issues. Resolving these items will elevate GoBSAseq into a robust, publication-grade genomic tool.
