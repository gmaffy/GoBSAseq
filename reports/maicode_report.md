# GoBSAseq code review report

## Scope reviewed
I reviewed the core implementation of the BSA-seq pipeline in the following files:
- [filter/filter.go](filter/filter.go)
- [stats/stats.go](stats/stats.go)
- [stats/smoothing.go](stats/smoothing.go)
- [stats/thresholds.go](stats/thresholds.go)
- [stats/detect.go](stats/detect.go)
- [stats/brm.go](stats/brm.go)
- [plots/plots.go](plots/plots.go)
- [run/run.go](run/run.go)
- [utils/config.go](utils/config.go)
- [filter/filter_test.go](filter/filter_test.go)

I also verified the current repository health by running:
- `go test ./...` → passed with exit code 0

---

## Executive summary
Overall, this is a strong and thoughtful BSA-seq prototype. The pipeline is structurally sound, the statistical workflow is mostly coherent, and the implementation shows clear attention to practical issues such as malformed AD/RO fields, multi-allelic variants, and pooled-sample filtering. The code is not merely a toy script; it is a serious attempt to build a reproducible analysis workflow in Go.

That said, the pipeline is not yet at a “10/10” level for publication-quality or production-scale BSA-seq analysis. The biggest limitations are not syntax or plumbing; they are biological and statistical assumptions. The methods are broadly sensible, but some of the current choices are still simplifications that may reduce accuracy in realistic datasets, especially when depth is uneven, parents are not perfectly homozygous, or marker density is low.

### Overall score: 8.3/10
- Correctness of implementation: 8.7/10
- Biological realism: 7.8/10
- Usability and reproducibility: 8.0/10

---

## What the pipeline does well

### 1. The overall workflow is well organized
The pipeline follows a sensible BSA-seq analytical sequence:
1. Determine analysis mode from user input/sample roles
2. Hard filter variants
3. Compute raw statistics
4. Smooth statistics
5. Estimate thresholds
6. Detect QTLs
7. Merge QTLs/BRM intervals
8. Plot outputs

This is the correct high-level architecture for a BSA-seq pipeline and the orchestration in [run/run.go](run/run.go) is clear and easy to follow.

### 2. The hard-filtering logic is thoughtful and adapted to pooled data
The filtering code in [filter/filter.go](filter/filter.go) is one of the strongest parts of the project. It does not blindly apply GATK best-practice filters to pooled BSA-seq data, which is important because many of those filters are tuned for single-sample variant calling and can discard real segregating alleles in bulked pools.

The introduction of a “light filter” mode is especially valuable. That is biologically sensible for BSA-seq, where read imbalance and allele-frequency skew are expected and should not be mistaken for technical artifact in the same way they would be in a standard germline variant-calling context.

### 3. The code handles real-world VCF complications well
The implementation shows good engineering judgment in handling:
- malformed AD/RO fields without panicking
- multi-allelic loci
- missing or low-information allele depth data
- coordinate-sorted output preservation during parallelism

This is above average for a bioinformatics pipeline and reflects real experience with messy genomics data.

### 4. The statistics are standard and broadly appropriate
The raw statistics in [stats/stats.go](stats/stats.go) are familiar and appropriate for BSA-seq:
- SNP index / allele-frequency based measures
- ΔSI
- G-statistic
- LOD
- Bayes factor

These are all reasonable choices for detecting segregation distortion and allele-frequency shifts between contrasting bulks.

### 5. The smoothing and thresholding approach is modern and strong
The Gaussian-kernel smoothing in [stats/smoothing.go](stats/smoothing.go) is a good choice for signal denoising, and the Monte Carlo-based thresholding in [stats/thresholds.go](stats/thresholds.go) is a meaningful improvement over purely analytical thresholding. The use of a two-stage null model is especially commendable and is biologically more realistic than a naive binomial model that ignores the finite number of individuals in each bulk.

---

## Major correctness strengths
The code appears correct in several important ways:
- The flow from hard filter → statistic calculation → smoothing → thresholding → QTL detection is logically consistent.
- The use of allele depths rather than genotype calls alone to infer parent polarity is a good improvement over simplistic genotype-only logic.
- The handling of multi-allelic sites is better than many pipelines that simply discard them.
- The implementation is not fragile in the face of common malformed VCF fields.
- The test suite passes, which is important evidence that the current implementation is at least internally consistent.

---

## Main concerns and areas for improvement

### 1. The biological model is still simplified
The biggest biological limitation is that the statistics are still fairly simplified relative to real segregating populations.

The current implementation assumes that the null expectation is governed mainly by the chosen population structure, such as F2, F3, RIL, or backcross classes. That is reasonable as a first approximation, but in reality BSA-seq data can be affected by:
- residual heterozygosity in parents
- contamination in pools
- uneven sequencing depth across samples
- allele-specific mapping bias
- repetitive or low-complexity regions
- structural variants and CNVs
- variable recombination rate across chromosomes

This means the pipeline will likely work well in clean synthetic or well-controlled datasets, but its biological realism may be lower in messy real-world datasets.

### 2. Parent polarity logic is sensible but still heuristic
The code in [stats/stats.go](stats/stats.go) tries to infer which allele corresponds to the “high” phenotype using allele-depth information. That is a good decision. However, the rule is still fairly simple and deterministic. For real data, a more robust approach would incorporate:
- confidence scores for parent allele assignment
- uncertainty when parent is not fully homozygous
- optional use of a probabilistic model rather than a hard classification

This is not a fatal flaw, but it is one of the main places where the workflow could be made more biologically faithful.

### 3. Thresholding is good but still not fully tailored to the data
The Monte Carlo thresholding in [stats/thresholds.go](stats/thresholds.go) is a major strength, but the current implementation still uses mostly global thresholds and does not fully account for local genomic context.

In practice, QTL detection benefits from considering:
- chromosome-specific marker density
- local recombination structure
- multiple testing burden
- locus-specific sequencing depth and genotype confidence

A purely global threshold can still be useful, but it is not the most statistically refined approach for a pipeline intended to produce publication-ready QTL calls.

### 4. Peak detection is useful but still somewhat coarse
The QTL detection logic in [stats/detect.go](stats/detect.go) is reasonable, but it is still fairly simple. It will detect broad peaks and intervals, but it may not optimally distinguish:
- adjacent but distinct QTLs
- weak broad signals from strong narrow peaks
- regions affected by sparse marker coverage

The BRM block logic in [stats/brm.go](stats/brm.go) is a good addition, but integration between Z-score peaks and BRM blocks is still somewhat heuristic.

### 5. The test coverage is adequate but not yet comprehensive
The tests in [filter/filter_test.go](filter/filter_test.go) are useful and cover a couple of important behaviors. However, the repository still lacks deeper tests for:
- raw statistic correctness
- smoothing outputs
- threshold calculation consistency
- QTL detection behavior
- end-to-end synthetic data analysis

For a pipeline of this complexity, more targeted tests would substantially improve confidence.

---

## What would be needed to reach 10/10
If the goal is a truly top-tier BSA-seq pipeline, I would prioritize the following improvements:

### High-priority improvements
1. Add explicit, validated analysis-mode selection
   - The current implementation infers mode from sample roles, which is convenient, but a clear user-facing mode flag would make the pipeline easier to reason about and less error-prone.

2. Add a probabilistic parent-polarity model
   - Replace the current hard allele polarity decision with a probability-based or confidence-weighted approach.

3. Improve the null model for thresholds
   - The current Monte Carlo approach is already better than naive methods, but it could be expanded to include overdispersion and technical variance more explicitly.

4. Add local and multiple-testing-aware QTL calling
   - Use chromosome-specific or local empirical thresholds and account for the number of tests or genomic regions assessed.

5. Improve robustness to low-depth and low-quality data
   - Add more graceful handling of sparse coverage, low-confidence genotype calls, and non-standard INFO/FORMAT fields.

### Medium-priority improvements
6. Expand unit and integration tests
   - Add synthetic VCF-based regression tests for the full workflow.

7. Improve documentation and example datasets
   - A small test dataset and documented expected outputs would make the pipeline far easier to evaluate and trust.

8. Add support for more population structures and ploidy systems
   - The current framework is strong for common diploid segregating populations, but broader biological use cases would benefit from a more general model.

---

## Bottom line
This is a good and credible BSA-seq analysis pipeline with a solid software foundation and a biologically sensible statistical workflow. It already does many things right, and the current implementation is certainly suitable for exploratory analysis and many standard BSA-seq use cases.

However, if the goal is a 10/10 pipeline that is both statistically rigorous and biologically faithful across diverse datasets, the project would benefit from more explicit probabilistic modeling, more adaptive thresholding, and broader validation. The code is already strong enough to be useful; it now needs refinement to become truly state-of-the-art.

### Final verdict
- Strong prototype: yes
- Production-ready for many use cases: probably
- Publication-grade without further refinement: not yet
- Likely score for a first serious implementation: 8.3/10
