# GoBSAseq

**A high-performance pipeline for Bulk Segregant Analysis (BSAseq) implemented in Go**

GoBSAseq is a modern, efficient implementation of BSAseq analysis designed to identify genomic regions associated with phenotypic traits through statistical analysis of segregating populations. The pipeline applies hard filtering to VCF files, computes allele frequency statistics, performs spatial smoothing, and detects QTL intervals using multiple statistical methods (Z-score thresholding and Bayesian regression model—BRM).

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
  - [Basic Command](#basic-command)
  - [CLI Arguments](#cli-arguments)
  - [Analysis Modes](#analysis-modes)
- [Input Files](#input-files)
  - [VCF Format Requirements](#vcf-format-requirements)
  - [Sample Naming](#sample-naming)
- [Output Files](#output-files)
- [Algorithm Overview](#algorithm-overview)
  - [Pipeline Steps](#pipeline-steps)
  - [Statistical Methods](#statistical-methods)
- [Population Structures](#population-structures)
- [Filtering Thresholds](#filtering-thresholds)
  - [SNP Filters](#snp-filters)
  - [INDEL Filters](#indel-filters)
- [Gene Space Analysis](#gene-space-analysis)
- [Examples](#examples)
  - [Two-Parent Two-Bulk Analysis](#two-parent-two-bulk-analysis)
  - [High-Parent High-Bulk Analysis](#high-parent-high-bulk-analysis)
- [Output Format Reference](#output-format-reference)
- [Troubleshooting](#troubleshooting)
- [Development](#development)
- [References](#references)

---

## Features

- **Modular pipeline**: Hard filtering → Raw statistics → Smoothing & normalization → Threshold calculation → QTL detection → Plotting
- **Multiple analysis modes**: Support for 2-parent 2-bulk (2p2b), 2-parent high-bulk (2phb), 2-parent low-bulk (2plb), high-parent 2-bulk (hp2b), low-parent 2-bulk (lp2b), single-parent single-bulk (hphb, hplb, lphb, lplb), and bulks-only (2b) designs
- **Robust statistics**: Gaussian kernel smoothing, median-absolute-deviation (MAD) Z-scores, Stouffer composite Z-scores, and Bayesian regression model (BRM) block detection
- **Adaptive thresholds**: Empirical per-variant significance thresholds based on population structure
- **VCF compression and indexing**: BGZF compression with tabix indexing for fast variant lookup
- **Interactive and scripted modes**: Sample selection via prompts or command-line flags
- **HTML visualization**: Interactive scatter plots and overlay charts via echarts
- **Gene space analysis**: Annotation of QTL regions using snpEff and custom GFF annotations (optional)

---

## Quick Start

1. **Clone and build**:
   ```bash
   cd GoBSAseq
   go build -o gobsaseq
   ```

2. **Run a 2-parent 2-bulk analysis**:
   ```bash
   ./gobsaseq \
     -V my_variants.vcf.gz \
     -P parent1,parent2 \
     -B bulk_high,bulk_low \
     -p 10,10 \
     -b 100,100 \
     -S 500,500 \
     -w 2000000 \
     -s 100000 \
     -o results
   ```

3. **Check output**:
   ```bash
   ls results/goBSAseqResults/<timestamp>/
   # plots/, qtls/, stats/
   ```

---

## Installation

### Requirements

- **Go 1.26.3** or later ([download](https://go.dev/dl/))
- **Standard build tools**: `make`, `git` (optional)
- **Python 3** (optional, for gene space analysis: snpEff)

### Build from Source

```bash
git clone https://github.com/gmaffy/GoBSAseq.git
cd GoBSAseq
go mod download
go build -o gobsaseq
```

### Verify Installation

```bash
./gobsaseq --help
# Shows usage and available flags
```

---

## Usage

### Basic Command

```bash
GoBSAseq -V <vcf> -P <parents> -B <bulks> [options]
```

### CLI Arguments

#### Required Arguments

| Flag | Long Form | Type | Description |
|------|-----------|------|-------------|
| `-V` | `--variant` | string | Path to VCF file (`.vcf` or `.vcf.gz`) |
| `-P` | `--parents` | string | Parent sample names, comma-separated (e.g., `parent1,parent2` or `None` to skip) |
| `-B` | `--bulks` | string | Bulk sample names, comma-separated (e.g., `bulk_high,bulk_low` or `None` to skip) |

#### Sequencing Parameters

| Flag | Long Form | Type | Default | Description |
|------|-----------|------|---------|-------------|
| `-p` | `--parents-depth` | string | `5,5` | Minimum read depth for parents (high,low) |
| `-b` | `--bulks-depth` | string | `40,40` | Minimum read depth for bulks (high,low) |
| `-S` | `--bulk-sizes` | string | `20,20` | Number of individuals in bulks (high,low) |

#### Smoothing Parameters

| Flag | Long Form | Type | Default | Description |
|------|-----------|------|---------|-------------|
| `-w` | `--window-size` | int64 | `2000000` | Gaussian kernel σ (bp); actual window size is 2×σ |
| `-s` | `--step-size` | int64 | `100000` | Sliding window step size (bp) |

#### Statistical Parameters

| Flag | Long Form | Type | Default | Description |
|------|-----------|------|---------|-------------|
| `-m` | `--population` | string | `F2` | Population structure: `F2`, `F3`, `RIL`, `BC1H`, `BC1L`, `BC2H`, `BC2L` |
| `--brm-alpha` | N/A | float64 | `0.05` | Significance level for BRM threshold calculation |
| `--rep` | N/A | int | `1000` | Number of simulations for threshold calculation |
| `--min-qtl-length` | N/A | int64 | (internal) | Minimum QTL interval length |
| `--merge-distance` | N/A | int64 | (internal) | Distance to merge nearby QTL intervals |

#### Variant Filtering (SNP)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--min-QD-SNP` | float64 | `2.0` | Minimum quality-by-depth (QualByDepth) |
| `--min-QUAL-SNP` | float64 | `30.0` | Minimum variant quality score |
| `--min-MQ-SNP` | float64 | `40.0` | Minimum RMS mapping quality |
| `--min-SOR-SNP` | float64 | `3.0` | Minimum strand odds ratio (higher = allow fewer strand bias) |
| `--min-FS-SNP` | float64 | `60.0` | Maximum FisherStrand bias (lower = more stringent) |
| `--min-MQRankSum-SNP` | float64 | `-12.5` | Minimum mapping quality rank sum test |
| `--min-ReadPosRankSum-SNP` | float64 | `-8.0` | Minimum read position rank sum |

#### Variant Filtering (INDEL)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--min-QD-INDEL` | float64 | `2.0` | Minimum quality-by-depth |
| `--min-QUAL-INDEL` | float64 | `30.0` | Minimum variant quality score |
| `--max-FS-INDEL` | float64 | `200.0` | Maximum FisherStrand bias |
| `--max-SOR-INDEL` | float64 | `10.0` | Maximum strand odds ratio |
| `--min-ReadPosRankSum-INDEL` | float64 | `-20.0` | Minimum read position rank sum |

#### Gene Space Analysis (Optional)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--snpEffDB` | string | `` | snpEff database ID (e.g., `GRCh38.p13`) |
| `--gff` | string | `` | GFF3 annotation file path |
| `--cds` | string | `` | CDS FASTA file path |
| `-r` / `--reference` | string | `` | Reference genome FASTA path |
| `--protein` | string | `` | Protein FASTA file path |
| `--gene-descriptions` | string | `` | Gene descriptions TSV (tab-separated: gene_id, description) |
| `--prg` | string | `` | PRG (pangenome graph) blast file path |

#### Output

| Flag | Long Form | Type | Default | Description |
|------|-----------|------|---------|-------------|
| `-o` | `--out` | string | `.` | Output directory |

---

## Analysis Modes

GoBSAseq detects the analysis mode automatically based on which samples are provided. Each mode has its own filtering logic and statistical calculations:

### 2p2b (Two Parent, Two Bulk)
- **When**: Both parents and both bulks are provided
- **Filtering**: Parents must be homozygous and divergent; bulks must show segregation
- **Stats**: ΔSI (allele frequency difference), G-statistic, LOD, Bayes factor
- **Use case**: Classic F2, RIL, or backcross design with distinct high and low phenotype bulks

### 2phb (Two Parent, High Bulk)
- **When**: Both parents + high bulk (no low bulk)
- **Stats**: Same as 2p2b but only for high bulk; one-bulk statistics also computed
- **Use case**: Phenotype extreme available; other bulk discarded or unavailable

### 2plb (Two Parent, Low Bulk)
- **When**: Both parents + low bulk (no high bulk)
- **Stats**: Same as 2phb but focusing on low bulk

### hp2b (High Parent, Two Bulk)
- **When**: Only high parent + both bulks
- **Use case**: Low parent genotype unknown or unavailable

### lp2b (Low Parent, Two Bulk)
- **When**: Only low parent + both bulks

### hphb (High Parent, High Bulk)
- **When**: Only high parent + high bulk
- **Stats**: Single-parent single-bulk tests (AF deviation, G-stat, LOD, BB log-BF)

### hplb (High Parent, Low Bulk)
- **When**: Only high parent + low bulk

### lphb (Low Parent, High Bulk)
- **When**: Only low parent + high bulk

### lplb (Low Parent, Low Bulk)
- **When**: Only low parent + low bulk

### 2b (Bulks Only)
- **When**: Only high and low bulks (no parents)
- **Filtering**: Bulks segregate (not fixed identical genotypes)
- **Stats**: Two-bulk statistics only

---

## Input Files

### VCF Format Requirements

GoBSAseq expects a **multi-sample VCF** (version 4.1 or 4.2) with:

1. **Sample columns**: One column per sample in the order they appear in the VCF header
2. **Genotype field (GT)**: Required
   - Format: `ref/alt` or `ref|alt` (phased or unphased)
   - Allowed alleles: `0` (ref), `1` (alt), or `.` (missing)
3. **Allele depth field (AD)**: Required for statistics
   - Format: `ref_depth,alt1_depth[,alt2_depth,...]`
4. **Read depth field (DP)**: Required for filtering
5. **Quality fields**: QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum (in INFO)
   - Used by hard filter (SNP/INDEL specific thresholds)

#### Example VCF Header

```vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand bias">
...
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	parent1	parent2	bulk_high	bulk_low
chr1	1000	.	A	T	60	PASS	QD=20.0;FS=0.0;SOR=0.8;MQ=60	GT:AD:DP	0/0:50,0:50	1/1:0,45:45	0/1:25,25:50	0/1:30,14:44
```

### Sample Naming

- Sample names are taken from the VCF header
- Provide sample names to GoBSAseq using comma-separated lists
- Use `None` as a placeholder to skip a role (e.g., `-P None -B bulk_high,bulk_low` for bulks-only analysis)
- If no names provided on command-line, GoBSAseq will prompt interactively

---

## Output Files

Results are written to `<out>/goBSAseqResults/<timestamp>/` with the following directory structure:

```
goBSAseqResults/
├── 10_06_2026_14_30_45/         # Timestamp: day_month_year_hour_minute_second
│   ├── stats/
│   │   ├── GoBSAseq.2p2b.raw.tsv
│   │   ├── GoBSAseq.2p2b.smoothed_and_normalised.tsv
│   │   ├── GoBSAseq.2p2b.brm_blocks.tsv
│   │   ├── GoBSAseq.2p2b.hardfiltered.vcf.gz
│   │   ├── GoBSAseq.2p2b.hardfiltered.vcf.gz.tbi
│   │   └── GoBSAseq.2p2b.lowqual.vcf.gz
│   ├── qtls/
│   │   ├── GoBSAseq.2p2b.qtls.tsv
│   │   └── GoBSAseq.2p2b.merged_qtls.tsv
│   ├── plots/
│   │   ├── GoBSAseq.2p2b.individual_plots.html
│   │   ├── GoBSAseq.2p2b.robust_z_overlay.html
│   │   └── GoBSAseq.2p2b.composite_signal.html
│   └── genespace/
│       ├── GoBSAseq.2p2b.genes_in_qtls.tsv    (if enabled)
│       └── GoBSAseq.2p2b.qtl_summary.txt      (if enabled)
```

### Key Output Files

#### Raw Statistics (`*.raw.tsv`)
Per-variant raw statistics before smoothing. Columns:
- `CHROM`, `POS`, `REF`, `ALT`: variant location and alleles
- `HighParGT`, `LowParGT`, `HighBulkGT`, `LowBulkGT`: genotypes (if applicable)
- `HighBulkAD`, `LowBulkAD`: allele depths (if bulk present)
- `HighSI`, `LowSI`: **selection index** (allele frequency in high/low bulk)
- `DeltaSI`: difference in allele frequencies
- `Gstat`, `ED4`, `LOD`, `BBLogBF`: statistical test values
- `Depth`: combined read depth

#### Smoothed & Normalised (`*.smoothed_and_normalised.tsv`)
Per-variant Gaussian smoothed values and robust Z-scores. Columns:
- `Sm_*`: smoothed raw values (Gaussian kernel, bandwidth = window-size / 2)
- `Z_*`: robust Z-scores (normalized by median ± MAD, global normalization)
- `CompositeZ`: Stouffer composite Z-score (average of all available Z-scores)
- `MaxAbsZ`: maximum absolute Z-score across all statistics

#### QTLs (`*.qtls.tsv`)
Z-score-based QTL regions (compos Z peaks ≥ 3.0 per chromosome). Columns:
- `CHROM`, `START`, `STOP`: QTL interval
- `PEAK`: z-score at the most extreme position

#### BRM Blocks (`*.brm_blocks.tsv`)
Bayesian regression model (BRM) detected blocks. Columns:
- `CHROM`, `START`, `STOP`: interval
- `PEAK_POS`: genomic position of peak
- `PEAK`: peak Δ SI or AF deviation
- `THRESHOLD`: BRM significance threshold at peak

#### Merged QTLs (`*.merged_qtls.tsv`)
Union of Z-score and BRM intervals per chromosome. Columns:
- `CHROM`, `START`, `STOP`: merged interval
- `SOURCE`: `"ZScore"`, `"BRM"`, or `"ZScore+BRM"`
- `Z_PEAK`, `BRM_PEAK`, `BRM_THRESHOLD`: respective peak values (NA if method not present)

#### Filtered VCF (`*.hardfiltered.vcf.gz` + `.tbi`)
- BGZF-compressed VCF with only variants passing hard filters
- Tabix index (`.tbi`) for fast region lookups
- Includes only samples relevant to the analysis mode

#### Plots (`*.html`)
Interactive echarts visualizations:
- `individual_plots.html`: Per-statistic line charts with threshold overlays and BRM shading
- `robust_z_overlay.html`: All Z-scores overlaid per chromosome with ±2/±3 reference lines
- `composite_signal.html`: CompositeZ and MaxAbsZ together; easier regional scanning

---

## Algorithm Overview

### Pipeline Steps

```
1. Hard Filter (filter.HardFilterVcf)
   • Remove variants with <1 real ALT
   • SNP/INDEL-specific QC thresholds (QD, FS, SOR, MQ, etc.)
   • Sample genotype filtering (homozygosity checks, depth filters)
   → Output: GoBSAseq.*.hardfiltered.vcf.gz (indexed)

2. Raw Statistics (stats.RawStats)
   • Per-variant allele depths and frequency calculations
   • SI (selection index) = freq(high-allele) in high/low bulk
   • Statistical tests: G, LOD, Bayes factor
   → Output: GoBSAseq.*.raw.tsv

3. Smoothing & Normalisation (stats.SmoothAndNormalise)
   • Gaussian kernel smoothing per chromosome
   • Depth-weighted averaging (w = depth × exp(-0.5 × d²/σ²))
   • Robust Z-score normalisation (global, median±MAD)
   • Stouffer composite Z-score
   → Output: GoBSAseq.*.smoothed_and_normalised.tsv

4. Threshold Calculation (stats.CalculateThresholds)
   • Per-variant empirical p99/p95 thresholds via simulation
   • Depends on population structure (expected allele frequency)
   → Internal use (stored with each variant)

5. BRM Block Detection (stats.RunBRM)
   • Sliding window significance tests
   • Thresholds scale with allele frequency and bulk size
   • Outputs contiguous regions above threshold
   → Output: GoBSAseq.*.brm_blocks.tsv

6. QTL Detection (stats.DetectQTLs)
   • Scan CompositeZ for |z| ≥ 3.0
   • Keep most extreme run per chromosome
   → Output: GoBSAseq.*.qtls.tsv

7. Merge QTLs & BRM (stats.MergeQTLsAndBRM)
   • Union intervals from both methods
   • Label source (ZScore, BRM, or both)
   → Output: GoBSAseq.*.merged_qtls.tsv

8. Plotting (plots.GeneratePlots)
   • Interactive HTML echarts per chromosome
   → Output: *.html

9. Gene Space (stats.GeneSpaceFromMerged)
   • Annotate QTL regions with genes
   • Optional; requires snpEff DB or GFF file
   → Output: GoBSAseq.*.genes_in_qtls.tsv (if enabled)
```

### Statistical Methods

#### Selection Index (SI)
$$\text{SI} = \frac{\text{freq}(\text{high allele})}{\text{freq}(\text{all alleles})}$$

In high bulk: high SI = allele enriched in phenotype extreme
In low bulk: low SI = allele depleted in phenotype extreme

#### Allele Frequency Difference (ΔSI)
$$\Delta\text{SI} = \text{SI}_{\text{high}} - \text{SI}_{\text{low}}$$

Large |ΔSI| suggests a QTL allele.

#### G-Statistic (Likelihood Ratio)
$$G = 2 \sum_i n_i \log\left(\frac{o_i}{e_i}\right)$$

Tests deviation from expected allele frequency ratio under null. Approximately χ²-distributed.

#### LOD Score
$$\text{LOD} = \log_{10}\left(\frac{L(\text{alt})}{L(\text{null})}\right)$$

Likelihood ratio of alternative (independent frequencies per bulk) vs. null (pooled frequency). LOD > 3 = significant, LOD < -2 = significant against.

#### Bayes Factor (Beta-Binomial)
$$\text{BB LogBF} = \log_{10}\left(\frac{P(\text{data}|\text{alt})}{P(\text{data}|\text{null})}\right)$$

Uses Beta(0.5, 0.5) conjugate prior. More conservative than LOD under model misspecification.

#### Robust Z-Score (MAD)
$$Z_i = \frac{x_i - \text{median}(x)}{1.4826 \times \text{median}(|x - \text{median}(x)|)}$$

Resistant to outliers; less affected by single extreme data points than classical Z = (x - μ) / σ.

#### Stouffer Composite Z-Score
$$Z_{\text{composite}} = \frac{\sum_i Z_i}{\sqrt{k}}$$

Averages all available Z-scores per variant, preserving sign (direction).

#### BRM (Bayesian Regression Model)
Two-bulk case:
$$\Delta\text{SI} \sim N(0, \sigma^2),\quad \sigma^2 = u_\alpha \sqrt{\frac{(n_1 + n_2) \cdot p(1-p)}{2 n_1 n_2}}$$

where $u_\alpha$ is the inverse normal CDF at 1 - α/2 (e.g., u₀.₀₅ ≈ 1.96 for p = 0.05).

One-bulk case:
$$\text{AF} \sim N(p_0, \sigma^2),\quad \sigma^2 = u_\alpha \sqrt{\frac{p_0(1-p_0)}{2 n}}$$

Blocks are regions where the statistic exceeds the local threshold.

---

## Population Structures

Expected allele frequencies (p₀) differ by population structure. These affect threshold calculations:

| Structure | p₀ | Use Case |
|-----------|-----|----------|
| **F2** | 0.5 | Classical F2 intercross; all allele frequencies segregate equally |
| **F3** | 0.5 | F3 intercross (same as F2 for allele frequencies) |
| **RIL** | 0.5 | Recombinant inbred lines |
| **BC1H** | 0.75 | One backcross to high parent; 3:1 ratio expected |
| **BC1L** | 0.25 | One backcross to low parent; 1:3 ratio expected |
| **BC2H** | 0.875 | Two backcrosses to high parent; strong bias toward high |
| **BC2L** | 0.125 | Two backcrosses to low parent; strong bias toward low |

Specify with `-m` / `--population` flag.

---

## Filtering Thresholds

### SNP Filters

Hard filters are applied **before** statistical analysis to remove low-quality variants:

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `--min-QD-SNP` | 2.0 | Quality ≥ 2× depth (rules out low-coverage noise) |
| `--min-QUAL-SNP` | 30.0 | Phred quality ≥ 30 (~0.1% error probability) |
| `--min-MQ-SNP` | 40.0 | RMS mapping quality ≥ 40 (high-confidence alignments) |
| `--min-SOR-SNP` | 3.0 | Strand odds ratio ≤ 3 (some strand bias OK; inverted logic) |
| `--min-FS-SNP` | 60.0 | FisherStrand ≤ 60 (moderate strand bias threshold) |
| `--min-MQRankSum-SNP` | -12.5 | RankSum test for mapping quality |
| `--min-ReadPosRankSum-SNP` | -8.0 | RankSum test for read position bias |

### INDEL Filters

INDELs typically tolerate higher strand bias and require slightly different thresholds:

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `--min-QD-INDEL` | 2.0 | Same as SNPs |
| `--min-QUAL-INDEL` | 30.0 | Same as SNPs |
| `--max-FS-INDEL` | 200.0 | Much higher tolerance than SNPs (INDELs accumulate strand bias) |
| `--max-SOR-INDEL` | 10.0 | Also more permissive |
| `--min-ReadPosRankSum-INDEL` | -20.0 | Looser than SNPs |

**Tip**: Defaults are conservative. To be more permissive, increase QD/QUAL thresholds and relax strand-bias filters. To be stricter, lower thresholds.

---

## Gene Space Analysis

If reference annotation files are provided, GoBSAseq can annotate QTL regions with overlapping genes:

```bash
./gobsaseq \
  -V variants.vcf.gz \
  -P parent1,parent2 \
  -B bulk_h,bulk_l \
  --gff annotations.gff3 \
  --protein proteins.fasta \
  --gene-descriptions genes.tsv \
  -o results
```

**Required files**:
- `--gff`: GFF3 file with gene features (CDS / gene spans)
- `--gene-descriptions`: TSV with columns `gene_id` and `description`

**Optional files**:
- `--reference`: Reference genome FASTA (for sequence extraction)
- `--protein`: Protein FASTA (for later analysis)
- `--snpEffDB`: snpEff DB ID (for variant annotation)

snpEff annotation can be slow for large VCFs. Consider running separately and merging results if needed.

---

## Examples

### Two-Parent Two-Bulk Analysis

**Command**:
```bash
./gobsaseq \
  -V data/crosses.vcf.gz \
  -P parent_a,parent_b \
  -B phenotype_high,phenotype_low \
  --parents-depth 10,10 \
  --bulks-depth 100,100 \
  --bulk-sizes 250,250 \
  -w 1000000 \
  -s 50000 \
  -m F2 \
  --brm-alpha 0.05 \
  --rep 1000 \
  -o results
```

**What it does**:
1. Loads VCF with 4 samples: two parents, two bulks (each ≥250 individuals)
2. Hard-filters variants (SNPs: QD≥2, QUAL≥30, FS≤60; INDELs: QD≥2, QUAL≥30, FS≤200)
3. Extracts allele depths per bulk; calculates SI per variant
4. Applies 1 Mb Gaussian window (σ) with 50 kb step
5. Computes threshold via 1000 simulations assuming F2 structure
6. Detects BRM blocks and Z-score QTLs
7. Writes HTML plots and TSV tables

**Output** (`results/goBSAseqResults/<timestamp>/`):
- `stats/GoBSAseq.2p2b.raw.tsv` — Per-variant statistics (HighSI, LowSI, DeltaSI, G, LOD, BBLogBF, Depth)
- `stats/GoBSAseq.2p2b.smoothed_and_normalised.tsv` — Smoothed + Z-scored
- `stats/GoBSAseq.2p2b.brm_blocks.tsv` — BRM-detected blocks
- `qtls/GoBSAseq.2p2b.qtls.tsv` — Z-score QTLs (|z| ≥ 3)
- `qtls/GoBSAseq.2p2b.merged_qtls.tsv` — Combined QTLs from both methods
- `plots/GoBSAseq.2p2b.*.html` — Interactive visualization

**Interpretation**:
- Peaks in `individual_plots.html` near or above p99 thresholds = potential QTL regions
- High `CompositeZ` in `composite_signal.html` = strong signal
- Compare `qtls.tsv` (Z-score) vs. `brm_blocks.tsv` (BRM) for consensus

### High-Parent High-Bulk Analysis

**Command**:
```bash
./gobsaseq \
  -V data/variants.vcf.gz \
  -P parent_high,None \
  -B bulk_high,None \
  --parents-depth 10,10 \
  --bulks-depth 150,150 \
  --bulk-sizes 500,500 \
  -w 2000000 \
  -s 100000 \
  -m F2 \
  -o results
```

**What it does**:
1. Uses only high parent + high bulk
2. Computes one-bulk statistics (AF deviation, G, LOD vs. expected p₀ = 0.5)
3. Same smoothing & detection pipeline
4. Thresholds adapted for single-bulk case

**Output**: Same structure, mode = `hphb`

---

## Output Format Reference

### TSV Column Descriptions

#### `*.raw.tsv`

Two-bulk mode columns:
```
CHROM     Chromosome name
POS       Position (1-based)
REF       Reference allele
ALT       Alternate allele
HighParGT Genotype of high parent (e.g., "[0, 0]" = homozygous ref)
LowParGT  Genotype of low parent
HighBulkGT Genotype of high bulk (e.g., "[0, 1]" = heterozygous)
HighBulkAD Allele depths in high bulk (e.g., "100,50")
HighBulkL  Count of recessive (lower-frequency) allele in high bulk
HighBulkH  Count of dominant (higher-frequency) allele in high bulk
HighSI    Selection index in high bulk: H / (H + L)
LowBulkGT  Genotype of low bulk
LowBulkAD  Allele depths in low bulk
LowBulkL   Count of recessive allele in low bulk
LowBulkH   Count of dominant allele in low bulk
LowSI     Selection index in low bulk
DeltaSI   HighSI - LowSI
Gstat     G-statistic (likelihood ratio)
ED4       Euclidean distance (HighSI - LowSI)⁴
LOD       Log odds ratio (log₁₀)
BBLogBF   Bayes factor beta-binomial (log₁₀)
Depth     Minimum depth across bulks
```

One-bulk mode adds:
```
P0          Expected allele frequency (from population structure)
AFDev       SI - P0 (allele frequency deviation)
Gstat1      One-bulk G-statistic
LOD1        One-bulk LOD
BBLogBF1    One-bulk Bayes factor
```

#### `*.smoothed_and_normalised.tsv`

Extends `raw.tsv` columns with:
```
Sm_HighSI   Gaussian smoothed HighSI
Z_HighSI    Robust Z-score of Sm_HighSI (global normalization)
Sm_LowSI    Smoothed LowSI
Z_LowSI     Robust Z-score of Sm_LowSI
Sm_DeltaSI  Smoothed DeltaSI
Z_DeltaSI   Robust Z-score of Sm_DeltaSI
Sm_Gstat    Smoothed G-statistic
Z_Gstat     Robust Z-score
Sm_ED4      Smoothed ED⁴
Z_ED4       Robust Z-score
Sm_LOD      Smoothed LOD
Z_LOD       Robust Z-score
Sm_BBLogBF  Smoothed BB log-BF
Z_BBLogBF   Robust Z-score
CompositeZ  Stouffer composite Z = (Σ Z_i) / √k (signed, direction = sign of DeltaSI / AFDev)
MaxAbsZ     Maximum |Z_i| across all statistics (unsigned)
```

#### `*.brm_blocks.tsv`

```
CHROM     Chromosome
START     Block start (bp, median between prior variant and current)
STOP      Block end (bp)
PEAK_POS  Genomic position of peak value
PEAK      Peak Δ SI or AF deviation
THRESHOLD Significance threshold at PEAK position
```

#### `*.qtls.tsv` (Z-Score QTLs)

```
CHROM Chromosome
START Most upstream variant position with |z| ≥ 3.0
STOP  Most downstream variant position with |z| ≥ 3.0
PEAK  Signed CompositeZ at the peak
```

#### `*.merged_qtls.tsv`

```
CHROM           Chromosome
START           Union interval start
STOP            Union interval stop
SOURCE          "ZScore" | "BRM" | "ZScore+BRM"
Z_PEAK          Peak CompositeZ from Z-score method (NA if not present)
BRM_PEAK        Peak Δ SI / AFDev from BRM method (NA if not present)
BRM_THRESHOLD   BRM threshold at peak position (NA if not present)
```

---

## Troubleshooting

### VCF Parsing Errors

**Error**: `VCF parse error at line X: bad sample string`

- **Cause**: Malformed genotype field (e.g., `0|` instead of `0/0`)
- **Solution**: Validate VCF with `bcftools check <file.vcf>` or re-generate from variant caller

### Missing Samples in VCF

**Error**: `HIGH PARENT <name> is not part of the VCF sample list`

- **Cause**: Sample name mismatch or typo
- **Solution**: Check VCF header with `bcftools query -l <file.vcf>` and use exact names

### No Variants Pass Filtering

**Error**: Hard filtering complete: 0 variant records read → 0 passed

- **Possible causes**:
  1. VCF is empty or contains only reference blocks (gVCF)
  2. Filtering thresholds too stringent
  3. VCF uses different INFO field names (non-standard)

- **Solution**:
  1. Relax filters: increase `--min-QUAL-SNP`, `--min-QD-SNP`; increase `--max-FS-SNP`, `--max-FS-INDEL`
  2. Check VCF INFO fields: `bcftools query -f '%CHROM\t%POS\t%QD\t%FS\n' file.vcf | head`
  3. Validate genotypes: `bcftools view -g ^m file.vcf | head` (show non-missing genotypes)

### Empty Output TSV Files

**Cause**: All variants removed after hard filter + BSAseq filter

- **Solution**:
  1. Reduce filtering stringency
  2. Verify samples are correctly specified
  3. Check that genotypes are segregating in bulks (not all ref or all alt)

### No QTLs Detected

**Cause**: No variants reach |CompositeZ| ≥ 3.0

- **Possible reasons**:
  1. Low statistical power (small bulk sizes, low depth)
  2. Weak or multiple QTLs (signal spread across genome)
  3. Artifact/noise (not biology)

- **Diagnostics**:
  1. Check smoothed TSV: search for maximum `CompositeZ` and `MaxAbsZ` values
  2. Lower detection threshold in post-processing (use `--brm-alpha 0.10` or visualize manually in HTML plots)
  3. Inspect individual chromosome plots for regional peaks

### HTML Plots Not Rendering

**Cause**: Browser security policy blocks local file echarts

- **Solution**:
  1. Open plots with a local HTTP server: `python3 -m http.server 8000` then `localhost:8000/plots/`
  2. Or use a tool like `live-server` (npm) to serve files locally

### Memory Usage Spike

**Cause**: Large VCF or many variants after filtering

- **Solution**:
  1. Process chromosomes separately (subset VCF by `--regions`)
  2. Increase window size (`-w`) to reduce per-chromosome variant count after smoothing
  3. Run on a machine with more RAM

### BRM Blocks vs. Z-Score QTLs Disagree

This is expected! The methods test different aspects:

- **Z-score (CompositeZ)**: Averaged normalized test statistics; sensitive to consistent (if mild) signals across multiple tests
- **BRM**: Regression model specifically for Δ SI; more parsimonious (fewer parameters)

Consensus intervals (SOURCE = "ZScore+BRM") are most robust.

---

## Development

### Building & Testing

```bash
# Build
go build -o gobsaseq

# Run linters
go vet ./...
golangci-lint run ./...

# Run tests (currently minimal; see "Contributing")
go test ./...

# Format code
gofmt -w .
```

### Project Structure

```
GoBSAseq/
├── main.go           # Entry point
├── cmd/root.go       # CLI setup (Cobra)
├── run/run.go        # Main pipeline orchestration
├── filter/filter.go  # VCF filtering + tabix indexing
├── stats/
│   ├── stats.go      # Raw statistics calculation
│   ├── smoothing.go  # Gaussian smoothing + normalization
│   ├── detect.go     # QTL detection + merging
│   ├── brm.go        # BRM block detection
│   ├── thresholds.go # Empirical threshold calculation
│   └── genespace.go  # Gene annotation (optional)
├── plots/plots.go    # HTML plotting via go-echarts
├── utils/
│   ├── config.go     # Configuration structs
│   └── utils.go      # Utility functions
├── go.mod            # Module definition
└── README.md         # This file
```

### Key Packages

- **brentp/vcfgo**: VCF parsing & writing
- **go-echarts**: Interactive plotting
- **spf13/cobra**: CLI framework
- **gonum**: Numerical utilities

### Contributing

1. Fork and clone the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Add tests: `*_test.go` files in relevant packages
4. Run `go test ./...` and `go vet ./...`
5. Submit a pull request

**High-priority areas for contribution**:
- Unit tests for statistical functions (`stats/stats.go`, `stats/smoothing.go`)
- Additional population structures (custom backcross schemes)
- Performance optimization (parallel per-chromosome processing, memory reduction)
- Extended gene annotation (additional file formats, integration with external DBs)

---

## References

1. **Bulk Segregant Analysis (BSA)**:
   - Young, N. D. (1996). Constructing a binned molecular map using molecular markers. In DNA-based Markers in Plants, pp. 49-60. Springer.
   - Michelmore, R. W., Paran, I., & Kesseli, R. V. (1991). Identification of markers linked to disease-resistance genes by bulked segregant analysis: a rapid method to detect markers in specific genomic regions using pooled DNA samples. PNAS, 88(21), 9828-9832.

2. **Selection Index (SI) & Δ SNPindex**:
   - Magwene, P. M., Broman, K. W., & Scan, J. H. (2011). Inferring the joint evolutionary history of multiple loci. G3: Genes, Genomes, Genetics, 1(5), 417-425.

3. **G-Statistic**:
   - Sokal, R. R., & Rohlf, F. J. (1995). Biometry: The Principles and Practice of Statistics in Biological Research (3rd ed.). Freeman.

4. **LOD Score**:
   - Morton, N. E. (1955). Sequential tests for the detection of linkage. American Journal of Human Genetics, 7(3), 277.

5. **Bayesian Methods (Bayes Factor)**:
   - Kass, R. E., & Raftery, A. E. (1995). Bayes factors. Journal of the American Statistical Association, 90(430), 773-795.

6. **Robust Statistics (MAD)**:
   - Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the median absolute deviation. Journal of the American Statistical Association, 88(424), 1273-1283.

7. **Stouffer's Z**:
   - Stouffer, S. A., et al. (1949). The American Soldier: Adjustment During Army Life (Vol. 1). Princeton University Press.

8. **VCF Format**:
   - Danecek, P., et al. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156-2158.

9. **Go Language**:
   - The Go Programming Language (https://golang.org/)

10. **Bioinformatics Tools**:
    - Htslib: https://github.com/samtools/htslib
    - BCFtools: https://github.com/samtools/bcftools
    - snpEff: http://snpeff.sourceforge.net/

---

## License

GoBSAseq is distributed under the [LICENSE](LICENSE) file in this repository.

---

## Acknowledgments

- Original BSAseq concept and statistical framework
- Contributors: See CONTRIBUTORS file (if available)
- Go community for excellent tooling and libraries

---

## Contact & Support

**Issues & Feature Requests**: [GitHub Issues](https://github.com/gmaffy/GoBSAseq/issues)

**Questions**: Open a discussion or check existing issues for common questions.

---

**Last Updated**: June 2026  
**Version**: 0.1.0
