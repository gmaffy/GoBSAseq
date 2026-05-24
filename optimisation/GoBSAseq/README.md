# GoBSAseq

GoBSAseq is a Cobra-based command line program for bulked segregant analysis sequencing (BSA-seq). It reads a VCF containing parent and/or bulk samples, calculates marker-level BSA statistics, smooths them across genomic windows, simulates significance thresholds, calls QTL intervals, and renders interactive HTML plots with go-echarts.

## What It Does

The pipeline supports these analysis layouts:

- Two parents + two bulks: resistant/high parent, susceptible/low parent, high bulk, and low bulk.
- Bulks only: high bulk and low bulk, without parent samples.
- Two parents + one bulk: parent samples plus either a high or low bulk.
- One parent + one bulk: one parent sample and one bulk sample.

When both parents are provided, the susceptible/low parent is treated as the reference background and the resistant/high parent allele is tracked. If parents are not provided, the VCF ALT allele is tracked.

For each informative marker, GoBSAseq calculates:

- SNP-index for each bulk
- Delta SNP-index
- G-statistic
- Euclidean distance
- LOD score
- Beta-binomial Bayes factor
- BRM, block regression mapping

It then applies Gaussian kernel smoothing, simulates null thresholds such as the 95th and 99th percentiles, identifies QTL intervals where smoothed statistics cross thresholds, normalises statistics with robust z-scores, and calls final QTLs from z-score threshold intersections.

## Installation

From the project directory:

```powershell
go build .
```

Or run directly:

```powershell
go run . --help
```

## Input Requirements

The primary input is a VCF file with sample columns and genotype fields containing:

- `GT`: genotype, for example `0/0`, `0/1`, `1/1`
- `AD`: allele depths, for example `12,28`
- `DP`: read depth

The VCF header must contain a standard sample header line:

```text
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  sample1  sample2 ...
```

SNP hard filters are applied from VCF `QUAL` and common INFO fields such as `QD`, `SOR`, `FS`, `MQ`, `MQRankSum`, and `ReadPosRankSum`. Missing INFO fields are allowed.

## Usage

### Interactive VCF-Only Mode

If you provide only a VCF, GoBSAseq reads the VCF header and prompts you to choose parents and bulks from the available samples. Each prompt includes `None`, so missing parents or bulks can be skipped.

```powershell
go run . input.vcf
```

Equivalent:

```powershell
go run . --variant input.vcf
```

You will be prompted for:

- resistant/high parent
- susceptible/low parent
- high bulk
- low bulk

At least one bulk must be selected.

### Two Parents + Two Bulks

```powershell
go run . `
  --variant input.vcf `
  --parents ResistantParent,SusceptibleParent `
  --bulks HighBulk,LowBulk `
  --parents-depth 5,5 `
  --bulks-depth 40,40 `
  --bulk-sizes 20,20 `
  --window-size 2000000 `
  --step-size 100000 `
  --rep 1000 `
  --alpha 0.05,0.01 `
  --out results
```

### Bulks Only

```powershell
go run . `
  --variant input.vcf `
  --bulks HighBulk,LowBulk `
  --bulks-depth 40,40 `
  --bulk-sizes 20,20 `
  --out results_bulks_only
```

### Two Parents + One Bulk

```powershell
go run . `
  --variant input.vcf `
  --parents ResistantParent,SusceptibleParent `
  --bulks HighBulk `
  --parents-depth 5,5 `
  --bulks-depth 40 `
  --bulk-sizes 20 `
  --out results_one_bulk
```

## Important Flags

| Flag | Meaning | Default |
| --- | --- | --- |
| `--variant`, `-V` | Input VCF file | required unless positional VCF is used |
| `--parents`, `-P` | Parent sample names, high/resistant then low/susceptible | none |
| `--bulks`, `-B` | Bulk sample names, high then low | none |
| `--parents-depth`, `-p` | Minimum parent depths | `5,5` |
| `--bulks-depth`, `-b` | Minimum bulk depths | `40,40` |
| `--bulk-sizes`, `-S` | Number of individuals in bulks | `20,20` |
| `--window-size`, `-w` | Smoothing window size in bp | `2000000` |
| `--step-size`, `-s` | Window step size in bp | `100000` |
| `--rep` | Number of threshold simulations | `1000` |
| `--alpha` | Significance levels | `0.05,0.01` |
| `--min-qtl-length` | Minimum QTL interval length | `100000` |
| `--merge-distance` | Merge nearby QTL intervals within this distance | `500000` |
| `--out`, `-o` | Output directory | current directory |

## Outputs

GoBSAseq writes the following files to the output directory:

| File | Contents |
| --- | --- |
| `variant_statistics.tsv` | Marker-level SNP-index and BSA statistics |
| `smoothed_statistics.tsv` | Smoothed statistics and robust z-scores per genomic window |
| `thresholds.tsv` | Simulated threshold values for each statistic and alpha |
| `qtls_by_simulated_thresholds.tsv` | QTL intervals from simulated threshold crossings |
| `final_qtls_by_robust_z.tsv` | Final QTL intervals from robust z-score thresholds |
| `plots/*.html` | Interactive go-echarts plots for each statistic |

Each statistic gets two plots:

- raw smoothed statistic with simulated threshold lines
- robust z-score plot with `z=2`, `z=-2`, `z=3`, and `z=-3` threshold lines

## Pipeline Stages

The pipeline logs major stages while running:

1. Read configuration and infer analysis mode.
2. Read and filter VCF markers.
3. Calculate marker-level statistics.
4. Smooth statistics across chromosomes.
5. Simulate null thresholds.
6. Normalise smoothed statistics with robust z-scores.
7. Detect QTL intervals.
8. Write result tables.
9. Render go-echarts plots.

## Core Functions

### CLI and Sample Selection

- `rootCmd` Run handler in `cmd/root.go`: Cobra command entry point. Reads flags, optionally prompts for VCF samples, builds analysis and filter configs, and starts the pipeline.
- `promptSamplesFromVCF`: Used when only a VCF is provided. It asks the user to select high/low parents and high/low bulks from the VCF header, with `None` allowed.
- `readVCFSamples`: Reads the VCF `#CHROM` header and extracts sample names.
- `splitCSV`: Parses comma-separated sample flags while ignoring empty values and `None`.

### Pipeline Orchestration

- `run.Run`: Main pipeline coordinator. It validates config, creates output folders, logs each stage, runs statistics, smoothing, simulation, QTL calling, table writing, and plotting.
- `inferMode`: Determines whether the run is two-parent/two-bulk, bulks-only, two-parent/one-bulk, or one-parent/one-bulk based on selected samples.

### VCF Parsing and Marker Filtering

- `readMarkers`: Streams the VCF, extracts sample depths from `GT:AD:DP`, applies hard filters, identifies informative markers, orients allele depths, and creates marker records.
- `parseSample`: Converts a sample FORMAT field into reference depth, alternate depth, total depth, and genotype.
- `homozygousAllele`: Identifies homozygous parent genotypes so parent-informative SNPs can be selected.
- `passesHardFilters`: Applies SNP and INDEL hard-filter settings from `utils.HardFilterConfig`.

### Marker Statistics

- `calcTwoBulkStats`: Calculates high/low SNP-index, delta SNP-index, G-statistic, Euclidean distance, LOD, beta-binomial BF, and BRM seed values for two-bulk comparisons.
- `calcOneBulkStats`: Calculates one-bulk equivalents against an expected allele frequency of 0.5 where paired-bulk statistics are not available.
- `gStatistic`: Computes the two-bulk likelihood-ratio G-statistic from allele counts.
- `betaBinomialBF`: Estimates support for separate bulk allele frequencies versus a shared allele frequency.

### Smoothing and BRM

- `smoothMarkers`: Creates sliding genomic windows and smooths each statistic with Gaussian kernel weights.
- `kernelMean`: Weighted mean smoother used for most statistics.
- `localLinear`: Local linear regression smoother used for BRM and marker-level BRM assignment.
- `assignBRM`: Applies block regression mapping values to markers before window smoothing.

### Thresholds and QTL Detection

- `simulateThresholds`: Simulates null data from observed depths, recomputes smoothed statistics, and estimates percentile thresholds for each statistic and alpha level.
- `attachRobustZ`: Converts smoothed values to robust z-scores using median and MAD scaling.
- `detectThresholdQTLs`: Calls QTL intervals where smoothed statistics cross simulated thresholds.
- `detectZQTLs`: Calls final QTL intervals where robust z-scores cross absolute z thresholds of 2 or 3.
- `detectIntervals`: Shared interval builder that merges nearby significant windows and records peak positions.

### Output and Plotting

- `writeMarkers`, `writeSmooth`, `writeThresholds`, and `writeQTLs`: Write tab-delimited result tables.
- `plotAll`: Creates HTML plots for every statistic.
- `plotStat`: Builds individual go-echarts line plots by chromosome, with simulated or z-score thresholds.

## Notes

- For two-parent analyses, parent genotypes must be homozygous and different for a marker to be informative.
- For bulks-only analyses, the VCF ALT allele is tracked because parent origin is unknown.
- For one-bulk analyses, delta-like statistics are measured against an expected allele frequency of `0.5`.
- Increasing `--rep` improves threshold stability but increases runtime.
- Window size and step size should be chosen according to marker density and genome size.
