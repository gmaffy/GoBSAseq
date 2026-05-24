# GoBSAseq Best Of Both

This directory contains a fresh, self-contained GoBSAseq app that combines the strongest parts of the original app with the cleaner CLI and pipeline shape of the optimized app.

The goal is practical: keep the older app's stronger BSA-seq biology and statistical machinery, while using the newer app's simpler orchestration, interactive VCF-only workflow, clearer logging, consistent output tables, and maintainable layout.

## What This Version Combines

From the older app:

- Structured VCF parsing with `vcfgo`
- Hard filtering from VCF `QUAL` and INFO annotations
- BSA-specific informative-marker filtering
- Parent-aware allele orientation
- ED4 rather than plain Euclidean distance
- Tricube spatial smoothing with depth weighting
- Depth-adaptive simulated thresholds
- Robust z-score normalisation with trimmed background estimation
- Consensus QTL calls
- Max absolute z-score composite QTL calls
- BRM-style block detection and high-confidence QTL intersection

From the optimized app:

- Cleaner single pipeline entry point
- Positional VCF support
- Interactive VCF-only sample selection with `None`
- Clear stage logging
- Unified outputs for all supported modes
- Simpler mode inference
- Self-contained README and app structure

## Supported Analysis Modes

The app automatically chooses the analysis path from the selected samples:

- Two parents + two bulks
- Bulks only
- Two parents + one bulk
- One parent + one bulk

When two parents are present, the resistant/high parent allele is tracked against the susceptible/low parent background. Without parents, the real VCF ALT allele is tracked.

## Build

From this directory:

```powershell
go build .
```

Run help:

```powershell
go run . --help
```

## Usage

### Interactive VCF-Only Mode

```powershell
go run . input.vcf
```

or:

```powershell
go run . --variant input.vcf
```

If no `--parents` or `--bulks` are supplied, the app reads the VCF header and prompts for:

- resistant/high parent
- susceptible/low parent
- high bulk
- low bulk

Each prompt includes `0) None`. At least one bulk must be selected.

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
| `--variant`, `-V` | Input VCF path | positional VCF also accepted |
| `--parents`, `-P` | Parent sample names: high/resistant,low/susceptible | none |
| `--bulks`, `-B` | Bulk sample names: high,low | none |
| `--parents-depth`, `-p` | Minimum parent depth | `5,5` |
| `--bulks-depth`, `-b` | Minimum bulk depth | `40,40` |
| `--bulk-sizes`, `-S` | Individuals per bulk | `20,20` |
| `--window-size`, `-w` | Smoothing window size in bp | `2000000` |
| `--step-size`, `-s` | Window step size in bp | `100000` |
| `--population`, `-m` | `F2`, `F3`, `BC`, or `RIL` | `F2` |
| `--rep` | Number of simulations | `1000` |
| `--alpha` | Significance levels | `0.05,0.01` |
| `--min-qtl-length` | Minimum QTL interval length | `100000` |
| `--merge-distance` | Merge nearby QTL windows within this distance | `500000` |
| `--out`, `-o` | Output directory | current directory |

## Input Requirements

The VCF should contain:

- a standard `#CHROM` header with sample columns
- `GT` genotype values
- `AD` allele depths
- `DP` read depth

The app uses `vcfgo` for parsing and handles `.vcf` and `.vcf.gz`.

Markers are retained when:

- there is exactly one real ALT allele
- hard filters pass
- required samples are present and deep enough
- parent genotypes are homozygous and different, when two parents are supplied
- samples contain only REF or the target ALT allele

## Outputs

The output directory contains:

| File | Contents |
| --- | --- |
| `markers.tsv` | Marker-level high/low allele counts and statistics |
| `smoothed.tsv` | Smoothed statistics, robust z-scores, and adaptive thresholds |
| `brm_blocks.tsv` | BRM-style significant blocks |
| `qtls.tsv` | QTL intervals from permutation, z-score, consensus, MaxZ, and BRM high-confidence methods |
| `plots/GoBSAseq_best_of_both.html` | Interactive go-echarts report |

## Statistics

Marker-level statistics include:

- `high_si`
- `low_si`
- `delta_si`
- `g_statistic`
- `ed4`
- `lod`
- `beta_binomial_bf`
- `brm`

`ed4` is used because it amplifies strong allele-frequency separation and suppresses low-amplitude noise better than a plain Euclidean distance in this setting.

## Pipeline Stages

1. Open VCF and resolve samples.
2. Infer analysis mode.
3. Apply hard filters and BSA-specific informative-marker filters.
4. Orient allele depths to the high/resistant parent allele when parents are available.
5. Calculate marker statistics.
6. Assign BRM seed values.
7. Smooth statistics using tricube spatial weights and depth weights.
8. Simulate depth-adaptive thresholds for each window.
9. Apply robust z-score normalisation.
10. Detect permutation, z-score, consensus, MaxZ, BRM, and high-confidence QTLs.
11. Write TSV outputs.
12. Render go-echarts plots.

## Core Function Map

### CLI

- `cmd/root.go`: Cobra command setup, flag parsing, positional VCF support, interactive sample selection.
- `promptSamplesFromVCF`: Reads sample names from the VCF header and prompts for parents/bulks.
- `configFromFlags`: Converts CLI flags and prompt selections into analysis and filter config structs.

### Pipeline

- `run.Run`: Main orchestrator. Handles logging, output directories, parsing, statistics, smoothing, thresholds, QTL detection, and plots.
- `inferMode`: Chooses the appropriate analysis path from selected parent and bulk samples.
- `readMarkers`: Opens the VCF with `vcfgo`, validates sample names, streams variants, and keeps informative markers.

### Filtering And Allele Orientation

- `singleRealAlt`: Requires exactly one real ALT allele.
- `passesHardFilter`: Applies SNP/INDEL hard filters from `HardFilterConfig`.
- `sampleHasOnlyRefOrAlt`: Ensures samples contain only REF or the selected ALT allele.
- `orientedDepth`: Converts REF/ALT depths into tracked/high-parent allele depth and other allele depth.

### Statistics

- `calcTwoBulkStats`: Calculates two-bulk statistics.
- `calcOneBulkStats`: Calculates one-bulk statistics against expected 0.5 allele frequency.
- `gStatistic`: Likelihood-ratio G-statistic.
- `lod`: Two-bulk LOD score.
- `betaBinomialBF`: Beta-binomial Bayes factor.

### Smoothing And Thresholds

- `smoothMarkers`: Creates windows and smooths all statistics.
- `tricubeWeight`: Spatial kernel inherited from the older app's stronger smoother.
- `weightedMean`: Combines tricube spatial weights with depth weights.
- `localLinear`: Used for BRM-style local regression.
- `attachAdaptiveThresholds`: Adds per-window, depth-adaptive simulated thresholds.
- `calcThresholdsCached`: Caches thresholds by depth and expected allele frequency.

### QTL Detection

- `attachRobustZ`: Computes robust z-scores with trimmed median/MAD background.
- `detectPermutationQTLs`: Calls QTLs from adaptive threshold crossings.
- `detectZQTLs`: Calls z=2 and z=3 robust-z QTLs.
- `detectConsensusQTLs`: Calls intervals where multiple statistics fire together.
- `detectMaxZQTLs`: Calls intervals from the maximum absolute z-score across statistics.
- `detectBRMBlocks`: Finds BRM-style significant blocks.
- `intersectWithBRM`: Marks QTLs overlapping BRM blocks as high-confidence.

## Notes

- The app is intentionally a fresh directory and does not modify the older app or the optimized app.
- Gene-space flags are accepted and logged. The QTL output is ready for downstream annotation; direct SnpEff/GeneSpace execution can be enabled after validating local database paths.
- For large VCFs, increase `--rep` for better thresholds and tune `--window-size` / `--step-size` to match marker density.

