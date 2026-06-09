# GoBSAseq

**GoBSAseq** is a fast, single-binary pipeline for **Bulk Segregant Analysis (BSA-seq)**, written in Go. It takes a multi-sample VCF produced from sequencing a mapping population's parents and phenotypic bulks, and runs the full analysis — variant filtering, per-variant association statistics, signal smoothing, significance thresholds, QTL/region detection and gene-space annotation — emitting result tables and interactive HTML plots.

BSA-seq maps the genomic regions (QTLs) underlying a trait by contrasting allele frequencies between two pools (bulks) of individuals selected from opposite tails of a phenotype distribution. Regions linked to the trait show a skew in allele frequency between the bulks; GoBSAseq quantifies and localizes that skew.

## Features

- **Single static binary** built on [Cobra](https://github.com/spf13/cobra) — no Python/R runtime required.
- **Flexible experimental designs.** Automatically dispatches based on which samples you supply: two parents + two bulks, two parents + one bulk, one parent + two bulks, single parent + single bulk, or bulks-only (11 combinations in total).
- **Multiple BSA statistics** computed per variant: SNP-index (per bulk), ΔSNP-index, Magwene **G-statistic**, **Euclidean distance (ED)**, **LOD**, and a **beta-binomial log Bayes factor**; one-bulk designs add an allele-frequency-deviation statistic against the expected frequency.
- **Depth-weighted Gaussian smoothing** per chromosome to suppress coverage-driven noise without the edge artefacts of hard sliding windows.
- **Robust normalisation and consolidation** — per-statistic robust Z-scores (median/MAD) combined into a single composite signal via Stouffer's Z.
- **Two complementary region callers** — simulation-based significance thresholds and the **Block Regression Method (BRM)** for QTL interval detection, with results merged into final candidate intervals.
- **Population-aware** expected allele frequencies for F2, F3, RIL, and backcross populations (BC1H, BC1L, BC2H, BC2L, …).
- **Interactive plots** (`go-echarts`): smoothed statistic tracks with p95/p99 reference lines and shaded BRM blocks, exported as standalone HTML.
- **Optional gene-space analysis** — annotate candidate intervals using snpEff, a GFF3/CDS/protein reference, gene descriptions, and BLAST hits.

## How it works

The pipeline (`run.Run`) executes these stages on the input VCF:

1. **Hard filtering** — GATK-style hard filters applied separately to SNPs and INDELs (`QD`, `QUAL`, `SOR`, `FS`, `MQ`, `MQRankSum`, `ReadPosRankSum`). Writes a hard-filtered, bgzipped VCF.
2. **Raw statistics** — per-variant allele depths are turned into SNP-indices and the association statistics listed above (`stats.RawStats`).
3. **Smoothing** — depth-weighted Gaussian kernel per chromosome with bandwidth `σ = WindowSize / 2` and a 3σ truncation; sparse windows (fewer than the minimum contributing variants) are skipped (`stats.SmoothAndNormalise`).
4. **Normalisation & consolidation** — robust Z-score `Z = (x − median) / (1.4826 × MAD)` per statistic, combined into a composite `Z = Σ Zᵢ / √k`.
5. **Thresholds** — simulation/permutation-based significance thresholds at level `--brm-alpha` over `--rep` iterations (`stats.CalculateThresholds`).
6. **BRM blocks** — contiguous windows exceeding the block threshold are flagged as candidate QTL blocks (`stats.RunBRM`).
7. **Plots** — interactive HTML charts of the smoothed/composite signals with thresholds and BRM overlays (`plots.GeneratePlots`).
8. **QTL detection & merging** — peaks are called and merged with BRM blocks into final intervals (`stats.DetectQTLs`, `stats.MergeQTLsAndBRM`).
9. **Gene-space annotation** *(optional)* — annotates merged intervals when snpEff/GFF/gene-description/BLAST inputs are provided (`stats.GeneSpaceFromMerged`).

Results are written to a timestamped directory under `<out>/goBSAseqResults/`.

## Installation

Requires [Go](https://go.dev/dl/) (see the version in [`go.mod`](go.mod)).

```bash
git clone https://github.com/gmaffy/GoBSAseq.git
cd GoBSAseq
go build -o gobsaseq .      # build a local binary
# or install onto your PATH:
go install github.com/gmaffy/GoBSAseq@latest
```

## Usage

The input is a single multi-sample VCF (`.vcf` or `.vcf.gz`) containing the parents and bulks. If sample names are not provided on the command line, GoBSAseq lists the samples in the VCF and prompts you to select them interactively.

```bash
gobsaseq \
  --variant joint_called.vcf.gz \
  --parents HighParent,LowParent \
  --bulks   HighBulk,LowBulk \
  --bulk-sizes 30,30 \
  --population F2 \
  --window-size 2000000 \
  --out results/
```

To run a design with fewer samples, pass `None` (or omit a name) for the missing role — e.g. `--parents HighParent,None --bulks HighBulk,LowBulk` runs the "high parent + two bulks" design.

### Key flags

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--variant` | `-V` | | Input VCF file (`.vcf` / `.vcf.gz`) |
| `--parents` | `-P` | | High,Low parent sample names (comma-separated) |
| `--bulks` | `-B` | | High,Low bulk sample names (comma-separated) |
| `--bulk-sizes` | `-S` | `20,20` | Number of individuals in each bulk |
| `--parents-depth` | `-p` | `5,5` | Minimum depth per parent |
| `--bulks-depth` | `-b` | `40,40` | Minimum depth per bulk |
| `--population` | `-m` | `F2` | Population type (`F2`, `F3`, `RIL`, `BC1H`, `BC1L`, `BC2H`, `BC2L`, …) |
| `--window-size` | `-w` | `2000000` | Smoothing window size (bp); sets the Gaussian bandwidth |
| `--step-size` | `-s` | `100000` | Step size (bp) |
| `--brm-alpha` | | `0.05` | Significance level for thresholds / BRM |
| `--rep` | | `1000` | Number of simulations for threshold estimation |
| `--min-qtl-length` | | | Minimum QTL interval length (bp) |
| `--merge-distance` | | | Maximum gap (bp) when merging nearby intervals |
| `--out` | `-o` | `.` | Output directory |

SNP and INDEL hard-filter cutoffs are exposed as flags too (e.g. `--min-QD-SNP`, `--min-QUAL-SNP`, `--min-SOR-SNP`, `--min-FS-SNP`, `--min-MQ-SNP`, `--max-FS-INDEL`, `--max-SOR-INDEL`, …). Run `gobsaseq --help` for the full list and current defaults.

Gene-space annotation is enabled by additionally passing `--snpEffDB`, `--gff`, `--cds`, `--reference`, `--protein`, `--gene-descriptions` and `--prg`. If these are omitted you are prompted to continue without annotation.

## Output

A run creates `<out>/goBSAseqResults/<DD_MM_YYYY_HH_MM_SS>/` containing:

- `stats/` — the hard-filtered VCF (`GoBSAseq.<type>.hardfiltered.vcf.gz`), raw and smoothed/normalised per-variant statistic tables (TSV), threshold tables, and detected BRM blocks, QTLs and merged intervals.
- `plots/` — interactive HTML charts of the smoothed and composite signals with significance thresholds and BRM block overlays.
- Gene-space annotation tables for the candidate intervals (when annotation inputs are supplied).

## References

The statistical methods draw on established BSA-seq literature, including QTL-seq / ΔSNP-index, the Magwene *et al.* (2011) G-statistic, Euclidean-distance mapping, and block-regression-style interval detection. See [`docs.txt`](docs.txt) and [`brm_and_plot.md`](brm_and_plot.md) for implementation notes.

## License

See [`LICENSE`](LICENSE).
