# GoBSAseq

**A high-performance Bulk Segregant Analysis (BSAseq) pipeline in Go**

GoBSAseq identifies genomic regions associated with phenotypic traits by applying hard filtering to VCF files, computing allele frequency statistics, performing Gaussian smoothing, and detecting QTL intervals using Z-score thresholding, Monte Carlo simulations, and Bayesian regression modeling (BRM).

## Getting Started

### Prerequisites

- **Go 1.26+** — [download](https://go.dev/dl/)
- A **multi-sample VCF** (`.vcf` or `.vcf.gz`) with `GT`, `AD`, `DP`, and GATK-style INFO fields (`QD`, `FS`, `SOR`, `MQ`, …), **or** bulk BAM/CRAM files plus a reference FASTA for variant calling

### Install & build

```bash
git clone https://github.com/gmaffy/GoBSAseq.git
cd GoBSAseq
go build -o gobsaseq .
```

### Run (VCF input)

The simplest path is an existing multi-sample VCF. Assign samples with `-P` (parents) and `-B` (bulks); use `None` for a missing role (e.g. bulks-only: `-P None,None`).

```bash
./gobsaseq \
  -V my_variants.vcf.gz \
  -P parent_high,parent_low \
  -B bulk_high,bulk_low \
  -p 10,10 \
  -b 100,100 \
  -S 250,250 \
  -w 2000000 \
  -m F2 \
  -o results
```

If `-P` / `-B` names are omitted or not found in the VCF header, GoBSAseq prompts you interactively to pick samples (enter `0` for **None**).

Results are written to `<out>/goBSAseqResults/`:

```
stats/   # filtered VCF, raw & smoothed TSVs, BRM blocks
qtls/    # individual, composite-Z, and merged QTL Excel files
plots/   # interactive HTML genome plots
```

---

## Usage

```bash
gobsaseq [options]
```

### Input modes

Provide **one** input type only (they cannot be combined):

| Mode | Flags | Notes |
|------|-------|-------|
| **VCF** | `-V file.vcf.gz` | Most common. No BAM or read flags. |
| **BAM** | `--parents-bams hp.bam,lp.bam` and/or `--bulks-bams hb.bam,lb.bam` | Variant calling runs first; requires `-r` (reference) and `-o`. Parent BAMs are optional. |
| **Reads** | `--hp-reads`, `--lp-reads`, `--hb-reads`, `--lb-reads` | Not yet implemented. |

### Sample roles

`-P` and `-B` take comma-separated names: `high,low`. Use `None` to skip a role.

| Example | Mode | Description |
|---------|------|-------------|
| `-P a,b -B x,y` | `2p2b` | Two parents, two bulks |
| `-P a,None -B x,y` | `hp2b` | High parent + both bulks |
| `-P None,None -B x,y` | `2b` | Bulks only |

Analysis mode is inferred automatically from which roles are set.

### Essential flags

| Flag | Default | Description |
|------|---------|-------------|
| `-V` | — | Input VCF (VCF mode) |
| `-P` | `,` | Parent sample names (`high,low`) |
| `-B` | `,` | Bulk sample names (`high,low`) |
| `-p` | `5,5` | Min depth in parents |
| `-b` | `40,40` | Min depth in bulks |
| `-S` | `20,20` | Individuals per bulk (for Monte Carlo thresholds) |
| `-w` | `2000000` | Gaussian smoothing window (bp); σ = window / 2 |
| `-m` | `F2` | Population: `F2`, `F3`, `RIL`, `BC1H`, `BC1L`, `BC2H`, `BC2L` |
| `-o` | `.` | Output directory |
| `-r` | — | Reference FASTA (required for BAM mode) |
| `--rep` | `1000` | Monte Carlo simulations for thresholds |
| `--light-filtering` | off | Skip GATK hard-filter thresholds |

### BAM example

```bash
./gobsaseq \
  --parents-bams parent_a.bam,parent_b.bam \
  --bulks-bams bulk_high.bam,bulk_low.bam \
  -r reference.fa \
  -p 10,10 -b 100,100 -S 250,250 \
  -w 2000000 -m F2 \
  -o results
```

### Gene space (optional)

Add annotation files to map QTLs to genes. If omitted, the pipeline skips gene-space analysis (or prompts to continue without it):

```bash
  --snpEffDB GRCh38.99 --gff genes.gff3 \
  --gene-descriptions gene_desc.tsv --prg pangenome.prg
```

See [CLI Flags](#cli-flags) below for the full option list.

---

#### Input & Samples

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-V` / `--variant` | string | — | Path to VCF file (`.vcf` or `.vcf.gz`) |
| `-P` / `--parents` | string | — | Parent sample names: `high,low` (use `None` to skip) |
| `-B` / `--bulks` | string | — | Bulk sample names: `high,low` (use `None` to skip) |
| `-p` / `--parents-depth` | string | `5,5` | Min read depth for parents (high,low) |
| `-b` / `--bulks-depth` | string | `40,40` | Min read depth for bulks (high,low) |
| `-S` / `--bulk-sizes` | string | `20,20` | Number of individuals in bulks (high,low) - **Important for deep sequencing: used in Monte Carlo null model** |
| `-o` / `--out` | string | `.` | Output directory |

#### BAM & Read Inputs (Optional)

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--parents-bams` | string | — | Parent BAM files: `high.bam,low.bam` |
| `--bulks-bams` | string | — | Bulk BAM files: `high.bam,low.bam` |
| `--hp-reads` | string | — | High parent reads: `fwd.fq,rev.fq` |
| `--lp-reads` | string | — | Low parent reads: `fwd.fq,rev.fq` |
| `--hb-reads` | string | — | High bulk reads: `fwd.fq,rev.fq` |
| `--lb-reads` | string | — | Low bulk reads: `fwd.fq,rev.fq` |

#### Smoothing & Statistics

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-w` / `--window-size` | int64 | `2000000` | Gaussian smoothing window (bp); σ = window / 2 |
| `-s` / `--step-size` | int64 | `100000` | BRM sliding-window step size (bp) |
| `-m` / `--population` | string | `F2` | Population: `F2`, `F3`, `RIL`, `BC1H`, `BC1L`, `BC2H`, `BC2L` |
| `--brm-alpha` | float64 | `0.05` | BRM significance level |
| `--rep` | int | `1000` | Number of simulations for thresholds |

#### Variant Filtering

**SNP** (defaults follow GATK best practices):

| Flag | Default | Description |
|------|---------|-------------|
| `--min-QD-SNP` | `2.0` | Min quality-by-depth |
| `--min-QUAL-SNP` | `30.0` | Min variant quality |
| `--min-MQ-SNP` | `40.0` | Min RMS mapping quality |
| `--min-SOR-SNP` | `3.0` | Max strand odds ratio |
| `--min-FS-SNP` | `60.0` | Max FisherStrand bias |
| `--min-MQRankSum-SNP` | `-12.5` | Min MQ rank sum |
| `--min-ReadPosRankSum-SNP` | `-8.0` | Min read position rank sum |

**INDEL** (more permissive than SNP defaults):

| Flag | Default | Description |
|------|---------|-------------|
| `--min-QD-INDEL` | `2.0` | Min quality-by-depth |
| `--min-QUAL-INDEL` | `30.0` | Min variant quality |
| `--max-FS-INDEL` | `200.0` | Max FisherStrand bias |
| `--max-SOR-INDEL` | `10.0` | Max strand odds ratio |
| `--min-ReadPosRankSum-INDEL` | `-20.0` | Min read position rank sum |

#### Gene Space Analysis (Optional)

| Flag | Description |
|------|-------------|
| `--snpEffDB` | snpEff database ID (e.g., `GRCh38.p13`) |
| `--gff` | GFF3 annotation file |
| `--cds` | CDS FASTA file |
| `-r` / `--reference` | Reference genome FASTA |
| `--protein` | Protein FASTA file |
| `--gene-descriptions` | Gene descriptions TSV (`gene_id\tdescription`) |
| `--prg` | PRG (pangenome graph) blast file |

---

## Analysis Modes

The mode is auto-detected based on which samples are provided (`None` = not provided):

| Mode | Parents | Bulks | Description |
|------|---------|-------|-------------|
| `2p2b` | both | both | Classic two-parent two-bulk (F2, RIL, backcross) |
| `2phb` | both | high only | Two parents, high bulk only |
| `2plb` | both | low only | Two parents, low bulk only |
| `hp2b` | high only | both | High parent with both bulks |
| `lp2b` | low only | both | Low parent with both bulks |
| `hphb` | high only | high only | Single parent, single bulk |
| `hplb` | high only | low only | Single parent, single bulk |
| `lphb` | low only | high only | Single parent, single bulk |
| `lplb` | low only | low only | Single parent, single bulk |
| `2b` | none | both | Bulks only, no parents |

**Two-bulk modes** compute: ΔSI, G-statistic, ED⁴, LOD, Bayes factor.
**One-bulk modes** compute: AF deviation, one-bulk G-stat, LOD, Bayes factor.

---

## Input: VCF Requirements

GoBSAseq expects a **multi-sample VCF** (v4.1/4.2) with:

- **GT** (genotype): `0/0`, `0/1`, `1/1`, or `./.`
- **AD** (allele depth): `ref_depth,alt_depth`
- **DP** (read depth)
- **INFO fields**: `QD`, `FS`, `SOR`, `MQ`, `MQRankSum`, `ReadPosRankSum` (for hard filtering)

```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  parent1  parent2  bulk_high  bulk_low
chr1    1000 .   A    T    60    PASS    QD=20.0;FS=0.0;SOR=0.8;MQ=60  GT:AD:DP  0/0:50,0:50  1/1:0,45:45  0/1:25,25:50  0/1:30,14:44
```

Use `None` for missing roles: `-P None,None -B bulk_high,bulk_low` for bulks-only analysis.

---

## Output

Results are written to `<out>/goBSAseqResults/`:

```
stats/
  *.raw.tsv                        # Per-variant raw statistics
  *.smoothed_and_normalised.tsv    # Gaussian smoothed + robust Z-scores
  *.brm_blocks.tsv                 # BRM-detected blocks
  *.hardfiltered.vcf.gz (+.tbi)    # Filtered VCF with tabix index
  *.lowqual.vcf.gz                 # Variants failing hard filter
qtls/
  *.individual.qtl.xlsx            # Per-statistic QTL peaks
  *.compositez.qtl.xlsx            # CompositeZ QTL peaks
  *.final.qtl.xlsx                 # Merged CompositeZ + BRM intervals
plots/
  *.individual_plots.html          # Per-statistic charts with thresholds
  *.robust_z_overlay.html          # All Z-scores overlaid per chromosome
  *.composite_signal.html          # CompositeZ + MaxAbsZ overview
genespace/  (if annotation enabled)
  *.genes_in_qtls.tsv              # Genes overlapping QTL regions
```

Primary QTL output is `qtls/*.final.qtl.xlsx` (merged CompositeZ and BRM intervals).

### Key Columns

**Raw TSV**: `CHROM`, `POS`, `REF`, `ALT`, genotypes, `HighSI`/`LowSI` (selection index), `DeltaSI`, `Gstat`, `ED4`, `LOD`, `BBLogBF`, `Depth`. One-bulk modes add `P0`, `AFDev`, `Gstat1`, `LOD1`, `BBLogBF1`.

**Smoothed TSV**: Adds `Sm_*` (smoothed values), `Z_*` (robust Z-scores), `CompositeZ` (Stouffer), `MaxAbsZ`.

**Merged QTLs**: `CHROM`, `START`, `STOP`, `SOURCE` (`ZScore`/`BRM`/`ZScore+BRM`), `Z_PEAK`, `BRM_PEAK`, `BRM_THRESHOLD`.

---

## Pipeline Overview

```
1. Hard Filter         → Remove low-quality variants (SNP/INDEL-specific thresholds)
2. Raw Statistics      → Per-variant SI, ΔSI, G-stat, LOD, Bayes factor
3. Smooth & Normalise  → Gaussian kernel smoothing, robust Z-scores (MAD), Stouffer composite Z
4. Thresholds          → **Two-stage Monte Carlo simulation** for empirical thresholds per depth + bulk size
5. BRM Detection       → Bayesian regression model block detection (incorporates bulk size)
6. QTL Detection       → Regions with |CompositeZ| ≥ 3.0 (**primary method**)
7. MC QTL Detection    → **QTLs using per-variant MC thresholds** (fully sound for deep sequencing)
8. Merge               → Union of Z-score and BRM intervals
9. Plots               → Interactive HTML charts (echarts) with all thresholds
10. Gene Space         → Annotate QTLs with genes (optional, requires --gff)
```

### Statistical Methods

| Method | Formula / Description |
|--------|----------------------|
| **Selection Index (SI)** | `freq(high allele) / freq(all alleles)` per bulk |
| **ΔSI** | `SI_high − SI_low`; large values suggest QTL |
| **G-statistic** | `2 Σ nᵢ log(oᵢ/eᵢ)`; likelihood ratio test (~χ²) |
| **LOD** | `log₁₀(L_alt / L_null)`; LOD > 3 = significant |
| **Bayes Factor** | Beta-binomial log₁₀ BF with Beta(0.5, 0.5) prior |
| **Robust Z-score** | `(x − median) / (1.4826 × MAD)`; outlier-resistant |
| **Composite Z** | Stouffer: `Σ Zᵢ / √k`; combines all Z-scores |
| **BRM** | Threshold = `u_α √((n₁+n₂)·p(1−p) / (V_scale·n₁n₂))` for two-bulk |
| **MC Thresholds** | Two-stage: `Binomial(bulk_size, p₀) → Binomial(depth, realized_af)`; empirical per-depth thresholds |

### Population Structures

| Structure | Expected p₀ | Description |
|-----------|-------------|-------------|
| F2 / F3 / RIL | 0.5 | Equal segregation |
| BC1H | 0.75 | Backcross to high parent |
| BC1L | 0.25 | Backcross to low parent |
| BC2H | 0.875 | Two backcrosses to high parent |
| BC2L | 0.125 | Two backcrosses to low parent |

---

## Deep Sequencing & Monte Carlo Thresholds

### Two-Stage Null Model (For Accurate Deep Sequencing)

For deeply sequenced datasets, GoBSAseq uses a **two-stage Monte Carlo null model** that properly accounts for the finite population of individuals in each bulk:

**Stage 1**: Sample the realized allele count in the bulk population
```
Alt Alleles ~ Binomial(n = bulk_size, p = p₀)
```

**Stage 2**: Sample the observed reads from that realized frequency
```
Reads ~ Binomial(n = depth, p = realized_allele_frequency)
```

**Why This Matters**: In deep sequencing, many reads may come from the same individuals. The old model (directly sampling reads from p₀) overestimated variance, leading to excessive false positives. The two-stage model correctly models this structure.

**When to Use**: This model is automatically used when `--bulk-sizes` is provided. For best results with deep sequencing (>50x coverage per bulk), always specify accurate bulk sizes.

### QTL Detection Methods

GoBSAseq provides **two complementary QTL detection approaches**:

| Method | Threshold | Primary/Secondary | Output File |
|--------|-----------|----------------|-------------|
| **Composite Z-score** | |CompositeZ| ≥ 3.0 | **Primary** | `*.qtls.tsv` |
| **Monte Carlo** | Empirical per-variant thresholds | **Secondary** | `*.mc_qtls.tsv` |
| **BRM Blocks** | Analytical window-based | **Secondary** | `*.brm_blocks.tsv` |

**Both CompositeZ ≥ 3.0 and Monte Carlo threshold-based detection are run by default.** CompositeZ is the primary method, while Monte Carlo detection provides a fully sound alternative for deep sequencing that uses the simulated thresholds directly.

---

## Examples

### Two-Parent Two-Bulk (F2)

```bash
./gobsaseq \
  -V data/crosses.vcf.gz \
  -P parent_a,parent_b \
  -B phenotype_high,phenotype_low \
  -p 10,10 -b 100,100 -S 250,250 \
  -w 1000000 -s 50000 \
  -m F2 --brm-alpha 0.05 --rep 1000 \
  -o results
```

### Single-Parent Single-Bulk

```bash
./gobsaseq \
  -V data/variants.vcf.gz \
  -P parent_high,None \
  -B bulk_high,None \
  -p 10,10 -b 150,150 -S 500,500 \
  -w 2000000 -s 100000 -m F2 \
  -o results
# Mode: hphb — computes AF deviation, one-bulk G/LOD/BF
```

### Deep Sequencing (>50x Coverage)

```bash
# For deep sequencing, specify bulk sizes for accurate Monte Carlo thresholds
./gobsaseq \
  -V deep_seq.vcf.gz \
  -P parent_a,parent_b \
  -B bulk_high,bulk_low \
  -p 20,20 -b 200,200 \
  -S 150,150  # <-- Bulk sizes: 150 individuals each
  -w 2000000 -s 100000 \
  -m F2 --rep 2000  # More simulations for deep data
  -o results_deep

# The two-stage Monte Carlo model will automatically use bulk sizes
# to reduce false positives from deep coverage
```

---

## Troubleshooting

| Problem | Likely Cause | Solution |
|---------|-------------|----------|
| `bad sample string` | Malformed GT field in VCF | Validate with `bcftools check` |
| `not part of the VCF sample list` | Sample name mismatch | Check names: `bcftools query -l file.vcf` |
| 0 variants pass filtering | Thresholds too strict or non-standard INFO fields | Relax `--min-QUAL-*`, `--min-QD-*`; check INFO fields |
| Empty output TSVs | All variants removed by filters | Reduce stringency; verify samples are segregating |
| No QTLs detected | Low power or weak signal | Check max CompositeZ in smoothed TSV; try `--brm-alpha 0.10` |
| HTML plots won't render | Browser blocks local files | Serve with `python3 -m http.server 8000` |
| BRM and Z-score disagree | Expected — different methods | Consensus intervals (`ZScore+BRM`) are most robust |

---

## Development

### Project Structure

```
main.go              Entry point
cmd/root.go          CLI setup (Cobra)
run/run.go           Pipeline orchestration
filter/filter.go     VCF hard filtering + tabix indexing
stats/
  stats.go           Raw statistics
  smoothing.go       Gaussian smoothing + normalization
  detect.go          QTL detection + merging
  brm.go             BRM block detection
  thresholds.go      Empirical threshold calculation
  genespace.go       Gene annotation (optional)
plots/plots.go       HTML plotting (go-echarts)
utils/
  config.go          Configuration structs
  utils.go           Utility functions
```

### Build & Test

```bash
go build -o gobsaseq
go test ./...
go vet ./...
```

### Contributing

1. Fork → feature branch → add tests → `go test ./... && go vet ./...` → PR

**Priority areas**: unit tests for `stats/`, additional population structures, parallel per-chromosome processing, extended gene annotation formats.

---

## References

1. Michelmore et al. (1991). Identification of markers linked to disease-resistance genes by bulked segregant analysis. *PNAS*, 88(21), 9828–9832.
2. Magwene et al. (2011). Inferring the joint evolutionary history of multiple loci. *G3*, 1(5), 417–425.
3. Sokal & Rohlf (1995). *Biometry* (3rd ed.). Freeman. (G-statistic)
4. Morton (1955). Sequential tests for the detection of linkage. *AJHG*, 7(3), 277. (LOD)
5. Kass & Raftery (1995). Bayes factors. *JASA*, 90(430), 773–795.
6. Rousseeuw & Croux (1993). Alternatives to the median absolute deviation. *JASA*, 88(424), 1273–1283.
7. Stouffer et al. (1949). *The American Soldier* (Vol. 1). Princeton University Press.

---

## License

See [LICENSE](LICENSE).

**Issues & Feature Requests**: [GitHub Issues](https://github.com/gmaffy/GoBSAseq/issues)
