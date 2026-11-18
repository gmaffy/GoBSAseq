# GoBSAseq

**A high-performance Bulk Segregant Analysis (BSAseq) pipeline in Go**

GoBSAseq identifies genomic regions associated with phenotypic traits by applying hard filtering to VCF files, computing allele frequency statistics, performing Gaussian smoothing, and detecting QTL intervals using Z-score thresholding, Monte Carlo simulations, and Bayesian regression modeling (BRM).

## Quick Start

```bash
# Build
git clone https://github.com/gmaffy/GoBSAseq.git
cd GoBSAseq
go build -o gobsaseq

# Run (two-parent, two-bulk F2 example)
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

# Check output
ls results/goBSAseqResults/<timestamp>/
# plots/  qtls/  stats/
```

**Requirements**: Go 1.26.3+ ([download](https://go.dev/dl/)). Python 3 optional (for snpEff gene annotation).

---

## Usage

```bash
gobsaseq -V <vcf> -P <parents> -B <bulks> [options]

# Or with BAMs/reads instead of (or alongside) a VCF:
gobsaseq -V <vcf> -P <parents> -B <bulks> \
  --parents-bams hp.bam,lp.bam --bulks-bams hb.bam,lb.bam [options]

gobsaseq -V <vcf> -P <parents> -B <bulks> \
  --hp-reads fwd.fq,rev.fq --lp-reads fwd.fq,rev.fq \
  --hb-reads fwd.fq,rev.fq --lb-reads fwd.fq,rev.fq [options]
```

### CLI Flags

#### Input & Samples

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-V` / `--variant` | string | — | Path to VCF file (`.vcf` or `.vcf.gz`) |
| `-P` / `--parents` | string | — | Parent sample names: `high,low` (use `None` to skip) |
| `-B` / `--bulks` | string | — | Bulk sample names: `high,low` (use `None` to skip) |
| `-p` / `--parents-depth` | string | `5,5` | Min read depth for parents (high,low) |
| `-b` / `--bulks-depth` | string | `40,40` | Min read depth for bulks (high,low) |
| `-S` / `--bulk-sizes` | string | `20,20` | Number of individuals in bulks (high,low) - **Used by the null-genome simulation that sets thresholds; set it accurately** |
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
| `-w` / `--window-size` | int64 | `2000000` | Gaussian kernel width (bp); σ = width / 2 |
| `-s` / `--step-size` | int64 | `100000` | Reserved for future regular-grid summaries; marker-centred smoothing does not currently use it |
| `-m` / `--population` | string | `F2` | Population: `F2`, `F3`, `RIL`, `BC{n}H`, `BC{n}L`, `BC{n}F2H`, `BC{n}F2L` |
| `--brm-alpha` | float64 | `0.01` | BRM genome-wide significance level: chance that a genome with no QTL gets any BRM block |
| `--rep` | int | `1000` | Number of simulated null genomes for thresholds |
| `--seed` | int | `1` | Random seed for the null simulation (same seed + same input = identical thresholds) |

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

Supply all of `--snpEffDB`, `--gene-descriptions`, `--prg` and `--gff` to run it, or none to skip it. The run never prompts; a partial set is an error.

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

Use `None` for missing roles: `-P None -B bulk_high,bulk_low` for bulks-only analysis.

---

## Output

Results are written to `<out>/goBSAseqResults/<timestamp>/`:

```
stats/
  *.raw.tsv                        # Per-variant raw statistics
  *.smoothed_and_normalised.tsv    # Gaussian smoothed + robust Z-scores
  *.brm_blocks.tsv                 # BRM-detected blocks
  *.hardfiltered.vcf.gz (+.tbi)    # Filtered VCF with tabix index
  *.lowqual.vcf.gz                 # Variants failing hard filter
qtls/
  *.compositez.qtl.xlsx            # CompositeZ QTL intervals (empirical threshold)
  *.final.qtl.xlsx                 # Interval-level union of CompositeZ + BRM evidence
plots/
  *.individual_plots.html          # Per-statistic charts with thresholds
  *.robust_z_overlay.html          # All Z-scores overlaid per chromosome
  *.composite_signal.html          # CompositeZ + MaxAbsZ overview
genespace/  (if annotation enabled)
  *.genes_in_qtls.tsv              # Genes overlapping QTL regions
```

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
4. Thresholds          → **Null-genome simulation**: genome-wide thresholds from simulated bulks with recombination
5. BRM Detection       → Bayesian regression model block detection (incorporates bulk size)
6. QTL Detection       → Regions exceeding the genome-wide CompositeZ P99 threshold
7. Individual QTLs     → Per-statistic intervals exceeding genome-wide P99 thresholds
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
| **BRM** | Block where \|ΔSI\| ≥ `u·σ`, σ = `√(p(1−p)·(1/(V·n₁) + 1/(V·n₂)))` (pool sampling only; one bulk: `√(p₀(1−p₀)/(V·n))`); `u` = (1−α) quantile of the per-genome max of \|ΔSI\|/σ in the null simulation, so α is genome-wide |
| **Thresholds** | Quantiles of the genome-wide maximum of each smoothed track across simulated null genomes (see below) |

### Population Structures

| Structure | Expected p₀ | Description |
|-----------|-------------|-------------|
| F2 / F3 / RIL | 0.5 | Equal segregation |
| BC1H | 0.75 | Backcross to high parent |
| BC1L | 0.25 | Backcross to low parent |
| BC2H | 0.875 | Two backcrosses to high parent |
| BC2L | 0.125 | Two backcrosses to low parent |

---

## Significance Thresholds (Null-Genome Simulation)

Smoothed statistics are driven, under the null, mostly by **genetic drift**: each bulk holds a limited number of individuals, so long parental segments are over- or under-represented by chance and linked markers deviate together. Thresholds therefore come from simulating whole null genomes, not independent SNPs.

For each of `--rep` null genomes:

1. Every bulk individual is simulated through the cross (`-m`: F2/F3/F4/RIL selfing or BCn backcrossing) chromosome by chromosome, with Haldane crossovers. Each chromosome is assumed to be **100 cM long** (about one crossover per meiosis), a typical value for plants.
2. At every real marker the pool allele frequency is the simulated dosage, and reads are drawn as `Binomial(observed depth at that marker, pool frequency)`, so the real marker density and depth profile are used.
3. The pipeline's own statistics, depth-weighted Gaussian smoothing, robust-Z normalisation and CompositeZ formula are applied.
4. The genome-wide maximum (and minimum, for signed tracks) of each track is recorded.

**P99** is the 99th percentile of those per-genome maxima: only 1% of genomes without any QTL produce a call anywhere, so it controls the **genome-wide** false-positive rate. Signed tracks (SI, ΔSI, AFDev) split it across both tails. QTLs are called at P99 on every chromosome; P95 values are reported in the thresholds TSV. The per-genome extremes are written to `*.null_genome_extremes.tsv`.

Accurate `--bulk-sizes` matter: smaller bulks mean more drift and higher thresholds. Weak QTLs in small bulks can be genuinely indistinguishable from drift at a 1% genome-wide error rate.

### QTL Detection Methods

GoBSAseq provides **two complementary QTL detection approaches**:

| Method | Threshold | Primary/Secondary | Output File |
|--------|-----------|----------------|-------------|
| **Composite Z-score** | Genome-wide null-simulation threshold | **Primary** | `*.compositez.qtl.xlsx` |
| **Individual statistics** | Genome-wide null-simulation thresholds | **Secondary** | `*.individual.qtl.xlsx` |
| **BRM Blocks** | Genome-wide null-simulation multiplier × drift σ | **Secondary** | `*.brm_blocks.tsv` |

**CompositeZ and individual-statistic detection are run by default.** The threshold TSV records the exact threshold used for every retained marker.

### Localising a QTL: Support Intervals

The threshold decides **whether** a region holds a QTL. It does not say **where**: with a strong QTL in an early generation (F2), linkage spans tens of Mb and the region above the threshold can be the whole chromosome. Each QTL is therefore also given a **support interval**: from the highest smoothed LOD in the region, it extends in both directions while LOD stays within **1.5** of the peak (a 1.5-LOD drop).

- QTL tables list the support interval first (`Support_Start`, `Support_End`, `Support_PeakPos`, `Support_PeakLOD`), then the above-threshold region (`Region_Start`, `Region_End`).
- `*.brm_blocks.tsv` has `SUPPORT_*` columns appended.
- Plots shade support intervals red; BRM blocks stay amber.
- Gene-space analysis searches the support interval.

On `Data/sim_F2_strongQTL.vcf.gz` (true QTL at chr2:15.0 Mb) the above-threshold region is 0–30 Mb and the support interval 14.60–16.29 Mb. The 1.5 drop has not yet been calibrated for smoothed BSA-seq LOD.

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
# Specify accurate bulk sizes: they drive the null-genome simulation
./gobsaseq \
  -V deep_seq.vcf.gz \
  -P parent_a,parent_b \
  -B bulk_high,bulk_low \
  -p 20,20 -b 200,200 \
  -S 150,150  # <-- Bulk sizes: 150 individuals each
  -w 2000000 -s 100000 \
  -m F2 --rep 2000  # More null genomes = more precise P99
  -o results_deep

# Thresholds account for drift in bulks of this size
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
  thresholds.go      Genome-wide thresholds from null-genome extremes
  nullsim.go         Null-genome simulation (cross, recombination, reads, smoothing)
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
