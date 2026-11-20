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
| `-w` / `--window-size` | int64 | `2000000` | Gaussian kernel width (bp); σ = width / 2 |
| `-s` / `--step-size` | int64 | `100000` | Reserved for future regular-grid summaries; marker-centred smoothing does not currently use it |
| `-m` / `--population` | string | `F2` | Population: `F2`, `F3`, `RIL`, `BC{n}H`, `BC{n}L`, `BC{n}F2H`, `BC{n}F2L` |
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

GATK semantics apply: a record fails only when an annotation is **present** and beyond its threshold (absent annotations pass, as in `VariantFiltration`). Records whose FILTER is not `PASS`/`.` are rejected. `--light-filtering` keeps FS/SOR but drops the two rank-sum tests, which are unreliable near allele fixation in pooled samples.

**BSA-seq filters** (applied per biallelic REF/ALT view, after multi-allelic splitting):

| Flag | Default | Description |
|------|---------|-------------|
| `--bulks-max-depth` / `--parents-max-depth` | off | Max depth caps (`hb,lb` / `hp,lp`) against collapsed repeats; uses the larger of FORMAT/DP and allele depth |
| `--min-GQ` | `0` | GQ floor for **parents only** (bulk GQ tracks allele frequency) |
| `--no-split-multiallelic` | off | Keep multi-allelic records whole; the first passing ALT represents the record |
| `--require-nonparent-allele` | off | One parent + one bulk modes: drop sites where the bulk has < 2 reads of the allele the parent does not carry |
| `--one-bulk-p0` | `empirical` | One-bulk modes: null frequency of the high-parent allele. `empirical` corrects the population's expected value for reference mapping bias, estimated from the data: the bias shifts the high-allele fraction up where the high allele is REF and down where it is ALT, so half the difference between those two groups of sites (interquartile means; bulk with ≥2 reads of each allele) is the bias. A QTL shifts both groups equally and cancels. Falls back to `theoretical` (the expected value) with < 200 sites per group. The estimate is printed in the log and the per-site value is in the `P0` column |
| `--cm-per-mb` | `0` | Recombination rate for the HMM LOD in cM/Mb; `0` assumes 1 Morgan (one crossover per meiosis) over each chromosome |
| `--null-reps` | `200` | Simulated no-QTL genomes for the genome-wide HMM LOD and CompositeZ thresholds |
| `--ref-parent` | `none` | `2b` only: which parent the reference genome is (`high`, `low` or `none`). `high`/`low` fixes allele polarity so ΔSI is signed. `none` makes `2b` use \|ΔSI\| (see below). Use `high`/`low` only if the reference is that parent or the same inbred line. Ignored when a parent sample is given |
| `--min-pooled-maf` | `-1` (auto) | Segregation floor for `2b`/`hp2b`/`lp2b`: pooled minor allele frequency of both bulks must reach this. Auto = 0.4 × the population's expected minor AF (0.2 for F2, 0.1 for BC1). `0` disables |

- **Bulk depth** (`-b`) is the REF + ALT allele depth from AD (the reads the SNP-index uses), not FORMAT/DP. Reads on other ALT alleles may be at most 5% of a sample's allele depth.
- **Parents** must be homozygous: ≥ 85% of allele depth on one allele (other-ALT reads count against both). Two-parent modes also require the parents to differ. Without AD, only a homozygous REF or target-ALT genotype counts.
- **One parent + one bulk** (`hphb`, `hplb`, `lphb`, `lplb`): a bulk fixed for the parent's allele cannot be told apart from a non-segregating site. Such sites are kept by default (they may be QTL markers); `--require-nonparent-allele` drops them.
- **Rejection report**: at the end of filtering, the number of records removed by each rule is printed, so the cost of every filter is visible.
- **Inputs**: `.vcf`, `.vcf.gz` (bgzip is decompressed on several threads; plain gzip on one thread) and FreeBayes-style RO/AO when AD is absent. For pooled bulks called with GATK HaplotypeCaller, set `--max-reads-per-alignment-start 0`; the default (50) caps AD in deep bulks. Site-level annotations (QD, FS, …) describe all samples in the VCF, so in a large cohort VCF they mostly reflect the other samples.
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

**Two-bulk modes** compute: ΔSI, G-statistic, ED⁴ and the pool-HMM LOD (see below).

**`2b` polarity.** Without a parent sample, SI is the REF-allele fraction. Unless the reference genome is one of the parents, REF comes from the high parent at some SNPs and the low parent at others, so the sign of ΔSI is arbitrary per SNP and would cancel when smoothed. With `--ref-parent none` (the default), `2b` therefore uses \|ΔSI\| (columns `AbsDeltaSI`, `Sm_AbsDeltaSI`, `Z_AbsDeltaSI`): its thresholds are upper-tail only, there is no `DeltaSI-` track, CompositeZ uses Z(\|ΔSI\|), and BRM uses the polarity-free form described below. QTLs can be located, but not which bulk carries which parent's allele. G and ED⁴ do not depend on polarity. The HMM LOD needs polarity and is not computed in this case; CompositeZ then uses \|ΔSI\| and G′ only. BRM uses a polarity-free form. It smooths the bias-corrected ΔSI², that is (SI_H − SI_L)² minus each bulk's read-sampling variance SI(1−SI)/(d−1), which estimates the squared pool frequency difference whatever allele is REF. BRM then tests its square root (\|ΔAF\|) against the usual threshold, with af at the population's expected frequency. Blocks have no direction. On simulated no-QTL genomes its false-block rate matched the signed BRM's (24 vs 22 of 100 genomes). If the reference is a parent, pass `--ref-parent high` or `low` to get signed ΔSI and BRM.
**One-bulk modes** compute: AF deviation, one-bulk G-stat and the pool-HMM LOD.

**HMM LOD** (after MULTIPOOL, Edwards & Gifford 2012). Each bulk is modelled as a pool whose hidden high-parent allele count changes along the chromosome only through crossovers, so all linked SNPs inform the same pool frequency, and pool sampling of the plants is part of the model. The pool prior and linkage depend on the population type (F2, F3, F4, RIL, BCn, BCnF2); recombination defaults to 1 Morgan per chromosome (`--cm-per-mb` sets a genome-wide rate); markers within 0.01 cM share one HMM step. Two-bulk LOD tests whether the bulks differ, with their shared frequency free, so segregation distortion that shifts both bulks is not scored as a QTL. One-bulk LOD tests the bulk against P0 (with the mapping-bias correction). It is computed in the raw-statistics stage and written to the raw TSV. The LOD is already linkage-aware, so it is not smoothed or normalised and is not part of CompositeZ: it stays in LOD units, as its own track with its own thresholds. Its thresholds are genome-wide: the 95%/99% quantiles of the maximum LOD over `--null-reps` simulated no-QTL genomes at the real markers and depths, run through the same scoring code. The same simulated genomes give the CompositeZ thresholds. The HMM LOD replaces the earlier per-SNP LOD and beta-binomial Bayes factor, which ignored pool sampling and linkage.

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
  *.lowqual.vcf.gz                 # Records where no ALT passed (kept samples only)
qtls/
  *.compositez.qtl.xlsx            # CompositeZ QTL intervals (empirical threshold)
  *.final.qtl.xlsx                 # Interval-level union of CompositeZ + BRM evidence
plots/
  *.individual_plots.html          # Per-statistic charts with thresholds
  *.composite_signal.html          # CompositeZ + MaxAbsZ overview
genespace/  (if annotation enabled)
  *.genes_in_qtls.tsv              # Genes overlapping QTL regions
```

### Key Columns

**Raw TSV**: `CHROM`, `POS`, `REF`, `ALT`, genotypes, bulk allele depths (`*BulkAD` = REF,ALT of the biallelic site; `*BulkOther` = reads on any other ALT of a split multi-allelic record), `HighSI`/`LowSI` (selection index), `DeltaSI`, `Gstat`, `ED4`, `HMMLOD` (when computed), `HighDepth`/`LowDepth` (two-bulk), `Depth` (one-bulk: the bulk's depth; two-bulk: effective depth 2·dH·dL/(dH+dL), used for smoothing weights). Two-bulk thresholds are simulated at each site's own (HighDepth, LowDepth) pair. One-bulk modes add `P0`, `AFDev`, `Gstat1`, `HMMLOD`. `HMMLOD` is absent in `2b` without `--ref-parent`.

**Smoothed TSV**: Adds `Sm_*` (smoothed values), `Z_*` (robust Z-scores), `CompositeZ` (Stouffer), `MaxAbsZ`. `HMMLOD` is carried over from the raw TSV unchanged (no Z).

**Merged QTLs**: `CHROM`, `START`, `STOP`, `SOURCE` (`ZScore`/`BRM`/`ZScore+BRM`), `Z_PEAK`, `BRM_PEAK`, `BRM_THRESHOLD`.

---

## Pipeline Overview

```
1. Hard Filter         → Remove low-quality variants (SNP/INDEL-specific thresholds)
2. Raw Statistics      → Per-variant SI, ΔSI, G-stat, ED⁴ (one-bulk: AFDev, G1); pool-HMM LOD along
                          each chromosome (not in 2b without --ref-parent)
3. Smooth & Normalise  → Kernel smoothing, robust Z-scores (MAD), Stouffer composite Z
                          (|ZΔSI| + ZG′)/√2; one-bulk (|ZAFDev| + ZG1)/√2; the HMM LOD is carried, not smoothed
                          SI, ΔSI and AFDev use a Gaussian kernel (σ = window/2), switching to a local linear fit
                          within 2σ of chromosome ends (blended over 2–3σ; only on chromosomes spanning ≥ 6σ)
                          to remove edge bias; markers sharing a position (split multi-allelic views) count as
                          one marker; G and G1 use a tricube kernel
                          over ±window (Magwene's G′); ED⁴ uses loess (local quadratic, tricube weights, ±window,
                          as in MMAPPR; clamped at 0), which has no edge bias at chromosome ends.
                          Markers are weighted by the inverse variance of their allele-frequency estimate,
                          1 / Σ_bulks [1/(c·n) + 1/depth] (c = population variance scale, n = bulk size), so
                          deep markers gain little once pool sampling dominates. A marker is dropped when the
                          effective number of markers in its kernel, (Σw)²/Σw² (`NEff` column), is below 5.
                          On chromosomes with more markers than grid points, kernels are evaluated on a grid
                          (--step-size, capped at σ/20) and interpolated to markers; markers in sparse grid
                          intervals (≤ 4 markers) are evaluated exactly, and dense stretches enter the kernel as
                          σ/100-wide bins at their weighted centroids (≤ 1% of a track's range from exact on
                          pure-noise input). Cost then grows with chromosome length, not marker density
                          (2M markers: smoothing 114 s → ~1 s). Chromosomes are processed one at a time.
                          In 2b without --ref-parent, SI is the REF fraction: its plots are labelled as such,
                          drawn without thresholds. In two-bulk modes HighSI/LowSI are smoothed but not
                          normalised (no Z_HighSI/Z_LowSI): ΔSI carries the contrast.
4. Thresholds          → **Two-stage Monte Carlo simulation** for per-statistic thresholds per depth + bulk size;
                          genome-wide null simulation (--null-reps genomes) for the HMM LOD and CompositeZ
5. BRM Detection       → Bayesian regression model block detection (incorporates bulk size)
6. QTL Detection       → Regions exceeding the empirical CompositeZ threshold
7. Individual QTLs     → Per-statistic intervals using depth-specific thresholds
8. Merge               → Union of Z-score and BRM intervals
9. Plots               → Interactive HTML charts (echarts) with all thresholds
10. Gene Space         → Annotate QTLs with genes (optional, requires --gff)
```

### Statistical Methods

| Method | Formula / Description |
|--------|----------------------|
| **Selection Index (SI)** | `freq(high allele) / freq(all alleles)` per bulk |
| **ΔSI** | `SI_high − SI_low`; large values suggest QTL |
| **G-statistic** | `2 Σ oᵢ ln(oᵢ/eᵢ)` over the 2×2 allele-count table, no pseudo-counts (Magwene et al. 2011) |
| **ED⁴** | `4·ΔSI⁴`, the fourth power of the Euclidean distance between the bulks' allele-frequency vectors (Hill et al. 2013) |
| **HMM LOD** | Pool hidden-Markov-model LOD: log₁₀ of the likelihood of all reads on the chromosome with the pool frequency at the site free vs. under no QTL; genome-wide thresholds from simulated no-QTL genomes |
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

**Stage 1**: Sample realised population-specific diploid genotype dosage in the bulk
```
Pool AF ~ cross-specific genotype sampling (F2/F3/F4/RIL/backcross)
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
| **Composite Z-score** | Empirical Monte Carlo threshold | **Primary** | `*.compositez.qtl.xlsx` |
| **Individual statistics** | Empirical depth-specific thresholds | **Secondary** | `*.individual.qtl.xlsx` |
| **BRM Blocks** | Analytical window-based | **Secondary** | `*.brm_blocks.tsv` |

**CompositeZ and individual-statistic detection are run by default.** The threshold TSV records the exact threshold used for every retained marker.

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
# Mode: hphb — computes AF deviation, one-bulk G and the HMM LOD
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
4. Edwards & Gifford (2012). High-resolution genetic mapping with pooled sequencing. *BMC Bioinformatics*, 13(Suppl 6), S8. (pool HMM LOD, MULTIPOOL)
6. Rousseeuw & Croux (1993). Alternatives to the median absolute deviation. *JASA*, 88(424), 1273–1283.
7. Stouffer et al. (1949). *The American Soldier* (Vol. 1). Princeton University Press.

---

## License

See [LICENSE](LICENSE).

**Issues & Feature Requests**: [GitHub Issues](https://github.com/gmaffy/GoBSAseq/issues)
