# GoBSAseq

**Bulk segregant analysis (BSA-seq) for QTL mapping, in Go**

GoBSAseq takes a multi-sample VCF of parents and phenotypic bulks, hard-filters it, computes allele-frequency statistics per marker, smooths them along each chromosome, and calls QTLs against **genome-wide thresholds from simulated QTL-free genomes**. Each QTL is reported with a calibrated **support interval** that localises it, and Block Regression Mapping (BRM) is run alongside as an independent method.

## Quick Start

```bash
# Build (Go 1.26.4+)
git clone https://github.com/gmaffy/GoBSAseq.git
cd GoBSAseq
go build -o gobsaseq

# Two parents, two bulks of 50 F2 plants each
./gobsaseq \
  -V my_variants.vcf.gz \
  -P parent_high,parent_low \
  -B bulk_high,bulk_low \
  -S 50,50 \
  -m F2 \
  -o results

# Try it on the bundled simulated data (true QTL at chr2:15.0 Mb)
./gobsaseq -V Data/sim_F2_strongQTL.vcf.gz -P P1,P2 -B HB,LB -S 50,50 -m F2 -o results_strong
```

Results go to `results/goBSAseqResults/` (`stats/`, `qtls/`, `plots/`, and `gene_space/` if annotation is enabled).

Two settings matter most:

- **`-S` / `--bulk-sizes`**: the number of plants in each bulk. It drives the null simulation that sets every threshold, so set it accurately.
- **`-b` / `--bulks-depth`**: the minimum read depth per bulk, default `40,40`. It drops every marker below that depth, so for bulks sequenced below ~60× lower it (e.g. `-b 10,10`). Low-depth markers are handled correctly because the simulation uses each marker's observed depth.

Python 3 is only needed for `Data/simulate_bsa.py`; snpEff is only needed for gene-space annotation.

---

## Usage

```bash
gobsaseq -V <vcf> -P <high_parent,low_parent> -B <high_bulk,low_bulk> -S <n_high,n_low> [options]

# From BAMs instead of a VCF (variants are called first with genome-whisperer):
gobsaseq -P <parents> -B <bulks> --parents-bams hp.bam,lp.bam --bulks-bams hb.bam,lb.bam -r ref.fa [options]

# Settings can also come from a config file (key: value lines);
# `gobsaseq CreateTemplate` writes a template to config.txt.
gobsaseq -c config.txt [options]
```

The run is fully non-interactive.

### CLI Flags

#### Input and samples

| Flag | Default | Description |
|------|---------|-------------|
| `-V` / `--variant` | — | VCF file (`.vcf`, `.vcf.gz` or BGZF) |
| `-P` / `--parents` | — | Parent sample names `high,low` (`None` for a missing parent) |
| `-B` / `--bulks` | — | Bulk sample names `high,low` (`None` for a missing bulk) |
| `-S` / `--bulk-sizes` | `20,20` | Plants per bulk `high,low`. **Drives the thresholds; set it accurately** |
| `-m` / `--population` | `F2` | Cross: `F2`, `F3`, `F4`, `RIL`, `BC{n}H`, `BC{n}L`, `BC{n}F2H`, `BC{n}F2L` |
| `-o` / `--out` | `.` | Output directory |
| `-c` / `--config` | — | Config file (see `CreateTemplate`) |
| `-t` / `--threads` | `4` | Threads for variant calling from BAMs/reads |
| `-v` / `--verbose` | off | Verbose output |

The parent and bulk sample names are also selected interactively if they are not given.

#### Smoothing, thresholds and QTL calling

| Flag | Default | Description |
|------|---------|-------------|
| `-w` / `--window-size` | `2000000` | Gaussian smoothing window (bp); σ = window / 2 |
| `--rep` | `1000` | Simulated QTL-free genomes used for thresholds |
| `--seed` | `1` | Random seed (same seed + same input = identical results) |
| `--brm-alpha` | `0.01` | BRM genome-wide error rate: chance that a genome with no QTL gets any BRM block |
| `--min-qtl-width` | `0` | Minimum width of a called region (bp); `0` = σ, negative disables |
| `--min-qtl-markers` | `5` | Minimum markers inside a called region; `0` disables |
| `--merge-distance` | `0` | Maximum gap (bp) bridged when merging regions; `0` = σ |
| `-s` / `--step-size` | `100000` | Not used by the current per-marker smoothing |

#### Variant filtering

Filtering always requires `FILTER=PASS` and a QUAL floor, splits multi-allelic records into biallelic ones (unless `--no-split-multiallelic`), and, in two-parent modes, keeps only markers where the parents carry different alleles. Annotation filters apply only when the annotation is present, so the same settings work for GATK and DeepVariant (the caller type is detected and reported). Defaults follow GATK hard-filtering recommendations.

| Flag | Default | Description |
|------|---------|-------------|
| `-p` / `--parents-depth` | `5,5` | Minimum read depth, parents |
| `-b` / `--bulks-depth` | `40,40` | Minimum read depth, bulks (lower it for bulks below ~60×) |
| `--parents-max-depth`, `--bulks-max-depth` | off | Maximum depth caps to drop coverage outliers (`a,b`) |
| `--min-GQ` | `0` | Minimum genotype quality for parents (not applied to bulks) |
| `--light-filtering` | off | BSA-tuned mode: keeps FS/SOR, drops MQRankSum/ReadPosRankSum (unreliable near fixation, i.e. at QTLs) |
| `--min-QUAL-SNP` / `--min-QUAL-INDEL` | `30` / `30` | Minimum QUAL |
| `--min-QD-SNP` / `--min-QD-INDEL` | `2` / `2` | Minimum QualByDepth |
| `--max-FS-SNP` / `--max-FS-INDEL` | `60` / `200` | Maximum FisherStrand |
| `--max-SOR-SNP` / `--max-SOR-INDEL` | `3` / `10` | Maximum StrandOddsRatio |
| `--min-MQ-SNP` | `40` | Minimum RMS mapping quality |
| `--min-MQRankSum-SNP` | `-12.5` | Minimum MQRankSum |
| `--min-ReadPosRankSum-SNP` / `-INDEL` | `-8` / `-20` | Minimum ReadPosRankSum |

`--min-SOR-SNP` and `--min-FS-SNP` are deprecated aliases of `--max-SOR-SNP` and `--max-FS-SNP`.

#### Variant calling from BAMs (optional)

Uses genome-whisperer's pipeline (gVCFs → joint genotyping → hard filter); the resulting VCF then goes through the normal analysis. The read (FASTQ) flags are accepted but read-based runs are not implemented yet.

| Flag | Default | Description |
|------|---------|-------------|
| `--parents-bams`, `--bulks-bams` | — | BAM files `high,low` |
| `--hp-reads`, `--lp-reads`, `--hb-reads`, `--lb-reads` | — | Read pairs `fwd.fq,rev.fq` (not implemented yet) |
| `-r` / `--reference` | — | Reference FASTA |
| `--caller` | `deepvariant` | `deepvariant` or `gatk` |
| `--deepvariant-version` | `1.10.0` | DeepVariant version |
| `--model-type` | `WGS` | DeepVariant model: `WGS`, `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA` |
| `--merger` | `glnexus` | GVCF merger: `glnexus` or `gatk` |
| `--no-merging` | off | Do not merge GVCFs |

#### Gene-space analysis (optional)

Give all of `--snpEffDB`, `--prg` and `--gff` to run it, or none to skip it; a partial set is an error. `--gene-descriptions` is optional (without it the Description column is `NA`). It annotates variants with snpEff and scores every gene in each QTL's **support interval** (genome-whisperer `genespace`), writing an Excel report per interval to `gene_space/`.

| Flag | Default | Description |
|------|---------|-------------|
| `--snpEffDB` | — | snpEff database ID |
| `--gff` | — | GFF3 annotation |
| `--prg` | — | PRG BLAST top-hits table |
| `--gene-descriptions` | — | Optional gene descriptions TSV (`gene_id<TAB>description`) |
| `--trait` | `disease-resistance` | Trait profile used to rank candidate genes |
| `--min-concordance` | `0.9` | Genotype concordance (0–1] for a variant to count as following the high/low split |
| `--promoter-length` | `2000` | Bases upstream of a gene whose variants count toward it (0 disables) |
| `--cds`, `--protein` | — | CDS and protein FASTA |

---

## Analysis Modes

The mode is detected from which samples are given (`None` = not provided):

| Mode | Parents | Bulks | Statistics |
|------|---------|-------|------------|
| `2p2b` | both | both | Two-bulk: ΔSI, G, ED⁴, LOD, Bayes factor |
| `hp2b` / `lp2b` | one | both | Two-bulk |
| `2b` | none | both | Two-bulk; allele direction is unknown per marker, so ΔSI-based statistics lose power (G and LOD do not) |
| `2phb` / `2plb` | both | one | One-bulk: deviation from the expected frequency, one-bulk G, LOD, Bayes factor |
| `hphb` / `hplb` / `lphb` / `lplb` | one | one | One-bulk |

---

## Input: VCF Requirements

A multi-sample VCF (v4.1/4.2), plain, gzip or BGZF, sorted by position, with per-sample **GT** and **AD** (`ref,alt` depths). INFO annotations (`QD`, `FS`, `SOR`, `MQ`, `MQRankSum`, `ReadPosRankSum`) are used when present. Only the parent and bulk columns are read, so large multi-sample VCFs are fine.

```
#CHROM  POS   ID  REF  ALT  QUAL  FILTER  INFO                  FORMAT    parent1      parent2      bulk_high     bulk_low
chr1    1000  .   A    T    60    PASS    QD=20.0;FS=0.0;SOR=0.8  GT:AD:DP  0/0:50,0:50  1/1:0,45:45  0/1:25,25:50  0/1:30,14:44
```

---

## Output

Files are named `GoBSAseq.<mode>.*` under `<out>/goBSAseqResults/`:

```
stats/
  *.hardfiltered.vcf.gz (+ .tbi)   Variants passing the filter, with tabix index
  *.lowqual.vcf.gz                 Variants failing it
  *.raw.tsv                        Per-marker raw statistics
  *.smoothed_and_normalised.tsv    Smoothed statistics, Z-scores, CompositeZ
  *.thresholds.tsv                 Genome-wide thresholds (repeated per marker)
  *.null_genome_extremes.tsv       Per-genome maxima/minima from the null simulation
  *.null_scales.tsv                Null median and scale used for each Z-score
  *.brm_blocks.tsv                 BRM blocks
qtls/
  *.final.qtl.xlsx                 CompositeZ and BRM calls merged (main result)
  *.compositez.qtl.xlsx            CompositeZ calls
  *.individual.qtl.xlsx            Calls from the individual statistics
plots/
  *.individual_plots.html          Each statistic with its thresholds
  *.robust_z_overlay.html          All Z-scores per chromosome
  *.composite_signal.html          CompositeZ and max |Z|
gene_space/                        Gene annotation of QTLs (optional)
```

### Key columns

- **QTL tables**: `Chrom`, then the **support interval** (`Support_Start`, `Support_End`, `Support_PeakPos`, `Support_PeakLOD`), then the above-threshold region (`Region_Start`, `Region_End`/`Region_Stop`), then per-method peaks. The final table adds `Source` (`CompositeZ`, `BRM` or `CompositeZ+BRM`), `CompositeZ_Pos`/`_Value` and `BRM_Peak`/`BRM_Threshold`.
- **BRM blocks**: `CHROM`, `START`, `STOP`, `PEAK_POS`, `PEAK`, `THRESHOLD`, `SUPPORT_START`, `SUPPORT_END`, `SUPPORT_PEAK_POS`, `SUPPORT_PEAK_LOD`.
- **Raw TSV**: `CHROM`, `POS`, `REF`, `ALT`, parent/bulk genotypes and allele depths, `HighSI`/`LowSI`, `DeltaSI`, `Gstat`, `ED4`, `LOD`, `BBLogBF`, `Depth`. One-bulk modes add `P0`, `AFDev`, `Gstat1`, `LOD1`, `BBLogBF1`.
- **Smoothed TSV**: `Sm_*` (smoothed values), `Z_*` (Z-scores against the null), `CompositeZ`, `MaxAbsZ`.

---

## Pipeline Overview

```
1. Hard filter      Column-subset the VCF (parallel BGZF decompression), split multi-allelics,
                    filter, write filtered VCF + tabix index
2. Raw statistics   Per marker: SI per bulk (polarised by the parents' allele depths), ΔSI, G, ED⁴, LOD, Bayes factor
3. Smoothing        Depth-weighted Gaussian kernel per chromosome (σ = window/2, 3σ cutoff)
4. Thresholds       Null-genome simulation → Z-score scales, CompositeZ, genome-wide P99/P95,
                    BRM multiplier, support-interval drop
5. BRM              Block Regression Mapping blocks
6. Individual QTLs  Regions where any statistic exceeds its genome-wide P99
7. CompositeZ QTLs  Regions where CompositeZ exceeds its genome-wide P99
8. Merge            Union of CompositeZ and BRM regions → final table
9. Plots            Interactive HTML (ECharts): thresholds, BRM blocks (amber), support intervals (red)
10. Gene space      Annotate support intervals (optional)
```

### Statistics

| Statistic | Definition |
|-----------|------------|
| **SI** (SNP-index) | Fraction of a bulk's reads carrying the high-parent allele |
| **ΔSI** | SI(high bulk) − SI(low bulk) |
| **G** | Likelihood-ratio G-statistic for a difference in allele counts between bulks |
| **ED⁴** | Euclidean distance between the bulks' allele frequencies, to the 4th power |
| **LOD** | log₁₀ likelihood ratio that the two bulks differ in allele frequency |
| **Bayes factor** | Beta-binomial log₁₀ BF with Beta(0.5, 0.5) priors |
| **Z-score** | `(x − median) / (1.4826 × MAD)`, with median and MAD from 50 simulated QTL-free genomes (not from the data, so a QTL cannot shrink its own score) |
| **CompositeZ** | Stouffer combination `(|Z_ΔSI| + Z_G + Z_LOD) / √3` (one bulk: AF deviation, one-bulk G, one-bulk LOD) |

### Populations

Expected high-parent allele frequency p₀: `F2`, `F3`, `F4`, `RIL` = 0.5; `BC{n}H` = 1 − 0.5ⁿ⁺¹ (BC1H 0.75, BC2H 0.875); `BC{n}L` = 0.5ⁿ⁺¹; `BC{n}F2H`/`BC{n}F2L` (backcross then selfed) have the same p₀ as `BC{n}H`/`BC{n}L`.

---

## Significance Thresholds (Null-Genome Simulation)

Under the null, smoothed statistics are driven mostly by **drift**: each bulk holds a limited number of plants, so long parental segments are over- or under-represented by chance and linked markers deviate together. Thresholds therefore come from simulating whole QTL-free genomes, not independent SNPs. For each of `--rep` genomes:

1. Every bulk plant is simulated through the cross (`-m`) chromosome by chromosome, with Haldane crossovers. Each chromosome is assumed to be **100 cM** long.
2. At every real marker the pool allele frequency is the simulated dosage, and reads are drawn as `Binomial(observed depth, pool frequency)`, so the real marker density and depths are used.
3. The pipeline's own statistics, smoothing, Z-score transform and CompositeZ formula are applied.
4. The genome-wide maximum (and minimum, for signed tracks) of each statistic is recorded.

**P99** is the 99th percentile of those per-genome maxima: only 1% of QTL-free genomes produce a call anywhere, so it controls the **genome-wide** false-positive rate. Signed statistics (SI, ΔSI, AF deviation) split it across both tails. QTLs are called at P99; P95 is reported for reference.

Smaller bulks mean more drift and higher thresholds, and some weak QTLs in small bulks are genuinely indistinguishable from drift. Larger bulks are the most effective way to gain power.

On simulated F2 data (bulks of 50, 60×), CompositeZ, the individual statistics and BRM made **no false calls in 30 QTL-free genomes**, found 20 of 20 strong and 8 of 10 weak QTLs, and the threshold step simulates 1,000 genomes in about a minute for a tomato-sized genome (900k markers, 16 cores). See `CHANGES.md` for the benchmarks.

### QTL detection methods

| Method | Rule | Output |
|--------|------|--------|
| **CompositeZ** | CompositeZ above its genome-wide P99 | `*.compositez.qtl.xlsx` |
| **Individual statistics** | Any of ΔSI, G, ED⁴, LOD, Bayes factor above its genome-wide P99 | `*.individual.qtl.xlsx` |
| **BRM** | Block where \|ΔSI\| ≥ u·σ (see below) | `*.brm_blocks.tsv` |

The final table (`*.final.qtl.xlsx`) merges CompositeZ and BRM regions. Regions narrower than `--min-qtl-width` or with fewer than `--min-qtl-markers` markers are dropped by every method.

### BRM (Block Regression Mapping)

Following Huang et al. (2020), a block is significant where |ΔSI| ≥ u·σ (one bulk: |SI − p₀| ≥ u·σ), with σ the standard deviation caused by sampling plants into bulks:

- σ² = p(1−p)·[1/(V·n₁) + 1/(V·n₂)] for two bulks, or p₀(1−p₀)/(V·n) for one bulk;
- V is the population variance factor (F2 = 2, F3 = 4/3, F4 = 8/7, RIL = 1, backcrosses by formula);
- u is taken from the null simulation as the (1 − `--brm-alpha`) quantile of each genome's largest |ΔSI|/σ, so `--brm-alpha` is a genome-wide error rate and read noise at your depths is accounted for.

### Localising a QTL: support intervals

The threshold decides **whether** a region holds a QTL, not **where**: with a strong QTL in an F2, linkage spans tens of Mb and the above-threshold region can be the whole chromosome. Each QTL therefore also gets a **support interval**. Starting from the highest smoothed LOD in the region, it extends in both directions while LOD stays within a fixed drop of the peak.

The drop is **6.5 × the null LOD scale** (the `LOD` row of `*.null_scales.tsv`; the value in LOD units is printed at the threshold step). Smoothed LOD grows with read depth, so a fixed drop in raw LOD units does not carry over between depths. Measured against the null scale it does.

Calibration used 216 simulated F2 datasets (depth 30/60/150×, bulks of 30/50/100, strong to weak QTLs, default window). The interval contained the true QTL in **95%** of them (92–98% per condition), with a median width of 2.3 Mb; a separate held-out set gave 93%. On `Data/sim_F2_strongQTL.vcf.gz` the above-threshold region is 0–30 Mb and the support interval 14.93–15.94 Mb (truth 15.0 Mb). Other crosses and window sizes have not been calibrated.

**Worked example** (`Data/sim_F2_strongQTL.vcf.gz`, chr2, true QTL at 15.0 Mb):

1. **Peak.** The highest smoothed LOD in the called region is 32.817, at 15,436,548.
2. **Drop.** The null LOD scale is 0.0879 (`*.null_scales.tsv`), so the drop is 6.5 × 0.0879 = 0.571 LOD (printed as `support-interval drop 0.571 LOD`).
3. **Floor.** 32.817 − 0.571 = 32.246.
4. **Walk outwards.** From the peak, move marker by marker in each direction and stop at the last marker whose LOD is still ≥ 32.246. On the left that is 14,929,613 (LOD 32.254); the next marker, 14,916,794, is at 32.225, below the floor. On the right it is 15,942,346.
5. **Result.** Support interval 14,929,613–15,942,346 (1.0 Mb), containing the true position; the above-threshold region was the whole chromosome.

**How to read it.** Treat the support interval roughly as a **95% interval for the QTL position**: in about 1 case in 20 the true position lies just outside it. Its width reflects how sharply the signal peaks, so strong QTLs give narrow intervals (about 1 Mb in the simulations) and weak QTLs wide ones (several Mb). The interval is centred on the LOD peak, which can differ slightly from the CompositeZ or BRM peak positions also reported.

### Reading the plots

Each chart (individual statistics, Z-score overlay, composite signal) shows one chromosome:

| Element | Meaning |
|---------|---------|
| Solid line | The smoothed statistic (or CompositeZ) along the chromosome |
| Dashed red / orange horizontal lines | Genome-wide **P99** / **P95** thresholds; a QTL is called where the line crosses P99 |
| **Amber** band | **BRM block**: region where Block Regression Mapping finds the allele-frequency difference significant. Like the above-threshold region, it can span most of a chromosome; it says *that* the chromosome carries a QTL |
| **Red** band, labelled `QTL <start>–<end> Mb` | **Support interval** from the final QTL table: *where* the QTL most likely is (≈95% interval, see above). Gene-space analysis uses this interval |

Red is drawn on top of amber. A strong QTL typically looks like a wide amber band with a narrow red band at the peak. Red without amber means CompositeZ called the QTL but BRM did not. The slider under each chart zooms into a region, and hovering shows values (★ marks points above P99, ● above P95; values below the lower thresholds are not marked).

---

## Test Data

`Data/` holds three simulated F2 datasets (strong QTL, weak QTL, no QTL) with their expected results, and `simulate_bsa.py`, which generates more with a known truth (bulk size, depth, marker spacing, QTL position, and strength `strong`/`moderate`/`weak`/`none`). See `Data/README.md`.

---

## Troubleshooting

| Problem | Likely cause | Solution |
|---------|--------------|----------|
| Very few markers pass filtering | Bulk depth below `--bulks-depth` (default 40) | Lower `-b` (e.g. `-b 10,10` for ~30× bulks) |
| 0 variants pass filtering | Strict thresholds or non-standard INFO | Relax `--min-QUAL-*`, `--min-QD-*`; check FILTER/INFO |
| `not part of the VCF sample list` | Sample name mismatch | `bcftools query -l file.vcf` |
| No QTLs detected | Weak QTL relative to drift | Check `CompositeZ` and `Sm_Gstat` in the smoothed TSV against the thresholds; larger bulks give more power |
| QTL region spans a whole chromosome | Expected for strong QTLs in early generations | Use the support interval (`Support_*` columns, red shading) |
| Thresholds vary between runs | Different `--seed` or `--rep` | Same seed + same input gives identical results; more `--rep` gives more stable P99 |
| HTML plots won't render | Browser blocks local files | `python3 -m http.server 8000` in the plots folder |

---

## Development

### Project structure

```
main.go                 Entry point
cmd/root.go             CLI (Cobra); cmd/CreateTemplate.go writes a config template
run/run.go              Pipeline orchestration
filter/
  filter.go             Hard filtering, compact variant records
  vcfreader.go          Sample-column subsetting, gzip/BGZF decompression
  split.go              Multi-allelic splitting
  tabix.go              Tabix (.tbi) index writer
stats/
  stats.go              Raw statistics
  smoothing.go          Gaussian smoothing, Z-scores, CompositeZ
  nullsim.go            Null-genome simulation (cross, recombination, reads, smoothing)
  thresholds.go         Genome-wide thresholds and null scales
  brm.go                Block Regression Mapping
  support.go            Support intervals
  regionfilter.go       Minimum region width / marker count
  detect.go             QTL detection, merging, Excel output
  genespace.go          Gene annotation (optional)
plots/plots.go          HTML plots (go-echarts)
utils/                  Configuration and helpers
Data/                   Simulated test data and simulator
```

### Build and test

```bash
go build -o gobsaseq
go test ./...
go vet ./...
go test ./stats -run XXX -bench NullGenome   # null-simulation benchmark at plant genome sizes
```

`CHANGES.md` documents recent changes with before/after benchmarks.

---

## References

1. Michelmore RW, Paran I, Kesseli RV (1991). Identification of markers linked to disease-resistance genes by bulked segregant analysis. *PNAS* 88(21):9828–9832.
2. Magwene PM, Willis JH, Kelly JK (2011). The statistics of bulk segregant analysis using next generation sequencing. *PLoS Comput Biol* 7(11):e1002255. doi:10.1371/journal.pcbi.1002255
3. Takagi H et al. (2013). QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. *Plant J* 74(1):174–183. doi:10.1111/tpj.12105
4. Huang L, Tang W, Bu S, Wu W (2020). BRM: a statistical method for QTL mapping based on bulked segregant analysis by deep sequencing. *Bioinformatics* 36(7):2150–2156. doi:10.1093/bioinformatics/btz861
5. Haldane JBS (1919). The combination of linkage values and the calculation of distances between the loci of linked factors. *J Genet* 8:299–309.
6. Kass RE, Raftery AE (1995). Bayes factors. *JASA* 90(430):773–795.
7. Rousseeuw PJ, Croux C (1993). Alternatives to the median absolute deviation. *JASA* 88(424):1273–1283.
8. Stouffer SA et al. (1949). *The American Soldier*, Vol. 1. Princeton University Press.

---

## License

See [LICENSE](LICENSE).

**Issues and feature requests**: [GitHub Issues](https://github.com/gmaffy/GoBSAseq/issues)
