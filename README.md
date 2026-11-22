# GoBSAseq

**A high-performance Bulk Segregant Analysis (BSAseq) pipeline in Go**

GoBSAseq identifies genomic regions associated with phenotypic traits by applying hard filtering to VCF files, computing allele frequency statistics and kernel smoothing them, and calling QTLs with a pool hidden-Markov-model LOD, block regression mapping (BRM; Huang et al. 2020) and DsiGconsensusZ (a consensus of ΔSI and G′), each against genome-wide thresholds from simulated no-QTL genomes.

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
| `--brm-alpha` | float64 | `0.01` | BRM genome-wide error rate: the chance that a genome without QTLs gets any BRM block. Also BRM's level in the final QTL overlap, matching the HMM LOD and DsiGconsensusZ p99 |
| `--merge-distance` | int | `0` | Gap in bp within which same-direction regions are joined (0 = smoothing σ, window/2) |
| `--min-qtl-markers` | int | `5` | Fewest markers beyond the threshold line a region needs (individual statistics, DsiGconsensusZ, BRM) |
| `--rep` | int | — | Deprecated and ignored: all thresholds come from the genome-wide null (`--null-reps`) |

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
| `--null-reps` | `200` | Simulated no-QTL genomes for every genome-wide threshold (HMM LOD, BRM, DsiGconsensusZ, smoothed tracks); DsiGconsensusZ's null model uses another `max(50, null-reps/2)` |
| `--no-null-cache` | off | Recompute the genome-wide null even if `stats/*.null_cache.gob` holds one for the same markers, depths and settings (a re-run that only changes QTL-calling, plotting or gene-space options reuses it) |
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
| `--gene-space-span` | Gene space analysis on each final QTL's span instead of its consensus core |
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

**`2b` polarity.** Without a parent sample, SI is the REF-allele fraction. Unless the reference genome is one of the parents, REF comes from the high parent at some SNPs and the low parent at others, so the sign of ΔSI is arbitrary per SNP and would cancel when smoothed. With `--ref-parent none` (the default), `2b` therefore uses \|ΔSI\| (columns `AbsDeltaSI`, `Sm_AbsDeltaSI`, `NZ_AbsDeltaSI`): its thresholds are upper-tail only, there is no `DeltaSI-` track, DsiGconsensusZ's ΔSI component is the polarity-free BRM signal, and BRM uses the polarity-free form described below. QTLs can be located, but not which bulk carries which parent's allele. G and ED⁴ do not depend on polarity. The HMM LOD needs polarity and is not computed in this case, so final QTLs rest on BRM and DsiGconsensusZ, which both test the ΔSI contrast: one line of evidence (the log says so). BRM uses a polarity-free form. It smooths the bias-corrected ΔSI², that is (SI_H − SI_L)² minus each bulk's read-sampling variance SI(1−SI)/(d−1), which estimates the squared pool frequency difference whatever allele is REF. BRM then tests its square root (\|ΔAF\|) against the usual threshold, with af at the population's expected frequency. Blocks have no direction. On simulated no-QTL genomes its false-block rate matched the signed BRM's (24 vs 22 of 100 genomes). If the reference is a parent, pass `--ref-parent high` or `low` to get signed ΔSI and BRM.
**One-bulk modes** compute: AF deviation, one-bulk G-stat and the pool-HMM LOD.

**HMM LOD** (after MULTIPOOL, Edwards & Gifford 2012). Each bulk is modelled as a pool whose hidden high-parent allele count changes along the chromosome only through crossovers, so all linked SNPs inform the same pool frequency, and pool sampling of the plants is part of the model. The pool prior and linkage depend on the population type (F2, F3, F4, RIL, BCn, BCnF2); recombination defaults to 1 Morgan per chromosome (`--cm-per-mb` sets a genome-wide rate); markers within 0.05 cM share one HMM step (against 0.01 cM: 3× faster on dense markers, the same peak and genome-maximum LOD on simulated data, support-interval ends moved by ≤ 10 kb for strong QTLs and up to 0.8 Mb on the flat plateaus of moderate ones). Two-bulk LOD tests whether the bulks differ, with their shared frequency free, so segregation distortion that shifts both bulks is not scored as a QTL. One-bulk LOD tests the bulk against P0 (with the mapping-bias correction). It is computed in the raw-statistics stage and written to the raw TSV. The LOD is already linkage-aware, so it is not smoothed or normalised and is not part of DsiGconsensusZ: it stays in LOD units, as its own track with its own thresholds. Its thresholds are genome-wide: the 95%/99% quantiles of the maximum LOD over `--null-reps` simulated no-QTL genomes at the real markers and depths, run through the same scoring code. The same simulated genomes give the DsiGconsensusZ thresholds, also genome-wide (the quantiles of each genome's maximum DsiGconsensusZ), so a 99% threshold means a 1% chance that a genome without QTLs has any marker above it. Every scanned statistic gets its threshold the same way (see Thresholds below). Quantiles come from a generalized extreme-value fit to the genome maxima (probability-weighted moments), which is stable at 200 genomes and allows tails heavier or lighter than Gumbel's. The simulated reads carry the read-level overdispersion ρ estimated from the data (beta-binomial; from tightly linked neighbouring markers, whose difference has no pool or QTL component; single sites in `2b` without `--ref-parent`), so overdispersed data are not scored against a binomial null. ρ is printed in the log. The HMM LOD replaces the earlier per-SNP LOD and beta-binomial Bayes factor, which ignored pool sampling and linkage.

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
  *.smoothed_and_normalised.tsv    # Smoothed tracks + null-standardised NZ_* columns and DsiGconsensusZ
  *.thresholds.tsv                 # Genome-wide thresholds (one row per statistic and tail), with ρ, --null-reps and DsiGconsensusZ's r
  *.null_cache.gob                 # Cached genome-wide null, reused by a re-run on the same data and settings
  *.brm_blocks.tsv                 # BRM blocks: CHROM, START, STOP, PEAK_POS, PEAK, THRESHOLD (u·σ at the peak), PEAK_Z (signal/σ)
  *.hardfiltered.vcf.gz (+.tbi)    # Filtered VCF with tabix index
  *.lowqual.vcf.gz                 # Records where no ALT passed (kept samples only)
qtls/
  *.individual.qtl.xlsx / .tsv     # Individual-statistic QTLs (by direction)
  *.dsigconsensusz.qtl.xlsx / .tsv # DsiGconsensusZ QTL intervals (genome-wide threshold)
  *.final.qtl.xlsx / .tsv          # Final QTLs: consensus of HMM LOD, BRM and DsiGconsensusZ (plus loci with partial support)
plots/
  *.individual_plots.html          # Per-statistic charts with thresholds (both tails), BRM blocks shaded
  *.dsigconsensusz.html            # DsiGconsensusZ with its genome-wide thresholds
genespace/  (if annotation enabled)
  *.genes_in_qtls.tsv              # Genes overlapping QTL regions
```

### Key Columns

**Raw TSV**: `CHROM`, `POS`, `REF`, `ALT`, genotypes, bulk allele depths (`*BulkAD` = REF,ALT of the biallelic site; `*BulkOther` = reads on any other ALT of a split multi-allelic record), `HighSI`/`LowSI` (selection index), `DeltaSI`, `Gstat`, `ED4`, `HMMLOD` (when computed), `HighDepth`/`LowDepth` (two-bulk), `Depth` (one-bulk: the bulk's depth; two-bulk: effective depth 2·dH·dL/(dH+dL), used for smoothing weights). Two-bulk thresholds are simulated at each site's own (HighDepth, LowDepth) pair. One-bulk modes add `P0`, `AFDev`, `Gstat1`, `HMMLOD`. `HMMLOD` is absent in `2b` without `--ref-parent`.

**Smoothed TSV**: Adds `Sm_*` (smoothed values); `NZ_DeltaSI` (`NZ_AbsDeltaSI`, one-bulk `NZ_AFDev`), the ΔSI-family z = signal/σ with unit null scale; `NZ_Gstat` (one-bulk `NZ_OneBulkGstat`), G′ standardised against the null at the marker, as a normal score; and `DsiGconsensusZ`. `HMMLOD` is carried over from the raw TSV unchanged.

**DsiGconsensusZ** replaces the old robust-Z CompositeZ, whose scale came from the genome's own median and MAD (genome-wide line about 20–25, different in every dataset) and which vetoed QTLs that the HMM LOD and BRM both called. Each component is now standardised against the simulated no-QTL genomes at every marker: the ΔSI family by its model null SD (BRM's z), G′ by its null mean and SD there. Each becomes a one-sided normal score through its pooled null distribution, and the two are combined by Stouffer's method corrected for their null correlation r: (s₁ + s₂)/√(2 + 2r). The genome-wide line is about 4. G′ and ΔSI come from the same reads, so r is about 0.9 (printed in the log and the thresholds TSV): the consensus carries little beyond ΔSI, and the two agreeing is not independent confirmation. On simulated F2 data (20-plant bulks) a moderate QTL reached 0.72 of the DsiGconsensusZ line, against 0.36 of the CompositeZ line.

**Individual QTLs**: each statistic is scanned across the genome at its p99 threshold; only if it has no p99 region anywhere are its p95 regions reported, labelled `p95 (suggestive)`. Regions are joined into QTLs by direction: ΔSI/AFDev upper-tail (`+`, high-parent allele enriched in the high bulk) and lower-tail (`−`) regions are never joined, even when they overlap; unsigned statistics (G′, ED⁴, |ΔSI|, G1, HMM LOD) join the signed QTL they overlap most, or form `none` QTLs. Each QTL lists each statistic's strongest peak.

**Final QTLs**: a locus is formed by HMM LOD regions (above the genome-wide p99), BRM blocks (`--brm-alpha`, default 0.01) and DsiGconsensusZ regions (above the genome-wide p99) that overlap or lie within `--merge-distance`. It is a **consensus** QTL when every available method supports it: all three, or BRM and DsiGconsensusZ in `2b` without `--ref-parent` (no HMM LOD). Columns:
- `Start`/`Stop`: the core, where all methods are above threshold. `Span_Start`/`Span_Stop`: the union of their regions.
- `Consensus` (yes/no), `Support` (e.g. `HMMLOD+BRM+DsiGconsensusZ`; loci with partial support are listed but are not final).
- `Evidence`: the independent lines of evidence behind the locus, `HMM` and/or `ΔSI-family` (BRM and DsiGconsensusZ test the same allele-frequency contrast, so together they are one line).
- `Direction`: `+`/`−` from the BRM peak, `none` without polarity, `mixed` if BRM blocks of both signs fall in the locus.
- `Lead_Pos`: the midpoint of the HMM LOD's 1.5-LOD support interval (else the BRM peak). Then each method's peak, value and threshold (BRM also z), the HMM support interval, and the individual statistics overlapping the locus.

Gene space analysis runs on the consensus QTLs' cores (`--gene-space-span` for spans). On simulated data no QTL-free genome out of 100 produced a consensus QTL.

---

## Pipeline Overview

```
1. Hard Filter         → Remove low-quality variants (SNP/INDEL-specific thresholds)
2. Raw Statistics      → Per-variant SI, ΔSI, G-stat, ED⁴ (one-bulk: AFDev, G1); pool-HMM LOD along
                          each chromosome (not in 2b without --ref-parent)
3. Thresholds          → Genome-wide null simulation (--null-reps no-QTL genomes at the real markers,
                          depths and bulk sizes, with the estimated read overdispersion), scored by the
                          same code as the data. Every scanned statistic (smoothed ΔSI, G′, ED⁴, AFDev, G1,
                          SI; HMM LOD; BRM; DsiGconsensusZ) gets genome-wide thresholds: quantiles of its genome
                          maximum (and minimum for lower tails; two-tailed tracks split the level). A further
                          max(50, null-reps/2) genomes fit DsiGconsensusZ's null standardisation. Replicates run
                          in parallel (results do not depend on the core count); the null is cached.
4. Smooth & Normalise  → Kernel smoothing; ΔSI and G′ (one-bulk AFDev and G1) standardised against the null
                          at every marker and combined into DsiGconsensusZ; the HMM LOD is carried, not smoothed
                          SI, ΔSI and AFDev use a Gaussian kernel (σ = window/2), switching to a local linear fit
                          within 2σ of chromosome ends (blended over 2–3σ; only on chromosomes spanning ≥ 6σ)
                          to remove edge bias; markers sharing a position (split multi-allelic views) count as
                          one marker; G and G1 use a plain tricube kernel
                          over ±window (Magwene's G′, as QTLseqr: no depth weights, since per-site G already
                          grows with depth); ED⁴ uses loess (local quadratic, plain tricube weights, ±window,
                          as in MMAPPR; clamped at 0), which has no edge bias at chromosome ends.
                          SI, ΔSI and AFDev markers are weighted by the inverse variance of their allele-frequency estimate,
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
                          drawn without thresholds.
5. BRM Detection       → Block Regression Mapping (Huang et al. 2020): runs of markers where the smoothed
                          ΔSI / AFDev exceeds u·σ of the smoothed value (pool + read variance); u from the
                          genome-wide null, so --brm-alpha is a genome-wide error rate
6. Individual QTLs     → Each statistic's regions beyond its genome-wide p99 (p95 "suggestive" only if the
                          statistic has no p99 region anywhere), joined by direction (+/−/none)
7. DsiGconsensusZ QTLs → Regions above the genome-wide DsiGconsensusZ p99 (same p95 rule)
8. Final QTLs         → Consensus of HMM LOD, BRM and DsiGconsensusZ regions (core, span, lead position, evidence)
9. Plots               → Interactive HTML charts (echarts) with all thresholds: position on a linear bp axis,
                          tracks drawn as tested (no spline), ★/● for values beyond a p99/p95 line of either
                          tail. Each track is simplified to the points needed for the line to pass within
                          0.1% of its range of every marker (exact extremes kept), values are written to
                          tooltip precision and thresholds as mark lines: ~25× smaller pages
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
| **NZ components** | ΔSI family: `z = signal/σ` (BRM's z) rescaled to unit null RMS; G′: `(G′ − E₀)/SD₀` at the marker from the null genomes. Each → normal score `Φ⁻¹(F₀(v))` via its pooled null distribution |
| **DsiGconsensusZ** | Correlation-corrected Stouffer: `(s₁ + s₂)/√(2 + 2r)`, r = null correlation of the two scores; genome-wide thresholds from the null |
| **BRM** | Threshold = `u·σ`, u the (1 − `--brm-alpha`) quantile of the null genomes' maximum \|z\| = \|signal\|/σ (extreme-value fit; about 2.9 at α = 0.05, higher at the default 0.01), σ² = pool variance `p₀(1−p₀)(1/n₁+1/n₂)/V_scale` + the read-sampling variance of the smoothed value, `Σ l_m²·r_m` over the kernel's effective weights with beta-binomial read variance (estimated ρ) |
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
| **DsiGconsensusZ** | Genome-wide null simulation (extreme-value fit to genome maxima) | **Primary** | `*.dsigconsensusz.qtl.xlsx` |
| **Individual statistics** | Genome-wide null simulation (per statistic) | **Secondary** | `*.individual.qtl.xlsx` |
| **BRM Blocks** | Standardised smoothed signal, genome-wide null threshold | **Secondary** | `*.brm_blocks.tsv` |

**DsiGconsensusZ and individual-statistic detection are run by default.** The threshold TSV records every threshold used: they are genome-wide, so a p99 threshold has a 1% chance of being crossed anywhere in a genome without QTLs.

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
  -m F2 --brm-alpha 0.05 \
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
  -m F2 --null-reps 500  # More null genomes: tighter thresholds, longer run
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
| No QTLs detected | Low power or weak signal | Check max DsiGconsensusZ and HMMLOD in the smoothed TSV against the thresholds TSV; partial-support loci are listed in the final table |
| HTML plots won't render | Browser blocks local files | Serve with `python3 -m http.server 8000` |
| Methods disagree | Expected — different methods | Final QTLs need every method; loci with partial support are listed with their `Support` |

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
