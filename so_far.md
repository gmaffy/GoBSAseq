# GoBSAseq: work so far

Status: all changes are in the working tree and **nothing is committed**. `go build`, `go vet` and `go test ./...` pass.

## 1. Codebase review
`claude_review.md` covers biology, statistics and code. It cross-checks the three root reviews (`copilotfusion_review.md`, `gemini_review.md`, `grokreview.md`) and `indep_bsa.py`.

Main conclusions:
- **Population genetics is correct:** `ExpectedAF`, `PopulationVarianceScale` and the BRM variance formula all check out.
- **Significance is not controlled genome-wide.** In a simulated null F2 genome (20+20 plants, 50x, 10 × 30 Mb, σ = 1 Mb):
  - CompositeZ "p99": 30% of null genomes get at least one false QTL (40% at 150x).
  - CompositeZ "p95": 80%.
  - BRM at α = 0.05: 45%.
  - This contradicts gemini_review's claim that the thresholds are too strict.
- **Recommended fix:** a chromosome-scale null simulation, run through the real pipeline, with genome-wide maximum thresholds (the approach of `null_maxima` in `indep_bsa.py`).

Simulation scripts were in the session scratchpad (`czsim.py`, `brmsim.py`) and may no longer exist.

## 2. Fixes applied

### 2a. Bulk sizes and backcross-F2 simulation (`stats/`, `utils/config.go`)
1. **Bulk size per mode:** `cfg.BulkSizes(bsaType)` returns the size each mode uses. `OneBulkSize` was never set, so one-bulk modes had used depth as the number of plants. `lb` modes now use the low bulk's size.
2. **No silent fallbacks:** a missing bulk size is now an error. Previously the code fell back to depth, 100 plants, or the high bulk's size.
3. **BCnF2 simulation:** `samplePoolAF` now samples three genotypes (recurrent hom 1−q+q/4, het q/2, donor hom q/4, with q = 0.5ⁿ) instead of treating BCnF2 as BCn.

Tests in `stats/bulksize_test.go`: the sampled mean and variance match `ExpectedAF` and `PopulationVarianceScale` for 11 populations, and BC1F2 produces donor homozygotes at 1/8.

### 2b. Hard filtering rewritten (`filter/`)
The old vcfgo-based code is **commented out in place**, not deleted:
- `filter/filter.go`, `filter/split.go` and `filter/filter_test.go`;
- the vcfgo helpers in `stats/stats.go` and `utils/utils.go` (`GetFloat`);
- the old reader in `run/run.go` and the `Rdr` field in config.

**Before:** vcfgo parsed every sample of every record in one goroutine, and each passing site kept a full `*vcfgo.Variant`. Splitting multi-allelic records threw away reads for the other ALTs. The `.tbi` was broken (stale offsets, plus two bugs in biogo/hts/tabix).

**Now:** a streaming line parser.
- **`vcfheader.go`:** header parsing and multithreaded BGZF input.
- **`record.go`:** reads only the kept sample columns; splitting keeps the other-ALT depth.
- **`fastfilter.go`:** parallel workers and the filter rules.
- **`bgzfwriter.go`:** parallel BGZF output with exact virtual offsets.
- **`tabix.go`:** its own TBI writer.
- **Output:** a compact `filter.Site` per passing site.
- **Moved:** `ExpectedAF` to `utils/population.go`.

**Performance** (300k sites × 60 samples, 4 used):

| | Before | After |
|---|---|---|
| Time | 13.3 s | 0.58 s |
| Live heap | 1.68 GB | 108 MB |
| Peak RSS | 3.6 GB | ~0.5 GB |

**GATK compliance:**
- SNP and indel hard-filter thresholds as GATK recommends; a value exactly at a threshold passes.
- Missing, `.` or NaN annotations pass (as in `VariantFiltration`).
- FILTER must be PASS or `.`.
- Light mode drops only the rank-sum tests.
- It complied before as well.
- Documented deviations:
  - multi-nucleotide substitutions get SNP thresholds;
  - rejected records go to `lowqual.vcf.gz` rather than being flagged;
  - in a cohort VCF, site annotations cover all samples.

**BSA-seq rules, now:**
- **Bulk minimum depth** = REF + ALT reads from AD (previously FORMAT/DP).
- **Other-ALT reads ≤ 5%** of allele depth: this check now actually runs after splitting.
- **Parent homozygosity (≥ 85%)** counts other-ALT reads against both alleles. A parent with AD present but zero informative reads is not confident.
- **GT fallback** accepts only a homozygous REF or target-ALT genotype.
- **A lone parent must be homozygous.**
- **Segregation floor for 2b/hp2b/lp2b:** pooled minor allele frequency ≥ 0.4 × the expected value (0.2 for F2). `--min-pooled-maf`; `0` turns it off.
- **One parent + one bulk:** the non-parent-allele rule is **off by default** (`--require-nonparent-allele`).
- **Polarity** (`HighIsAlt`) is decided once in the filter.
- **Rejection report:** a count per rule is printed at the end, and the counts sum exactly to the rejected total.

**Verification:**
- 19 filter tests: GATK boundaries, the old behaviours, each fix, and a tabix query test.
- **htslib `tabix`** agreed with a full scan on 300 random regions, including indels.
- **Compared with the old filter** on 50,000 messy records in 4 modes:
  - no site added;
  - every dropped site explained by an intended fix;
  - SNP-index values identical on shared sites.
- **End-to-end CLI run** OK, and bcftools reads the outputs.

Sites kept on the test VCF, old → new:

| Mode | Old | New |
|---|---|---|
| 2p2b | 3,569 | 2,519 |
| 2b | 12,536 | 8,086 |
| hp2b | 11,307 | 5,774 |
| hphb | 11,860 | 7,917 |

The test VCF has exaggerated problem rates, so expect smaller drops on real data.

**Not yet done:**
- no run on a real VCF yet;
- `go test -race` couldn't run (no C compiler in this environment).

## 3. Raw statistics review (`stats.RawStats`): not yet fixed

| # | Finding | Severity |
|---|---|---|
| 1 | Rounding stored values to 6 dp makes ED⁴ = 0 whenever \|ΔSI\| < 0.027 (14% of null sites). ED⁴ also lacks the literature factor of 4. | Medium |
| 2 | G adds a 0.5 pseudo-count to all cells (not Magwene's standard G). | Low–Medium |
| 3 | "LOD" equals G/(2 ln 10); CompositeZ counts it twice. | Medium |
| 4 | The Bayes factor and G are read-level; only simulated thresholds are valid. | Docs |
| 5 | One-bulk p₀ ignores reference mapping bias; option: p₀ from the genome-wide median SI. | Medium (one-bulk) |
| 6 | 2b SI = REF fraction, so signed ΔSI cancels when smoothed. | High for 2b |
| 7 | Depth = the smaller of the two bulk depths. | Medium |

Code issues:
- `writeRawTSV` ignores `Flush` and `Close` errors.
- A site with the wrong sample count is skipped silently.
- Dead code and substring-based mode logic.

Speed (2M sites: 6.4 s, 4 GB allocated):
- the TSV writer (`fmt.Sprintf`) takes about 70% of the time;
- the progress bar is updated per site under a mutex, with an atomic counter per site and false sharing on `keep[]`;
- GT slices and AD strings are allocated per site.

Proposed next steps:
1. Stop rounding stored values; write ED⁴ with `%g`.
2. Faster TSV writer and chunked worker loop; check write errors.
3. Optional p₀ estimated from the data for one-bulk modes.
4. Standard G and/or drop LOD from CompositeZ (this changes results).

## 4. Open items from the review (not started)
- **Thresholds:** genome-wide null (B1), which also covers unequal depths (B3), tail labels (B4), BRM α (B6) and the seed (B8).
- **2b polarity** (A3).
- **One-bulk BRM:** treats `SmAFDev == 0` as missing (B5).
- **Detection:**
  - p95 peaks are dropped on chromosomes that have a p99 peak (C1);
  - when merging lower-tail peaks, the weaker summit is kept (C2);
  - when peaks merge, a QTL keeps the last or first peak instead of the strongest (C3, C4);
  - chromosome order in the output is random (C5);
  - the smoothed output isn't sorted (C6).
- **CLI:**
  - unused flags (C7);
  - the reads path exits 0 without doing anything (C8);
  - interactive prompts, including when `-P hp,None` is given (C9);
  - FS/SOR flags named `--min-*` (C10);
  - no way to exclude contigs (C11).
- **README:** the Magwene citation and the "LOD > 3" rule.
