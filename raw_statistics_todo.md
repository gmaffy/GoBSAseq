# Raw statistics (`stats.RawStats`): todo

Scope: `stats/stats.go` (per-site SI, ΔSI, G, ED⁴, LOD, BBLogBF and one-bulk equivalents) and the `raw.tsv` writer.
Sources: the code, `so_far.md` §3, `claude_review.md` (A3, A5, D), `grokreview.md` (#4, #5, ED⁴), `gemini_review.md` (2.5, 2.6, 3.6), `copilotfusion_review.md` (#2).

Order: do Part 1 (behaviour-preserving speed and safety) first and lock it with a golden-output test. Then do Part 2 (changes to results) one item at a time, each with its own before/after comparison.

---

## Part 0: test harness first

- [x] **Golden test.** Run `RawStats` on a fixed `[]filter.Site` (all 10 modes, including polarity flips, zero-depth bulks and an all-one-allele bulk). Compare every `BSAstats` field and the TSV byte for byte. Part 1 must not change this output.
- [x] **Benchmark.** `BenchmarkRawStats` on about 2M synthetic sites in `2p2b` and `hphb`, with `-benchmem`. Report the compute time and the TSV write time separately. No benchmark exists yet. The last measurement (`so_far.md`) was 6.4 s and 4 GB allocated for 2M sites, about 70% of it in the TSV writer.
- [x] **Known-value tests** (none exist yet; requested by every review):
  - G on 2×2 tables with known values, including tables with zero cells.
  - LOD = G_nopseudo / (2 ln 10) on the same tables.
  - ED⁴ for a given Δ.
  - One-bulk G1 = 2 ln 10 · LOD1 exactly.
  - BBLogBF against a hand-computed value.
  - Polarity: with `HighIsAlt` set, H and L swap and HighSI becomes the ALT fraction.

## Part 1: speed, memory and code health (no change to results)

**Done (1.1–1.4):** 2M sites, 20 cores — 2p2b 3.8 s → 0.24 s, 3.5 GB → 0.55 GB allocated, 70.8M → ~1k allocations; hphb 3.4 s → 0.28 s. `BSAstats` 352 → 224 bytes. Output byte-identical to the golden files for all modes. Workers take 8192-site chunks, compute and format their rows; the caller writes chunks in order (≤ 2 chunks per worker in flight). Not done: `raw.tsv.gz`/`--no-raw-tsv` (writer no longer dominates), `float32` (would change results), per-chromosome streaming (longer-term refactor, with B1).

### 1.1 TSV writer (`writeRawTSV`, `stats.go:361`)
- [x] Replace the per-field `fmt.Sprintf` + `[]string` + `strings.Join` with a single reused `[]byte` buffer: `strconv.AppendInt` / `AppendFloat(buf, v, 'f', 6, 64)` and `append(buf, '\t')`. Use a 1 MB `bufio.Writer`.
- [x] GT: write directly from `SampleAD` counts (e.g. `[0 1]`) instead of `fmt.Sprintf("%v", []int)`.
- [x] Check errors: currently `Flush` and `Close` are deferred and their errors dropped, and `Fprintln` errors are ignored. A full disk produces a truncated TSV with exit status 0. Use `if err := w.Flush(); …` and `if err := f.Close(); …`.
- [x] Optional: format rows in parallel (chunks of about 64k rows, one buffer per worker, written in order). Worth doing only if the benchmark still shows the writer dominating after the single-buffer change.
- [ ] Optional: write `raw.tsv.gz` (or make the TSV optional with `--no-raw-tsv`). It is the largest output and most users never open it.

### 1.2 Stop storing what only the TSV needs
`HighParGT`, `LowParGT`, `HighBulkGT`, `LowBulkGT` (4 heap slices) and `HighBulkAD`, `LowBulkAD` (2 `Sprintf` strings) are allocated per site, but nothing downstream reads them (`smoothing.go`, `brm.go` and `thresholds.go` use only the numeric fields, `REF`, `ALT` and `Depth`).
- [x] Drop these fields from `BSAstats`. Have `writeRawTSV` read GT and AD from the `filter.Site` (pass `[]filter.Site` plus the kept indices, or write the TSV inside the worker loop from the site). This removes about 6 allocations per site and roughly 130 B from each struct.
- [x] Check `plots/` and `README.md` for any use of those columns before removing them. The TSV columns themselves should stay.

### 1.3 Worker loop (`stats.go:242-341`)
- [x] Hand out chunks (e.g. 4096 sites) instead of one atomic increment per site.
- [x] Update the progress bar once per chunk, not per site. `progressbar.Add` takes a mutex, so with all cores updating it per site it becomes a contention point.
- [x] `keep[]` is a `[]bool` written by every worker per site, which causes false sharing. Chunking removes most of it. Better: drop `keep` entirely (see 1.4).
- [x] Compute `highTotal`/`lowTotal` once. They are recomputed at `stats.go:312-313` and then discarded with `_ =`.

### 1.4 Memory
- [x] `results` and the compacted `stats` slice are both alive at the same time (`stats.go:345`), so peak memory is 2× the slice. Compact in place: `out := results[:0]; for i := range results { if keep[i] { out = append(out, results[i]) } }`. Or, if 1.5 makes the skip a hard error, return `results` as is.
- [x] Measure `unsafe.Sizeof(BSAstats{})` before and after 1.2. Consider `float32` for SI and the per-site stats only if the golden test tolerates it; it probably won't, so this is low priority.
- [ ] Longer term (gemini 3.6): stream per chromosome through raw stats, smoothing and writing so that whole-genome `[]BSAstats` and `[]SmoothedStats` never coexist. This is a bigger refactor, and it fits with the in-memory pipeline that the genome-wide null (claude_review B1) needs anyway.

### 1.5 Correctness and safety in the code
- [x] A site with `NSamples != len(idxs)` is skipped silently (`stats.go:253`), and the progress bar is not advanced for it. The filter guarantees the count, so this should be an error (or at least a counted warning). As it stands, a bug would show up only as fewer sites.
- [x] Mode flags in `writeRawTSV` use substring matching (`strings.Contains(bsaType, "hp")` …). Derive them once from the role indices already computed in `RawStats` (or reuse `filter.modeRoles`) and pass them in. Two sources of truth for 10 modes will drift.
- [ ] Delete the commented-out vcfgo code (`stats.go:53-181`) once the new filter has been run on a real VCF. Git history keeps it.
- [x] Parameter names: `GStatistic(highBulkAlt, highBulkRef, …)` is called with `(H, L, …)`, and `lod(ref1, alt1, …)` with `(L, H, …)`. Both are symmetric, so the results are correct, but rename the parameters to `a1, b1, a2, b2` so the next reader is not misled.
- [x] `GStatistic`: the `if hba > 0` guards and the `total == 0` check are dead code, because every cell is ≥ 0.5. Clean this up together with 2.2.

## Part 2: changes to the statistics (results change; compare before and after)

### 2.1 Rounding to 6 decimal places (MEDIUM) — DONE
Every stored value is rounded with `math.Round(x*1e6)/1e6` (`stats.go:287-333`).
- ED⁴ = Δ⁴ becomes 0 whenever |Δ| < 0.027, which is about 14% of null sites. That flattens the low end of the ED⁴ distribution before smoothing and robust Z, so the median and MAD of ED⁴ are distorted.
- ED⁴ is computed from SIs that were already rounded.
- [x] Store full precision and round only when writing the TSV. Write ED⁴ (and BBLogBF/G at extreme values) with `%g` or `'g', 7`, not `%.6f`. This also removes about 15 `math.Round` calls per site.

- [x] AD columns (from the audit): `HighBulkAD`/`LowBulkAD` stay REF,ALT; new `HighBulkOther`/`LowBulkOther` columns hold reads on other ALTs, so AD + Other accounts for every read.

### 2.2 G-statistic pseudo-counts (LOW–MEDIUM; gemini 2.6, claude A5) — DONE
`GStatistic` adds 0.5 to all four cells. This is neither Magwene et al. (2011) G nor a Haldane correction (0.5 added only to zero cells). It pulls G toward 0 at low depth.
- [x] Now the standard G with 0·ln 0 = 0 (no pseudo-count), which matches Magwene and QTLseqr and lets G′ methods be applied later. Keep the current version behind a flag only if someone needs it.

### 2.3 LOD duplicates G (MEDIUM; claude A5, gemini, grok)
In two-bulk mode `lod` = G_nopseudo / (2 ln 10); one-bulk LOD1 = G1 / (2 ln 10) exactly. After robust Z they are the same track, so CompositeZ counts that statistic twice.
- [ ] Remove LOD from CompositeZ (k = 2), in `smoothing.go` and in the threshold null together. Note that CompositeZ may be replaced altogether by the genome-wide null (claude_review B1).
- [ ] Either drop the LOD column, or keep it as a convenience column documented as "G / 4.605". It is not a Morton LOD, and the README's "LOD > 3" rule does not apply.
- [ ] Once 2.2 is done, LOD becomes an exact rescaling of G, so it should be derived from G (one less log-likelihood computation per site) if it stays.

### 2.4 ED⁴ definition (cosmetic) — DONE (×4)
- [x] Either multiply by 4 (the literature ED⁴ for a biallelic site is 4Δ⁴) or rename the column to `DeltaSI4`. Published cutoffs do not transfer otherwise.
- [ ] Related (grok #4): the literature smooths the distance and then takes the fourth power. The current code smooths Δ⁴ per site, and by Jensen's inequality that keeps noise. This is a smoothing change, but it means raw stats should keep counts or SIs, not just per-site ED⁴. See 2.8.

### 2.5 2b polarity: ΔSI sign is arbitrary per SNP (HIGH for 2b; claude A3, gemini) — DONE
- [x] `--ref-parent high|low|none` (default `none`). `low` sets `HighIsAlt` in the filter; `high` keeps REF. Ignored, with a warning, when a parent sample is given.
- [x] With `none`, `2b` stores, smooths and thresholds |ΔSI| (`AbsDeltaSI` columns). The null uses |ΔSI| with P99 = q0.99 and P95 = q0.95 and no lower tail (NaN). The `DeltaSI-` track is dropped, CompositeZ uses Z(|ΔSI|) unfolded (matching its null), and two-bulk BRM is skipped.
- [ ] Remaining: smoothed |ΔSI| has a depth-dependent baseline. Per-depth single-SNP thresholds partly cover it; the genome-wide smoothed null (B1) is the full fix.

### 2.6 Depth = the smaller bulk depth (MEDIUM; gemini 2.5, copilotfusion #2) — DONE
`s.Depth = min(highTotal, lowTotal)` (`stats.go:306`) is used for smoothing weights, BRM and threshold lookup. With unequal libraries (e.g. 80x / 25x) the simulated null uses 25x/25x and overstates sequencing noise.
- [x] Store `HighDepth` and `LowDepth` separately. `Depth` is now the harmonic mean 2·dH·dL/(dH+dL) (equal depth with the same ΔSI read variance), used for smoothing weights and the CompositeZ null's median depth. BRM uses both depths. Keep `Depth` only as a derived value for smoothing weights. Choose the weight deliberately: inverse variance ≈ 1/(1/dH + 1/dL) is better than min.
- [x] Pass both depths into the threshold cache key and the simulation (`thresholds.go:892-893`).

### 2.7 One-bulk p₀ (MEDIUM for one-bulk modes) — DONE
`OneBulkP0` is the theoretical `ExpectedAF(population)` (0.5 for F2). Reference mapping bias (typically 0.52–0.55 REF) shifts every site, so AFDev, G1, LOD1 and BBLogBF1 carry a genome-wide offset that depends on `HighIsAlt`.
- [x] `--one-bulk-p0 empirical|theoretical` (default empirical). Estimates the REF mapping bias b as half the difference between the interquartile means of the high-allele fraction at REF-high and ALT-high sites (bulk with ≥2 reads of each allele), then p₀ = theory ± b by polarity. A QTL or wrong population setting shifts both classes equally and cancels. Falls back to theory with < 200 sites per class. Thresholds still simulate around the theoretical p₀ (the AFDev/G1 null shape changes negligibly for |b| of a few %).
- [x] `OneBulkP0` is a constant but is stored and written on every row. Write it once in the log or header, or keep it only if it becomes per-site.

### 2.8 What raw stats should hand to smoothing (design; grok #4/#5, claude B1)
Smoothing per-site G, LOD, ED⁴ and BF is not the same as computing them from smoothed frequencies, and the per-site values are read-level (they ignore pool size).
- [ ] Make `BSAstats` carry the counts (H/L per bulk) and the SIs as the primary data. Treat G, ED⁴ and BF as derived values that can be computed either per site (for the TSV) or from smoothed counts or frequencies. This also makes it straightforward to feed simulated null replicates through the same code (B1).
- [ ] Document that BBLogBF and G are read-level scores: Kass–Raftery and χ² cut-offs do not apply, and only the simulated thresholds are meaningful.

---

## Suggested sequence
1. Part 0 (golden test, benchmark, known-value tests).
2. 1.1 (writer and error checks), then 1.3 and 1.4. Measure after each.
3. 1.2 (drop GT/AD from the struct), then 1.5 (cleanup).
4. 2.1 (stop rounding). This is the smallest change to results.
5. 2.6 (separate depths) and 2.5 (2b polarity). These touch thresholds, so coordinate them with B1.
6. 2.2 / 2.3 / 2.4 (G, LOD and ED⁴ definitions), in one release with README changes.
7. 2.7, then 2.8 (alongside the genome-wide null work).
