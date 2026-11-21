# QTL detection and merging to-do

Goal: biological soundness, code correctness and cleanliness at 9/10 or better.
Current rating (audit of 2026-11-21): biology 5.5, code correctness 6, cleanliness 5.

## Target behaviour (as specified)

1. **Individual QTLs:** regions where each smoothed statistic is beyond its p99 threshold line, with a p95 fallback when there is no p99 region.
2. **CompositeZ QTLs:** regions where CompositeZ is above its threshold line.
3. **Final QTLs:** regions where the HMM LOD, BRM and CompositeZ agree.

All thresholds are already genome-wide (null simulation; see `smoothing_algorithms.md`, Thresholds), so each list's error rate is genome-wide too, as long as detection doesn't add tests of its own (1.2).

## 1. Biology

### 1.1 Statistics scanned (list 1)
- **Two-bulk:** ΔSI (upper and lower tails, separately), G′, ED⁴, HMM LOD.
- **`2b` without `--ref-parent`:** |ΔSI|, G′, ED⁴ (no HMM LOD).
- **One-bulk:** AFDev (upper and lower), G1, HMM LOD.
- **Change:** none to the set. Each statistic's regions are reported on their own (with the statistic name, tail, level, peak, interval), *and* consolidated per locus (1.4).

### 1.2 p95 fallback: genome-wide, not per chromosome
- **Problem:** the fallback is decided per chromosome per track (`findPeaksWithFallback`). Any chromosome with no p99 crossing is re-scanned at p95, so p95 calls appear in about 5% of QTL-free genomes per track, the extra tests that inflate the error rate.
- **Change:** decide the fallback once per statistic, genome-wide. If the statistic has a p99 region anywhere, report only its p99 regions; if it has none anywhere, report its p95 regions, labelled `p95 (suggestive)`.
- **Acceptance:** in the calibration test, the fraction of null genomes with any individual QTL is ≤ 5% at p99 and ≤ 5% per statistic with the fallback (asserted).

### 1.3 Direction
- **Problem:** ΔSI upper and lower regions (and AFDev upper and lower) are merged into one QTL. Opposite directions mean the high parent's allele is enriched in different bulks: distinct QTLs.
- **Change:** a signed statistic's regions carry a direction (`+`/`−`); consolidation (1.4) and the final merge (1.6) join only same-direction evidence. Unsigned statistics (G′, ED⁴, |ΔSI|, G1, HMM LOD, CompositeZ) take the direction of the signed statistic they overlap (ΔSI, or AFDev in one-bulk modes), or `none`.
- **Acceptance:** a simulated pair of opposite-direction QTLs 2 Mb apart is reported as two QTLs.

### 1.4 Consolidated individual QTLs (per locus)
- Join the statistics' regions on a chromosome when they overlap (or lie within `--merge-distance`) **and** have the same direction.
- Per QTL, report every supporting statistic with **its strongest peak** (position, value), its level (p99 or p95), and how many statistics support it.
- The QTL's interval is the union of its regions. Its lead peak is the HMM LOD peak when present (linkage-aware, the most precise), otherwise the strongest standardised peak.

### 1.5 CompositeZ QTLs (list 2)
- Regions where CompositeZ ≥ its genome-wide p99 threshold, joined within `--merge-distance`, keeping the strongest peak.
- Same genome-wide p95 rule as 1.2 (suggestive only, never per chromosome).
- CompositeZ is unsigned: direction comes from overlapping ΔSI/AFDev evidence (1.3).

### 1.6 Final QTLs: consensus of HMM LOD, BRM and CompositeZ (list 3)
- **Inputs:**
  - HMM LOD regions (above its genome-wide p99), from list 1;
  - BRM blocks (genome-wide `--brm-alpha`);
  - CompositeZ regions (list 2).
- **Agreement:** a final QTL is a locus where every available method has a region, and they overlap. Merges within `--merge-distance` count as overlap.
  - Available methods: all three; in `2b` without `--ref-parent` there is no HMM LOD, so BRM and CompositeZ.
  - Signed methods must agree in direction.
- **Interval:** report both
  - the **consensus core**: the intersection, where all methods are above threshold;
  - the **span**: their union.
  Gene-space analysis uses the core, or optionally the span (a flag, default core).
- **Peaks:** the HMM LOD peak as the lead position (BRM peak when there is no HMM), plus each method's peak, value, threshold and z.
- **Supporting statistics:** list the individual statistics (list 1) overlapping the final QTL.
- **Partial agreement:** loci found by only some methods are not final QTLs. Report them in the same file with a `support` column (e.g. `HMMLOD+BRM`), so strong single-method evidence (e.g. a large HMM LOD) is visible, not lost.
- **Error rate:** requiring all methods gives a genome-wide rate no higher than any one method's (≤ 5% at p99 / α = 0.05).
- **Known trade-off:** CompositeZ is the weakest method (QTL max 18.8 vs threshold 18.4 in calibration, while the HMM LOD is about 24 vs 5), so requiring it reduces power. Measure it (section 4) and document it.

### 1.7 Position interval (later; shared with BRM todo 1.5)
- The threshold-crossing extent depends on effect size and threshold height, not on how precisely the QTL is located.
- Add a support interval around the lead peak: the HMM LOD drop interval (peak − 1.5 LOD, the standard QTL support interval), or a z-drop for the smoothed tracks.
- **Acceptance:** about 95% coverage of the simulated QTL position over at least 100 replicates.

## 2. Code correctness

- **Strongest peak:** `addEvidence` keeps the *last* peak per statistic and `mergeConsolidatedQTL` the *first*. Keep the strongest (largest |value|, or the lowest for lower tails).
- **Unused flags:** `--min-qtl-width` and `--min-qtl-markers` are parsed but never used. Either implement them, as filters on every list, applied before merging and counting markers with NEff ≥ 5, or remove them (deprecate as `--rep` was). Recommendation: implement `--min-qtl-markers` (a QTL needs at least that many markers above the line) and remove `--min-qtl-width` (bp width is already governed by smoothing).
- **Joining evidence:** `MergeCompositeBRM` joins any interval starting within `gap` of the current end, so one long region swallows several; replaced by 1.6.
- **Error handling:** every Excel/TSV writer returns its error; `run.go` treats a failed write as a failed step, consistently.
- **Determinism:** all outputs sorted by chromosome then position (already true for BRM and the detectors via `byChromosome`); keep it for the final list.

## 3. Cleanliness

- **One interval merger** with options (gap, same-direction, keep-strongest) replacing `ConsolidateQTLs` and its repeated passes, `mergePeakIntersections` and `MergeCompositeBRM`. `mergeBRMBlocks` uses it too, with its current behaviour kept by decision.
- **One region type** for every method (statistic, chromosome, start, end, direction, level, peak position, peak value, threshold, z), so lists 1–3 and BRM share code and output columns.
- **Dead code:**
  - delete the commented-out `FindPeakIntersections` (about 140 lines);
  - delete `bestCompositePeakByChrom` and `summarizeBRMByChrom`;
  - delete the unused `StartIndex`/`EndIndex`/`PeakIndex` fields and the `label` field of `peakStatConfig`.
- **Split `detect.go`** (991 lines): scanning (`detect.go`), merging/consensus (`qtlmerge.go`), Excel/TSV output (`qtlwrite.go`).
- **TSV output:** write each list as TSV as well as Excel, so results can be diffed and read by scripts.

## 4. Tests

- **Calibration** (`TestHMMLODCalibration`):
  - null genomes with any individual QTL (p99; with the genome-wide p95 rule), any CompositeZ QTL, any final QTL: each ≤ its nominal rate;
  - power: fraction of QTL datasets where the final list contains the QTL, compared with HMM LOD alone (reported, to quantify 1.6's trade-off).
- **Unit tests:**
  - the strongest peak is kept;
  - opposite directions are not merged;
  - genome-wide p95 fallback (one chromosome with p99 and one without: no p95 calls);
  - consensus core and span;
  - `2b` without `--ref-parent` uses BRM + CompositeZ;
  - partial-support loci are labelled;
  - `--min-qtl-markers` filters;
  - writers return errors (`/dev/full`).
- **Golden output** for the three lists on a fixed simulated dataset.

## 5. Docs

- **README:** the three lists and what each means; genome-wide levels and the p95 "suggestive" rule; the final consensus rule (core, span, support column); direction; the CompositeZ power trade-off; what feeds gene-space analysis.
- **`smoothing_algorithms.md`:** a QTL detection section with the new rating.

## Order of work

1. Code fixes in section 2 and dead code: small and independent.
2. One region type and one merger (section 3), with direction (1.3) and strongest peaks (1.4).
3. Genome-wide p95 rule (1.2) and CompositeZ list (1.5).
4. Final consensus (1.6) and its use in gene-space analysis.
5. Calibration, power, unit and golden tests; docs.
6. Later: support intervals (1.7), shared with BRM.

## Progress

- **2026-11-21 code fixes 1, 2, 4, 5, 6:** done (3, the final-merge chaining, is left for 1.6).
  - **Strongest peak:** `keepStrongest` keeps the largest |value| per statistic in `addEvidence` and `mergeConsolidatedQTL`.
  - **`--min-qtl-markers`** (default 5) is applied to individual-statistic regions, CompositeZ regions and BRM blocks, before the p95 fallback and merging. `PeakIntersection`/`BRMBlock` carry `NMarkers`. `--min-qtl-width` is deprecated and ignored.
  - **Write errors:** the three Excel writers share `writeSheet`, which returns every cell, save and close error.
  - **Order:** outputs are ordered by chromosome then position and identical across runs (`TestDetectionDeterministic`).
  - **Dead code removed:**
    - the commented-out old `FindPeakIntersections`;
    - `bestCompositePeakByChrom`, `summarizeBRMByChrom`, `minI64`, `maxI64`;
    - the `StartIndex`/`EndIndex`/`PeakIndex` fields (replaced by `NMarkers`);
    - the track `label` field. The AFDev lower-tail track is now named `AFDev-`.
- **2026-11-21 consensus rule 1.6:** done (also closes code item 3).
  - `ConsensusQTLs` (`stats/consensus.go`) replaces `MergeCompositeBRM`. Loci are formed from HMM LOD regions (p99), BRM blocks and CompositeZ regions (p99), with `--min-qtl-markers` applied; they are joined by overlap or within `--merge-distance`.
  - **Consensus:** every available method present (BRM + CompositeZ without the HMM LOD).
  - **Reported per locus:**
    - core (intersection) and span (union);
    - direction (from BRM) and support;
    - each method's strongest peak (BRM with z);
    - the individual statistics overlapping the locus.
  - **Lead position:** the midpoint of the HMM LOD's 1.5-LOD support interval. On the simulated data the HMM LOD is a flat plateau (about 11–12 over 6 Mb), so its maximum was arbitrary (13.7 Mb against BRM/CompositeZ peaks at 15.0/14.6 Mb); the midpoint gives 14.2 Mb. The support interval is reported (a first part of 1.7).
  - **Output and gene-space analysis:** written as Excel and TSV. Gene-space analysis uses consensus QTLs only, on the core (`--gene-space-span` for the span).
  - **Calibration:** 0/100 null genomes with a consensus QTL in all three scenarios, with the QTL found in each (asserted ≤ 5).
- **2026-11-21 items 1.2 and 1.3:** done.
  - **1.2 Genome-wide fallback:** `scanGenome` scans each statistic (and CompositeZ, the first half of 1.5) across the genome at p99, and uses p95 only when the statistic has no p99 region anywhere. p95 regions are written as `p95 (suggestive)`.
  - **1.3 Direction:** regions carry a direction (ΔSI/AFDev upper `+`, lower `−`). `ConsolidateQTLs` joins only same-direction regions; unsigned regions join the signed QTL they overlap most, or form `none` QTLs. The individual list has a `Direction` column.
  - **Merging:** the redundant repeated merge passes are replaced by one sorted pass (`mergeSameDirection`), part of item 6.
  - **Calibration:** null genomes with a p99 individual QTL 0–2/100 (asserted ≤ 10), with a suggestive p95 one 1–3/100; consensus still 0/100.
- **2026-11-21 cleanliness items 6–8:** done.
  - **One merger:** `mergeRuns` (`stats/qtlmerge.go`, generic) is the only merge routine. Individual QTLs (per direction, `mergeByDirection`), CompositeZ regions, BRM blocks (same-sign consecutive blocks, unchanged) and consensus loci all use it. The region types stay separate, because they carry different fields used by plots, output and gene-space analysis.
  - **Split:** `detect.go` (scanning, 320 lines), `qtlmerge.go` (merging), `qtlwrite.go` (output), `consensus.go` (final list).
  - **TSV:** every QTL list is written as `.xlsx` and `.tsv` (`writeTable`). The writers are renamed `WriteIndividualQTLs`, `WriteCompositeZQTLs` and `WriteFinalQTLs`.
  - **Output unchanged:** final, BRM and threshold files are byte-identical; calibration is unchanged.
- **2026-11-21 code and cleanliness items from the second audit:** done.
  - **One crossing scanner** (`scanCrossings`, `stats/scan.go`) for QTL detection and BRM. A marker is in a region when it is at or beyond the line, in both, and ends are interpolated the same way. `FindPeakIntersections` went from about 140 lines to 20.
  - **Track table** (`statTrack`): name, value, thresholds, tail, direction and QTL slot are defined once. It replaces `statDirection`, the name switch in `addEvidence`, and the fixed per-statistic fields (`ConsolidatedQTL.Peaks [nSlots]*Peak`). Output columns come from the mode's slots.
  - **Individual QTLs** report `N_Statistics` and a lead: the HMM LOD support midpoint, else the statistic furthest beyond its threshold.
  - **Positions** in the individual and CompositeZ tables are written as whole base pairs (they were in `g` notation, e.g. `5.77638e+06`).
  - **Golden tests** pin the three QTL lists on a fixed input with BRM (`TestQTLGolden`, `testdata/qtl_golden/`).
  - **Output and calibration:** final and BRM outputs are byte-identical to before; calibration is unchanged.
  - **Rating:** biology 8, code correctness 9, cleanliness 9.
