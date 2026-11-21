# Plotting to-do

Goal: speed, accuracy and cleanliness at 9/10 or better.
Current rating (audit of 2026-11-21): speed 3, accuracy 5, cleanliness 6.

Scope: `plots/plots.go` (`GeneratePlots`, step 9 in `run/run.go`), which writes `plots/GoBSAseq.<bsaType>.individual_plots.html` and `plots/GoBSAseq.<bsaType>.composite_signal.html`.

## Audit measurements

Synthetic `2p2b` data, 10 chromosomes, no BRM blocks:

| Markers | Time | HTML (both pages) |
|---|---|---|
| 100k | 4.6 s | 58 MB |
| 500k | 22.7 s | 295 MB |
| 2M | 92 s | 1.19 GB |

- About 46 µs and 590 bytes per marker, growing linearly. Above a few hundred thousand markers the pages can't be opened in a browser.
- CPU profile: 63% in go-echarts' post-render `regexp.ReplaceAll` over the whole output (it rewrites the `FuncOpts` JS markers), 18% in JSON encoding. GoBSAseq's own code is under 5%. Time follows output size, so cutting output fixes both time and file size.

## 1. Speed

### 1.1 Threshold lines as mark lines
- **Problem:** each p99/p95 line (and valley line) is a full series repeating the threshold once per marker: 3–5 of the 5–7 series on a chart (`createRawValueChart`, `createCompositeChart` `mkRef`).
- **Change:** draw thresholds as `markLine` `yAxis` entries on the statistic series (one element per line). Keep them in the legend with a label giving the value, e.g. `p99 (0.412)`.
- **Acceptance:** output size per chart drops by about 60% before any other change.

### 1.2 Downsampling with peaks kept
- **Problem:** every marker goes into every chart, so output grows with genome size.
- **Change:** reduce each chromosome's track to at most `maxPoints` (default 5000; a `--plot-points` flag) with LTTB, always keeping every point where the track crosses a threshold and each region's peak, so no threshold crossing disappears from the plot. All tracks of a chromosome share the same x positions, so the tooltip still lines up.
- Alternative or addition: ECharts `sampling: 'lttb'` for drawing speed. It does not shrink the file, so it does not replace the Go-side downsampling.
- **Acceptance:** 2M markers in ≤ 5 s and ≤ 30 MB for both pages; size roughly flat from 500k to 5M markers (benchmark, 1.4).

### 1.3 Drawing cost
- Turn off `Smooth: true` (also accuracy 2.7). Set `showSymbol: false` and `animation: false` on large series.
- Build each chromosome's x data once, not once per chart and again in `brmMarkAreaOpts`; use `strconv.AppendInt` instead of `fmt.Sprintf`.

### 1.4 Benchmark
- `BenchmarkGeneratePlots` at 100k, 500k and 2M markers, reporting time and bytes written, so regressions show up.

## 2. Accuracy

### 2.1 Unsigned `2b` SI charts get false thresholds (bug)
- **Problem:** in `2b` without `--ref-parent`, `siCharts` passes 0 for the thresholds, and `createRawValueChart` always draws p99/p95. Dashed lines appear at y = 0, the subtitle says `p99 threshold: 0.0000`, and the tooltip marks every position `★ p99` (a REF fraction is always > 0).
- **Change:** a chart without thresholds draws no threshold lines and has no significance marks. Pass thresholds as an optional list (e.g. `[]threshLine{name, value, colour}`), not as four floats with 0 meaning "none".
- **Acceptance:** test: the unsigned `2b` SI charts contain no threshold series or mark lines.

### 2.2 NaN or ±Inf blanks a chart without an error (bug)
- **Problem:** go-echarts ignores the JSON error and writes `let option_X = ` with nothing after it. The chart doesn't render, and the CLI still prints `Plot written` (confirmed: one NaN in CompositeZ gives a 1.2 kB page with no data). Smoothed values are protected (`robustZScale`), but thresholds can be NaN (one-sided tracks' lower tail) or ±Inf (`extremeQuantile` with no scored markers).
- **Change:** check every value before rendering. Data points that aren't finite become gaps (`"-"` / null). A threshold that isn't finite means the line is not drawn, with a note in the subtitle (`no threshold`). After rendering, check that the page contains no empty `option_` and return an error if it does.
- **Acceptance:** test: NaN data points and ±Inf/NaN thresholds give a valid page with the other values plotted.

### 2.3 HMM LOD threshold shown as 0 when missing (bug)
- **Problem:** `finite()` turns a NaN/Inf HMM LOD threshold into 0. Detection treats that case as "no threshold" (`FindPeakIntersections` returns nil), but the plot draws p99 at 0 and flags everything above it.
- **Change:** use 2.2's rule: no line, subtitle note. Remove `finite()`.

### 2.4 MaxAbsZ compared with CompositeZ thresholds
- **Problem:** the composite tooltip marks MaxAbsZ `★ p99 significant` against the CompositeZ thresholds. MaxAbsZ (max |Z| over all statistics) has a different null distribution and no genome-wide threshold.
- **Change:** mark significance for CompositeZ only. Show MaxAbsZ as a context line, labelled `max |Z| (no threshold)`, or drop it.

### 2.5 Lower tails not marked
- **Problem:** the tooltip only checks `val > t99`. ΔSI, SI and AFDev have lower-tail thresholds, and detection calls lower-tail (`−`) QTLs, but valley crossings are never marked.
- **Change:** mark `val ≤ valley p99/p95` the same way, and use `≥`/`≤` to match detection ("at or beyond the line", `scanCrossings`). Better: generate the tooltip from the chart's threshold list (2.1), so the marks always follow the drawn lines.

### 2.6 Physical x-axis
- **Problem:** positions are passed as string labels, so ECharts uses a category axis: markers are evenly spaced whatever the distance between them. Marker-dense regions are stretched and sparse ones squashed, so QTL widths and BRM shading don't match their bp size or the intervals in the QTL tables.
- **Change:** numeric (`type: 'value'`) x-axis with `[pos, value]` pairs, min 0, max the chromosome's last position. BRM mark areas use bp coordinates directly (no snapping to markers; also fixes 2.9).
- **Acceptance:** a test with uneven marker spacing: the HTML's x data are the positions, and the BRM area ends equal the block's `Start`/`Stop`.

### 2.7 Curve smoothing overshoots
- **Problem:** `Smooth: true` draws Bézier curves through already Gaussian-smoothed data. Near a threshold the curve can cross the line where no data point does.
- **Change:** straight line segments (1.3).

### 2.8 Detected QTLs not drawn
- **Problem:** only BRM blocks are shaded. Individual, CompositeZ and final consensus QTLs are not drawn, so the tables can't be checked against the plots by eye.
- **Change:** pass the detection results into `GeneratePlots`:
  - each statistic's chart: its own regions (p99 solid, p95 hatched or lighter), with direction;
  - composite chart: CompositeZ regions;
  - every chart: the final consensus QTLs' core as a band and the lead position as a vertical mark line; partial-support loci in a lighter colour.
- Tooltip or label on each band: QTL id, support, interval, as in the final table.
- **Acceptance:** the QTL bands drawn in the HTML equal the intervals in `final.tsv` (test on the golden dataset).

### 2.9 BRM block between two markers
- With the category axis, a block lying entirely between two markers is shaded on the next marker. Fixed by 2.6.

### 2.10 Y-axis floor
- Unsigned charts clamp y at 0. Check that HMM LOD can't be negative; if it can, don't clamp it.

## 3. Viewing

- **Genome overview:** add a Manhattan-style page (all chromosomes side by side on one axis, chromosome colours alternating) for CompositeZ, HMM LOD and ΔSI, with the consensus QTLs marked. This is the first thing a user needs; the per-chromosome pages are for zooming in.
- **Offline:** the pages load ECharts from `go-echarts.github.io`, so they fail on machines without internet (HPC nodes, field laptops). Embed `echarts.min.js` and the theme once per page (go:embed), or add `--plot-offline`.
- **Chromosome order:** use natural order (chr1, chr2 … chr10) in the plots. Do it in `stats.byChromosome` too so plots and tables match; one shared helper.
- **Navigation:** a chromosome index at the top of each page with links to its charts.

## 4. Cleanliness

- **One chart builder:** `createRawValueChart` and `createCompositeChart` share most of their setup. One function taking a track spec: title, y label, series, threshold lines (2.1), regions (2.8), signed/unsigned.
- **Generated tooltip:** one JS tooltip built from the chart's threshold list instead of two hand-written ones with different rules.
- **Group by chromosome** with `stats.byChromosome` (export a wrapper), not a separate map-and-sort in `GeneratePlots`.
- **Dead code and comments:**
  - delete `absMax`;
  - move `createRawValueChart`'s doc comment back above it (it is now above `siCharts`);
  - fix the header (it names `plotting.go`, `thresholds.go`, and refers to "legacy" twoBulk.go and hard-coded ±3/±2, which are gone); delete the `zSig`/`zSugg` constants if unused.
- **Split** `plots.go` once 2.8 and 3 are in: page assembly, chart builders, downsampling.

## 5. Tests

The `plots` package has no tests.

- **Unit:**
  - unsigned `2b`: no thresholds on SI charts (2.1);
  - NaN data and NaN/Inf thresholds give a valid page (2.2, 2.3);
  - lower-tail threshold lines and marks present for ΔSI/AFDev (2.5);
  - downsampling keeps every threshold crossing and every peak (1.2);
  - x data are bp positions; BRM and QTL areas at their bp ends (2.6, 2.8);
  - every mode (`2p2b`, `hp2b`, `lp2b`, `2b` with and without `--ref-parent`, one-bulk modes) renders with the expected charts.
- **Page check:** parse each `option_…` object as JSON; fail on an empty one.
- **Golden:** on the QTL golden dataset, the drawn regions equal the QTL tables.
- **Benchmark:** 1.4.

## 6. Docs

- **README:** what each page and chart shows; thresholds (genome-wide, p99/p95, two-tailed tracks at 99.5/97.5 labelled p99/p95); what the bands mean; downsampling and `--plot-points`; offline viewing.

## Order of work

1. Accuracy bugs 2.1–2.5 and the NaN check: small, independent.
2. Threshold mark lines (1.1), numeric x-axis (2.6), no smoothing (1.3/2.7).
3. Downsampling (1.2) and benchmark (1.4).
4. QTL bands (2.8) and the single chart builder (section 4).
5. Genome overview, offline assets, natural chromosome order (section 3).
6. Tests (section 5) alongside each step; docs.

## Progress

- **2026-11-21:** audit done; nothing changed yet.
