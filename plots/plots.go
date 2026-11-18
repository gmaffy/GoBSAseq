package plots

// plotting.go — HTML chart generation for GoBSAseq.
//
// Adapted from the legacy twoBulk.go plotting code.  All chart functions now
// operate on the new package-level types:
//
//   • SmoothedStats   (smoothing.go)  — Gaussian-kernel smoothed values and
//                                       robust Z-scores produced by
//                                       SmoothAndNormalise.
//   • Thresholds      (thresholds.go) — per-variant empirical significance
//                                       thresholds produced by
//                                       CalculateThresholds.
//   • BRMBlock        (brm.go)        — BRM block intervals produced by
//                                       RunBRM.
//
// Three HTML pages are written per bsaType:
//
//   GoBSAseq.<bsaType>.individual_plots.html   — per-statistic raw-value charts
//   GoBSAseq.<bsaType>.robust_z_overlay.html   — all Z-scores overlaid per chrom
//   GoBSAseq.<bsaType>.composite_signal.html   — CompositeZ / MaxAbsZ chart
//
// The chart style (theme, colours, data-zoom, toolbox, BRM mark-areas) is
// preserved exactly from the legacy code.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
)

// ── chart-wide style constants ───────────────────────────────────────────────

const (
	chartTheme  = types.ThemeWesteros
	chartWidth  = "900px"
	chartHeight = "500px"

	// Z-score reference thresholds (identical to legacy values).
	zSig  = 3.0 // ~p99 equivalent — significant
	zSugg = 2.0 // ~p95 equivalent — suggestive
)

// posFormatter formats a genomic position (integer) as Mb / kb / bp.
const posFormatter = `function(value) {
	if (value >= 1000000) { return (value / 1000000).toFixed(2) + ' Mb'; }
	if (value >= 1000)    { return (value / 1000).toFixed(1)    + ' kb'; }
	return value;
}`

// ── public entry point ───────────────────────────────────────────────────────

func GeneratePlots(cfg utils.AnalysisConfig, bsaType string, smoothed []stats.SmoothedStats, th stats.Thresholds, brmBlocks []stats.BRMBlock, qtls []stats.MergedQTL) error {
	color.Cyan("\n================================ Generating plots for %s... ==========================================", bsaType)
	if len(smoothed) == 0 {
		return nil
	}

	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := stats.BulkFlags(bsaType)

	// ── group by chromosome, preserve sort order ─────────────────────────────

	// Group marker indices by chromosome (indices, not copies: each
	// SmoothedStats + Thresholds pair is several hundred bytes).
	byChr := make(map[string][]int)
	for i, s := range smoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], i)
	}

	// Display resolution: curves are drawn at ~1% of the smoothing σ, keeping
	// the minimum and maximum of every bin (see downsampleXY), so every peak
	// and dip survives while files stay small. TSV/Excel outputs keep every
	// marker.
	bw := math.Max(1000, float64(cfg.WindowSize)/2/100)
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// Index BRM blocks by chrom for fast lookup.
	shadeByChr := make(map[string]shading)
	for _, b := range brmBlocks {
		s := shadeByChr[b.Chrom]
		s.brm = append(s.brm, b)
		shadeByChr[b.Chrom] = s
	}
	for _, q := range qtls {
		if q.Support.Valid() {
			s := shadeByChr[q.Chrom]
			s.qtl = append(s.qtl, q.Support)
			shadeByChr[q.Chrom] = s
		}
	}

	// ── create page containers ───────────────────────────────────────────────

	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics"

	robustZPage := components.NewPage()
	robustZPage.SetLayout(components.PageFlexLayout)
	robustZPage.PageTitle = "GoBSAseq — Robust Z-score Overlay"

	compositePage := components.NewPage()
	compositePage.SetLayout(components.PageFlexLayout)
	compositePage.PageTitle = "GoBSAseq — Composite Signal"

	// ── per-chromosome chart generation ─────────────────────────────────────

	for _, chrom := range chroms {
		idxs := byChr[chrom]
		if len(idxs) == 0 {
			continue
		}

		n := len(idxs)
		chromShade := shadeByChr[chrom]

		// Extract position arrays and per-variant threshold arrays.
		x := make([]int64, n)

		// Two-bulk smoothed values.
		hi := make([]float64, n)
		li := make([]float64, n)
		dsi := make([]float64, n)
		gs := make([]float64, n)
		ed := make([]float64, n)
		lod := make([]float64, n)
		bbl := make([]float64, n)

		// One-bulk smoothed values.
		afDev := make([]float64, n)
		obG := make([]float64, n)
		obLod := make([]float64, n)
		obBB := make([]float64, n)

		// Z-scores (all available for every bsaType).
		hiZ := make([]float64, n)
		liZ := make([]float64, n)
		dsiZ := make([]float64, n)
		gsZ := make([]float64, n)
		edZ := make([]float64, n)
		lodZ := make([]float64, n)
		bblZ := make([]float64, n)
		afDevZ := make([]float64, n)
		obGZ := make([]float64, n)
		obLodZ := make([]float64, n)
		obBBZ := make([]float64, n)
		compositeZ := make([]float64, n)
		maxAbsZ := make([]float64, n)

		// Per-variant threshold arrays for two-bulk individual charts.
		hiT99, hiTM99 := make([]float64, n), make([]float64, n)
		hiT95, hiTM95 := make([]float64, n), make([]float64, n)
		liT99, liTM99 := make([]float64, n), make([]float64, n)
		liT95, liTM95 := make([]float64, n), make([]float64, n)
		dsiT99, dsiTM99 := make([]float64, n), make([]float64, n)
		dsiT95, dsiTM95 := make([]float64, n), make([]float64, n)
		gsT99, gsT95 := make([]float64, n), make([]float64, n)
		edT99, edT95 := make([]float64, n), make([]float64, n)
		lodT99, lodT95 := make([]float64, n), make([]float64, n)
		bblT99, bblT95 := make([]float64, n), make([]float64, n)

		// Per-variant threshold arrays for one-bulk individual charts.
		afDevT99, afDevTM99 := make([]float64, n), make([]float64, n)
		afDevT95, afDevTM95 := make([]float64, n), make([]float64, n)
		obGT99, obGT95 := make([]float64, n), make([]float64, n)
		obLodT99, obLodT95 := make([]float64, n), make([]float64, n)
		obBBT99, obBBT95 := make([]float64, n), make([]float64, n)

		for i, j := range idxs {
			s, t := &smoothed[j], &th
			x[i] = s.POS

			if hasBothBulks {
				hi[i], li[i], dsi[i] = s.SmHighSI, s.SmLowSI, s.SmDeltaSI
				gs[i], ed[i], lod[i], bbl[i] = s.SmGstat, s.SmED, s.SmLOD, s.SmBBLogBF
				hiZ[i], liZ[i], dsiZ[i] = s.ZHighSI, s.ZLowSI, s.ZDeltaSI
				gsZ[i], edZ[i], lodZ[i], bblZ[i] = s.ZGstat, s.ZED, s.ZLOD, s.ZBBLogBF

				tb := t.TwoBulk
				hiT99[i], hiTM99[i] = tb.HighSIP99, tb.HighSIMp99
				hiT95[i], hiTM95[i] = tb.HighSIP95, tb.HighSIMp95
				liT99[i], liTM99[i] = tb.LowSIP99, tb.LowSIMp99
				liT95[i], liTM95[i] = tb.LowSIP95, tb.LowSIMp95
				dsiT99[i], dsiTM99[i] = tb.DeltaSIP99, tb.DeltaSIMp99
				dsiT95[i], dsiTM95[i] = tb.DeltaSIP95, tb.DeltaSIMp95
				gsT99[i], gsT95[i] = tb.GstatP99, tb.GstatP95
				edT99[i], edT95[i] = tb.ED4P99, tb.ED4P95
				lodT99[i], lodT95[i] = tb.LODP99, tb.LODP95
				bblT99[i], bblT95[i] = tb.BBLogBFP99, tb.BBLogBFP95
			}

			if hasHighBulk && !hasBothBulks {
				// High-bulk-only: HighSI Z-score available.
				hi[i] = s.SmHighSI
				hiZ[i] = s.ZHighSI
				hiT99[i], hiTM99[i] = t.TwoBulk.HighSIP99, t.TwoBulk.HighSIMp99
				hiT95[i], hiTM95[i] = t.TwoBulk.HighSIP95, t.TwoBulk.HighSIMp95
			}
			if hasLowBulk && !hasBothBulks {
				li[i] = s.SmLowSI
				liZ[i] = s.ZLowSI
				liT99[i], liTM99[i] = t.TwoBulk.LowSIP99, t.TwoBulk.LowSIMp99
				liT95[i], liTM95[i] = t.TwoBulk.LowSIP95, t.TwoBulk.LowSIMp95
			}

			if hasOneBulk {
				afDev[i] = s.SmAFDev
				obG[i] = s.SmOneBulkG
				obLod[i] = s.SmOneBulkLOD
				obBB[i] = s.SmOneBulkBBLogBF
				afDevZ[i] = s.ZAFDev
				obGZ[i] = s.ZOneBulkG
				obLodZ[i] = s.ZOneBulkLOD
				obBBZ[i] = s.ZOneBulkBBLogBF

				ob := t.OneBulk
				afDevT99[i], afDevTM99[i] = ob.AFDevP99, ob.AFDevMp99
				afDevT95[i], afDevTM95[i] = ob.AFDevP95, ob.AFDevMp95
				obGT99[i], obGT95[i] = ob.OneBulkGstatP99, ob.OneBulkGstatP95
				obLodT99[i], obLodT95[i] = ob.OneBulkLODP99, ob.OneBulkLODP95
				obBBT99[i], obBBT95[i] = ob.OneBulkBBLogBFP99, ob.OneBulkBBLogBFP95
			}

			compositeZ[i] = s.CompositeZ
			maxAbsZ[i] = s.MaxAbsZ
		}

		// ── average per-variant thresholds for the reference lines ───────────
		// The individual raw-value charts draw flat average threshold lines
		// (as in the legacy code) to give the reader a quick visual anchor.
		// The adaptive per-variant thresholds are still used for QTL detection;
		// here we only need averages for drawing.
		nf := float64(n)
		avg := func(arr []float64) float64 {
			var s float64
			for _, v := range arr {
				s += v
			}
			return s / nf
		}

		// ── individual raw-value charts ──────────────────────────────────────

		if hasBothBulks {
			individualPage.AddCharts(
				createRawValueChart(bw, chrom+" HighSI", x, hi,
					avg(hiT99), avg(hiT95), avg(hiTM99), avg(hiTM95), true, chromShade),
				createRawValueChart(bw, chrom+" LowSI", x, li,
					avg(liT99), avg(liT95), avg(liTM99), avg(liTM95), true, chromShade),
				createRawValueChart(bw, chrom+" ΔSI (DeltaSI)", x, dsi,
					avg(dsiT99), avg(dsiT95), avg(dsiTM99), avg(dsiTM95), true, chromShade),
				createRawValueChart(bw, chrom+" G-statistic", x, gs,
					avg(gsT99), avg(gsT95), 0, 0, false, chromShade),
				createRawValueChart(bw, chrom+" ED⁴", x, ed,
					avg(edT99), avg(edT95), 0, 0, false, chromShade),
				createRawValueChart(bw, chrom+" LOD", x, lod,
					avg(lodT99), avg(lodT95), 0, 0, false, chromShade),
				createRawValueChart(bw, chrom+" BB log-BF", x, bbl,
					avg(bblT99), avg(bblT95), 0, 0, false, chromShade),
			)
		}

		if hasHighBulk && !hasBothBulks {
			individualPage.AddCharts(
				createRawValueChart(bw, chrom+" HighSI", x, hi,
					avg(hiT99), avg(hiT95), avg(hiTM99), avg(hiTM95), true, chromShade),
			)
		}
		if hasLowBulk && !hasBothBulks {
			individualPage.AddCharts(
				createRawValueChart(bw, chrom+" LowSI", x, li,
					avg(liT99), avg(liT95), avg(liTM99), avg(liTM95), true, chromShade),
			)
		}
		if hasOneBulk {
			individualPage.AddCharts(
				createRawValueChart(bw, chrom+" AF Deviation", x, afDev,
					avg(afDevT99), avg(afDevT95), avg(afDevTM99), avg(afDevTM95), true, chromShade),
				createRawValueChart(bw, chrom+" One-bulk G-statistic", x, obG,
					avg(obGT99), avg(obGT95), 0, 0, false, chromShade),
				createRawValueChart(bw, chrom+" One-bulk LOD", x, obLod,
					avg(obLodT99), avg(obLodT95), 0, 0, false, chromShade),
				createRawValueChart(bw, chrom+" One-bulk BB log-BF", x, obBB,
					avg(obBBT99), avg(obBBT95), 0, 0, false, chromShade),
			)
		}

		// ── robust Z-score overlay ───────────────────────────────────────────

		robustZPage.AddCharts(
			createZOverlayChart(
				bw, chrom, bsaType, x,
				hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ,
				afDevZ, obGZ, obLodZ, obBBZ,
				hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk,
				chromShade,
			),
		)

		// ── composite signal chart ───────────────────────────────────────────
		// Use the empirical CompositeZ thresholds from the first variant on this
		// chromosome (they are global — identical for every variant — so any entry
		// gives the same value).
		czThresh := th.Z

		compositePage.AddCharts(
			createCompositeChart(bw, chrom, x, compositeZ, maxAbsZ, chromShade, czThresh),
		)
	}

	// ── write HTML pages ─────────────────────────────────────────────────────

	outDir := filepath.Join(cfg.OutputDir, "plots")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return fmt.Errorf("GeneratePlots: mkdir %s: %w", outDir, err)
	}

	pages := []struct {
		page *components.Page
		name string
	}{
		{individualPage, fmt.Sprintf("GoBSAseq.%s.individual_plots.html", bsaType)},
		{robustZPage, fmt.Sprintf("GoBSAseq.%s.robust_z_overlay.html", bsaType)},
		{compositePage, fmt.Sprintf("GoBSAseq.%s.composite_signal.html", bsaType)},
	}
	for _, p := range pages {
		path := filepath.Join(outDir, p.name)
		if err := writeHTMLPage(p.page, path); err != nil {
			return err
		}
		color.Green("Plot written to %s", path)
	}
	return nil
}

// ── shared chart options ─────────────────────────────────────────────────────

// commonGlobalOpts returns the standard set of go-echarts global options used
// by every chart in GoBSAseq.  bidirectional = true leaves the Y-axis minimum
// unset (allowing negative values); false clamps to 0.
func commonGlobalOpts(title, subtitle, yLabel, width, height string, bidirectional bool) []charts.GlobalOpts {
	yMin := opts.Float(0.0)
	if bidirectional {
		yMin = nil
	}
	return []charts.GlobalOpts{
		charts.WithInitializationOpts(opts.Initialization{
			Theme:  chartTheme,
			Width:  width,
			Height: height,
		}),
		charts.WithTitleOpts(opts.Title{
			Title:    title,
			Subtitle: subtitle,
			Left:     "center",
			Top:      "1%",
		}),
		charts.WithXAxisOpts(opts.XAxis{
			Type:         "value",
			Min:          "dataMin",
			Max:          "dataMax",
			Name:         "Genomic Position",
			NameLocation: "middle",
			NameGap:      35,
			AxisLabel: &opts.AxisLabel{
				Rotate:    30,
				Formatter: opts.FuncOpts(posFormatter),
			},
		}),
		charts.WithYAxisOpts(opts.YAxis{
			Name:         yLabel,
			NameLocation: "middle",
			NameGap:      55,
			Min:          yMin,
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
		}),
		charts.WithDataZoomOpts(
			opts.DataZoom{Type: "slider", XAxisIndex: []int{0}, Start: 0, End: 100},
			opts.DataZoom{Type: "inside", XAxisIndex: []int{0}},
		),
		charts.WithLegendOpts(opts.Legend{
			Show:   opts.Bool(true),
			Top:    "9%",
			Left:   "center",
			Type:   "scroll",
			Orient: "horizontal",
		}),
		charts.WithToolboxOpts(opts.Toolbox{
			Show:  opts.Bool(true),
			Right: "2%",
			Feature: &opts.ToolBoxFeature{
				SaveAsImage: &opts.ToolBoxFeatureSaveAsImage{Show: opts.Bool(true), Title: "Save PNG"},
				DataZoom:    &opts.ToolBoxFeatureDataZoom{Show: opts.Bool(true), Title: map[string]string{"zoom": "Zoom", "back": "Reset"}},
				Restore:     &opts.ToolBoxFeatureRestore{Show: opts.Bool(true), Title: "Reset"},
			},
		}),
		charts.WithGridOpts(opts.Grid{
			Left:         "8%",
			Right:        "4%",
			Top:          "20%",
			Bottom:       "14%",
			ContainLabel: opts.Bool(true),
		}),
	}
}

// ── individual raw-value chart ───────────────────────────────────────────────

// createRawValueChart builds a single line chart showing one smoothed
// statistic against its empirical p99/p95 threshold lines and BRM block
// mark-areas.  This replaces the legacy createInteractiveLineChart.
//
// t99 / t95     — upper-tail average thresholds.
// tm99 / tm95   — lower-tail (negative) thresholds; ignored when
//
//	hasNegativeThresh is false.
func createRawValueChart(
	bw float64,
	title string,
	x []int64,
	y []float64,
	t99, t95, tm99, tm95 float64,
	hasNegativeThresh bool,
	shade shading,
) *charts.Line {

	subtitle := fmt.Sprintf("p99 threshold: %.4f  |  p95 threshold: %.4f  |  amber: BRM blocks · red: QTL support interval", t99, t95)

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Value", chartWidth, chartHeight, hasNegativeThresh)...)

	line.SetGlobalOptions(charts.WithTooltipOpts(opts.Tooltip{
		Show:        opts.Bool(true),
		Trigger:     "axis",
		AxisPointer: &opts.AxisPointer{Type: "cross"},
		// Threshold values are embedded: the threshold lines are two-point
		// series, so they are not sampled at the hovered position.
		Formatter: opts.FuncOpts(fmt.Sprintf(`function(params) {
			let pos = params[0].axisValue;
			let posFmt = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
			let result = '<strong>Position: ' + posFmt + '</strong><br/>';
			let t99val = %g, t95val = %g;
			params.forEach(function(item) {
				let val = parseFloat(item.value[1]);
				if (isNaN(val)) return;
				let sig = '';
				if (item.seriesName === 'Statistic') {
					if (val > t99val)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ p99</span>';
					else if (val > t95val) sig = ' <span style="color:#f39c12">● p95</span>';
				}
				result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(5) + sig + '<br/>';
			});
			return result;
		}`, t99, t95)),
	}))

	yData := downsampleXY(x, y, bw)
	y99 := constLine(x, t99)
	y95 := constLine(x, t95)
	var ym99, ym95 []opts.LineData
	if hasNegativeThresh {
		ym99 = constLine(x, tm99)
		ym95 = constLine(x, tm95)
	}

	statSeriesOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#1f77b4"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	statSeriesOpts = append(statSeriesOpts, shadingOpts(shade, x)...)

	line.AddSeries("Statistic", yData, statSeriesOpts...).
		AddSeries("p99", y99,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("p95", y95,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	if hasNegativeThresh {
		line.
			AddSeries("p99 valley", ym99,
				charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
				charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
			).
			AddSeries("p95 valley", ym95,
				charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
				charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
			)
	}

	return line
}

// ── robust Z-score overlay chart ─────────────────────────────────────────────

// statColor returns the hex colour assigned to each statistic series, matching
// the legacy palette exactly.
func statColor(name string) string {
	switch name {
	case "HighSI":
		return "#1f77b4"
	case "LowSI":
		return "#ff7f0e"
	case "DeltaSI":
		return "#2ca02c"
	case "Gstat":
		return "#17becf"
	case "ED":
		return "#d62728"
	case "LOD":
		return "#9467bd"
	case "BBLogBF":
		return "#8c564b"
	case "AFDev":
		return "#e377c2"
	case "OneBulkG":
		return "#17becf"
	case "OneBulkLOD":
		return "#9467bd"
	case "OneBulkBB":
		return "#8c564b"
	}
	return "#7f7f7f"
}

// createZOverlayChart draws all available robust Z-score series on one chart
// per chromosome, with ±2 / ±3 reference lines and BRM block shading.
// This replaces the legacy createRobustZOverlayChart; it is bsaType-aware so
// it only draws the series relevant to the current analysis mode.
func createZOverlayChart(
	bw float64,
	chrom, bsaType string,
	x []int64,
	hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ []float64,
	afDevZ, obGZ, obLodZ, obBBZ []float64,
	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool,
	shade shading,
) *charts.Line {

	title := chrom + " — Robust Z-score Overlay"
	subtitle := "Z vs simulated QTL-free genomes (median+MAD). z = ±2 / ±3 per-position guides. " +
		"Amber: BRM blocks · red: QTL support interval."

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Robust Z-score", chartWidth, chartHeight, true)...)

	line.SetGlobalOptions(
		charts.WithYAxisOpts(opts.YAxis{
			Name:         "Robust Z-score",
			NameLocation: "middle",
			NameGap:      55,
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
			AxisLabel: &opts.AxisLabel{
				Formatter: opts.FuncOpts(`function(v) {
					let m = {3:'z=3 ★', 2:'z=2 ●', 0:'0', '-2':'z=-2 ●', '-3':'z=-3 ★'};
					let k = parseFloat(v.toFixed(1));
					return m[k] !== undefined ? m[k] : v.toFixed(1);
				}`),
			},
		}),
		charts.WithTooltipOpts(opts.Tooltip{
			Show:        opts.Bool(true),
			Trigger:     "axis",
			AxisPointer: &opts.AxisPointer{Type: "cross"},
			Formatter: opts.FuncOpts(`function(params) {
				let pos = params[0].axisValue;
				let posStr = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
				let result = '<strong>' + posStr + '</strong><br/>';
				let refSeries = ['z=0','z=+2 (sugg.)','z=+3 (sig.)','z=-2 (sugg.)','z=-3 (sig.)'];
				params.forEach(function(item) {
					if (refSeries.indexOf(item.seriesName) !== -1) return;
					let val = parseFloat(item.value[1]);
					if (isNaN(val)) return;
					let sig = '';
					if (Math.abs(val) >= 3.0)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ significant</span>';
					else if (Math.abs(val) >= 2.0) sig = ' <span style="color:#f39c12">● suggestive</span>';
					result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
				});
				return result;
			}`),
		}),
	)

	mkRef := func(val float64) []opts.LineData { return constLine(x, val) }

	// Reference lines — z=0 carries the BRM mark-area.
	zeroOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#bdc3c7"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	zeroOpts = append(zeroOpts, shadingOpts(shade, x)...)

	line.AddSeries("z=0", mkRef(0), zeroOpts...).
		AddSeries("z=+2 (sugg.)", mkRef(zSugg),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=+3 (sig.)", mkRef(zSig),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=-2 (sugg.)", mkRef(-zSugg),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=-3 (sig.)", mkRef(-zSig),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	// Statistic Z-score series — only rendered when available for this bsaType.
	type zSeriesDef struct {
		name string
		data []float64
		w    float32
	}
	var series []zSeriesDef
	if hasBothBulks {
		series = append(series,
			zSeriesDef{"HighSI", hiZ, 2.0},
			zSeriesDef{"LowSI", liZ, 2.0},
			zSeriesDef{"DeltaSI", dsiZ, 3.0},
			zSeriesDef{"Gstat", gsZ, 2.0},
			zSeriesDef{"ED", edZ, 2.0},
			zSeriesDef{"LOD", lodZ, 2.0},
			zSeriesDef{"BBLogBF", bblZ, 2.0},
		)
	} else {
		if hasHighBulk {
			series = append(series, zSeriesDef{"HighSI", hiZ, 2.0})
		}
		if hasLowBulk {
			series = append(series, zSeriesDef{"LowSI", liZ, 2.0})
		}
	}
	if hasOneBulk {
		series = append(series,
			zSeriesDef{"AFDev", afDevZ, 3.0},
			zSeriesDef{"OneBulkG", obGZ, 2.0},
			zSeriesDef{"OneBulkLOD", obLodZ, 2.0},
			zSeriesDef{"OneBulkBB", obBBZ, 2.0},
		)
	}

	for _, s := range series {
		col := statColor(s.name)
		line.AddSeries(s.name, downsampleXY(x, s.data, bw),
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: s.w, Color: col}),
			charts.WithItemStyleOpts(opts.ItemStyle{Color: col, Opacity: opts.Float(0)}),
		)
	}

	return line
}

// ── composite signal chart ────────────────────────────────────────────────────

// createCompositeChart renders both the Stouffer CompositeZ (signed, from
// smoothing.go) and MaxAbsZ (unsigned envelope) on the same chart.
// Threshold lines are drawn from the empirical MC-derived ZThresholds (P99 / P95)
// rather than the legacy hard-coded ±3 / ±2 constants.
// BRM block shading is included.
func createCompositeChart(
	bw float64,
	chrom string,
	x []int64,
	compositeZ, maxAbsZ []float64,
	shade shading,
	czThresh stats.ZThresholds,
) *charts.Line {

	czP99 := czThresh.CompositeZP99
	czP95 := czThresh.CompositeZP95

	title := chrom + " — Composite Signal"
	subtitle := fmt.Sprintf(
		"CompositeZ (ΔSI + G + LOD, k=%d) and max |Z|. Genome-wide p99 %.2f · p95 %.2f. "+
			"Amber: BRM blocks · red: QTL support interval.",
		czThresh.NumStats, czP99, czP95,
	)

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Z-score", chartWidth, chartHeight, true)...)

	// Embed the empirical threshold values into the JS tooltip so the reader
	// sees "★ p99 significant" / "● p95 suggestive" relative to the actual MC
	// upper-tail thresholds, not the old hard-coded 3.0 / 2.0.
	tooltipJS := fmt.Sprintf(`function(params) {
		var czP99 = %f, czP95 = %f;
		let pos = params[0].axisValue;
		let posStr = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
		let result = '<strong>' + posStr + '</strong><br/>';
		params.forEach(function(item) {
			let val = parseFloat(item.value[1]);
			if (isNaN(val)) return;
			let sig = '';
			if (item.seriesName === 'CompositeZ' || item.seriesName === 'MaxAbsZ') {
				if (val >= czP99)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ p99 significant</span>';
				else if (val >= czP95) sig = ' <span style="color:#f39c12">● p95 suggestive</span>';
			}
			result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
		});
		return result;
	}`, czP99, czP95)

	line.SetGlobalOptions(
		charts.WithTooltipOpts(opts.Tooltip{
			Show:        opts.Bool(true),
			Trigger:     "axis",
			AxisPointer: &opts.AxisPointer{Type: "cross"},
			Formatter:   opts.FuncOpts(tooltipJS),
		}),
	)

	mkRef := func(val float64) []opts.LineData { return constLine(x, val) }

	compositeOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#2ca02c"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	compositeOpts = append(compositeOpts, shadingOpts(shade, x)...)

	p99Label := fmt.Sprintf("p99 (%.3f)", czP99)
	p95Label := fmt.Sprintf("p95 (%.3f)", czP95)

	line.AddSeries("CompositeZ", downsampleXY(x, compositeZ, bw), compositeOpts...).
		AddSeries("MaxAbsZ", downsampleXY(x, maxAbsZ, bw),
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: 1.8, Color: "#9467bd", Type: "dotted"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries(p95Label, mkRef(czP95),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries(p99Label, mkRef(czP99),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	return line
}

// ── BRM mark-area helpers ─────────────────────────────────────────────────────

// shading holds the intervals shaded on a chromosome's charts.
type shading struct {
	brm []stats.BRMBlock        // amber: regions significant by BRM
	qtl []stats.SupportInterval // red: LOD-drop support interval (where the QTL is)
}

// shadingOpts returns go-echarts SeriesOpts that shade BRM blocks (amber)
// and QTL support intervals (red, drawn on top). Safe with empty shading.
func shadingOpts(sh shading, x []int64) []charts.SeriesOpts {
	if len(x) == 0 {
		return nil
	}

	clip := func(start, stop int64) (int64, int64, bool) {
		startIdx := sort.Search(len(x), func(i int) bool { return x[i] >= start })
		stopIdx := sort.Search(len(x), func(i int) bool { return x[i] > stop }) - 1
		if startIdx >= len(x) || stopIdx < 0 {
			return 0, 0, false
		}
		if stopIdx < startIdx {
			stopIdx = startIdx
		}
		return x[startIdx], x[stopIdx], true
	}

	var areas [][]opts.MarkAreaData
	for _, b := range sh.brm {
		if s, e, ok := clip(b.Start, b.Stop); ok {
			areas = append(areas, []opts.MarkAreaData{
				{Name: fmt.Sprintf("BRM block %.4f", b.Peak), XAxis: s},
				{XAxis: e},
			})
		}
	}
	for _, q := range sh.qtl {
		if s, e, ok := clip(q.Start, q.End); ok {
			areas = append(areas, []opts.MarkAreaData{
				{
					Name:  fmt.Sprintf("QTL %.2f–%.2f Mb", float64(q.Start)/1e6, float64(q.End)/1e6),
					XAxis: s,
					MarkAreaStyle: opts.MarkAreaStyle{
						ItemStyle: &opts.ItemStyle{Color: "rgba(192, 57, 43, 0.30)"},
						Label:     &opts.Label{Show: opts.Bool(true), Position: "insideTop", Color: "#922b21"},
					},
				},
				{XAxis: e},
			})
		}
	}
	if len(areas) == 0 {
		return nil
	}

	return []charts.SeriesOpts{
		charts.WithMarkAreaData(areas...),
		charts.WithMarkAreaStyleOpts(opts.MarkAreaStyle{
			Label:     &opts.Label{Show: opts.Bool(false)},
			ItemStyle: &opts.ItemStyle{Color: "rgba(243, 156, 18, 0.22)"},
		}),
	}
}

// ── utility helpers ───────────────────────────────────────────────────────────

// downsampleXY converts a series to [position, value] points for a numeric
// x-axis, keeping only the minimum and maximum point of every bin of width
// bw (in their original order). Every local peak and dip wider than a bin is
// preserved exactly, while the point count drops to at most two per bin.
func downsampleXY(x []int64, y []float64, bw float64) []opts.LineData {
	out := make([]opts.LineData, 0, 2*int(float64(x[len(x)-1]-x[0])/bw+2))
	emit := func(i int) { out = append(out, opts.LineData{Value: [2]float64{float64(x[i]), y[i]}}) }
	for start := 0; start < len(x); {
		bin := int64(float64(x[start]) / bw)
		end, lo, hi := start, start, start
		for end < len(x) && int64(float64(x[end])/bw) == bin {
			if y[end] < y[lo] {
				lo = end
			}
			if y[end] > y[hi] {
				hi = end
			}
			end++
		}
		first, second := min(lo, hi), max(lo, hi)
		emit(first)
		if second != first {
			emit(second)
		}
		start = end
	}
	return out
}

// constLine is a horizontal line across the chromosome: two points instead of
// one per marker.
func constLine(x []int64, v float64) []opts.LineData {
	return []opts.LineData{
		{Value: [2]float64{float64(x[0]), v}},
		{Value: [2]float64{float64(x[len(x)-1]), v}},
	}
}

// writeHTMLPage renders a go-echarts Page to a file.
func writeHTMLPage(page *components.Page, path string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create %s: %w", path, err)
	}
	if err := page.Render(f); err != nil {
		_ = f.Close()
		return fmt.Errorf("render %s: %w", path, err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("close %s: %w", path, err)
	}
	return nil
}

func absMax(a, b float64) float64 {
	if math.Abs(a) > math.Abs(b) {
		return math.Abs(a)
	}
	return math.Abs(b)
}
