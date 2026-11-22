package plots

// plotting.go — HTML chart generation for GoBSAseq.
//
// Two HTML pages are written per bsaType:
//
//   GoBSAseq.<bsaType>.individual_plots.html   — per-statistic charts
//   GoBSAseq.<bsaType>.dsigconsensusz.html     — DsiGconsensusZ chart
//
// Every chart puts genomic position on a linear base-pair axis (a SNP desert
// takes its true width, and a QTL's width on the chart is its width in bp),
// draws each track as tested (straight segments between points, no spline on
// top of the kernel smoothing), its genome-wide p99/p95 lines (and the lower
// "valley" lines of two-tailed tracks), and shades the BRM blocks. The
// tooltip marks a value beyond either tail's line.
//
// Size. Each chromosome's track is simplified (simplify): only the markers
// needed for the line to pass within 0.1% of the track's range of every
// marker are kept, so peaks, valleys and threshold crossings are drawn as in
// the full track while stretches the line already follows take no points
// (the smoothed tracks are linear between smoothing-grid points wherever
// markers are dense). Every marker stays in the smoothed TSV. Points are
// written as [position, value] pairs rounded to the precision the tooltip
// shows, and each threshold is one mark line rather than one point per
// marker.
//
// The chart style (theme, colours, data-zoom, toolbox, BRM mark-areas) is
// kept from the legacy code.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"

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

	// Plotted values are rounded to the decimals their tooltip shows.
	statDecimals      = 5
	consensusDecimals = 3

	colorP99 = "#e74c3c"
	colorP95 = "#f39c12"
)

// posFormatter formats a genomic position (integer) as Mb / kb / bp.
const posFormatter = `function(value) {
	if (value >= 1000000) { return (value / 1000000).toFixed(2) + ' Mb'; }
	if (value >= 1000)    { return (value / 1000).toFixed(1)    + ' kb'; }
	return value;
}`

// thresholds are a chart's genome-wide lines: upper p99/p95 and lower
// ("valley") p99/p95; NaN = no line.
type thresholds struct{ p99, p95, m99, m95 float64 }

func noThresholds() thresholds { return thresholds{math.NaN(), math.NaN(), math.NaN(), math.NaN()} }

func upperOnly(p99, p95 float64) thresholds { return thresholds{p99, p95, math.NaN(), math.NaN()} }

func finite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

// ── public entry point ───────────────────────────────────────────────────────

func GeneratePlots(cfg utils.AnalysisConfig, bsaType string, smoothed []stats.SmoothedStats, th stats.Thresholds, brmBlocks []stats.BRMBlock) error {
	color.Cyan("\n================================ Generating plots for %s... ==========================================", bsaType)
	if len(smoothed) == 0 {
		return nil
	}

	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := stats.BulkFlags(bsaType)
	signedDSI := cfg.SignedDeltaSI(bsaType)
	hmmOn := stats.HMMAvailable(cfg, bsaType)

	// ── group by chromosome, in position order ──────────────────────────────

	byChr := make(map[string][]int)
	for i := range smoothed {
		byChr[smoothed[i].CHROM] = append(byChr[smoothed[i].CHROM], i)
	}
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// Index BRM blocks by chrom for fast lookup.
	brmByChr := make(map[string][]stats.BRMBlock)
	for _, b := range brmBlocks {
		brmByChr[b.Chrom] = append(brmByChr[b.Chrom], b)
	}

	// Thresholds are genome-wide: the same lines on every chromosome.
	tb, ob := th.TwoBulk, th.OneBulk
	valleys := func(p99, p95, m99, m95 float64) thresholds { return thresholds{p99, p95, m99, m95} }

	// ── create page containers ───────────────────────────────────────────────

	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics"

	consensusPage := components.NewPage()
	consensusPage.SetLayout(components.PageFlexLayout)
	consensusPage.PageTitle = "GoBSAseq — DsiGconsensusZ"

	// ── per-chromosome chart generation ─────────────────────────────────────

	for _, chrom := range chroms {
		idx := byChr[chrom]
		sort.SliceStable(idx, func(a, b int) bool { return smoothed[idx[a]].POS < smoothed[idx[b]].POS })
		chromBRM := brmByChr[chrom]

		x := make([]int64, len(idx))
		for i, j := range idx {
			x[i] = smoothed[j].POS
		}
		span := [2]float64{float64(x[0]), float64(x[len(x)-1])}
		track := func(get func(*stats.SmoothedStats) float64, decimals int) [][2]float64 {
			y := make([]float64, len(idx))
			for i, j := range idx {
				y[i] = get(&smoothed[j])
			}
			return simplify(x, y, decimals)
		}
		stat := func(get func(*stats.SmoothedStats) float64) [][2]float64 { return track(get, statDecimals) }

		if hasBothBulks {
			individualPage.AddCharts(siCharts(chrom, span,
				stat(func(s *stats.SmoothedStats) float64 { return s.SmHighSI }),
				stat(func(s *stats.SmoothedStats) float64 { return s.SmLowSI }),
				valleys(tb.HighSIP99, tb.HighSIP95, tb.HighSIMp99, tb.HighSIMp95),
				valleys(tb.LowSIP99, tb.LowSIP95, tb.LowSIMp99, tb.LowSIMp95), signedDSI, chromBRM)...)
			individualPage.AddCharts(
				dsiChart(chrom, span, stat(func(s *stats.SmoothedStats) float64 { return s.SmDeltaSI }),
					valleys(tb.DeltaSIP99, tb.DeltaSIP95, tb.DeltaSIMp99, tb.DeltaSIMp95), signedDSI, chromBRM),
				createRawValueChart(chrom+" G-statistic", span, stat(func(s *stats.SmoothedStats) float64 { return s.SmGstat }),
					upperOnly(tb.GstatP99, tb.GstatP95), false, chromBRM),
				createRawValueChart(chrom+" ED⁴", span, stat(func(s *stats.SmoothedStats) float64 { return s.SmED }),
					upperOnly(tb.ED4P99, tb.ED4P95), false, chromBRM),
			)
		}
		if hasHighBulk && !hasBothBulks {
			individualPage.AddCharts(
				createRawValueChart(chrom+" HighSI", span, stat(func(s *stats.SmoothedStats) float64 { return s.SmHighSI }),
					valleys(tb.HighSIP99, tb.HighSIP95, tb.HighSIMp99, tb.HighSIMp95), true, chromBRM),
			)
		}
		if hasLowBulk && !hasBothBulks {
			individualPage.AddCharts(
				createRawValueChart(chrom+" LowSI", span, stat(func(s *stats.SmoothedStats) float64 { return s.SmLowSI }),
					valleys(tb.LowSIP99, tb.LowSIP95, tb.LowSIMp99, tb.LowSIMp95), true, chromBRM),
			)
		}
		if hasOneBulk {
			individualPage.AddCharts(
				createRawValueChart(chrom+" AF Deviation", span, stat(func(s *stats.SmoothedStats) float64 { return s.SmAFDev }),
					valleys(ob.AFDevP99, ob.AFDevP95, ob.AFDevMp99, ob.AFDevMp95), true, chromBRM),
				createRawValueChart(chrom+" One-bulk G-statistic", span, stat(func(s *stats.SmoothedStats) float64 { return s.SmOneBulkG }),
					upperOnly(ob.OneBulkGstatP99, ob.OneBulkGstatP95), false, chromBRM),
			)
		}
		if hmmOn {
			individualPage.AddCharts(
				createRawValueChart(chrom+" HMM LOD (genome-wide thresholds)", span, stat(func(s *stats.SmoothedStats) float64 { return s.HMMLOD }),
					upperOnly(th.HMM.LODP99, th.HMM.LODP95), false, chromBRM),
			)
		}

		// ── DsiGconsensusZ chart ─────────────────────────────────────────────
		consensusPage.AddCharts(
			createConsensusChart(chrom, span, track(func(s *stats.SmoothedStats) float64 { return s.DsiGconsensusZ }, consensusDecimals),
				chromBRM, th.Consensus, hasOneBulk),
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
		{consensusPage, fmt.Sprintf("GoBSAseq.%s.dsigconsensusz.html", bsaType)},
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

// ── simplification ───────────────────────────────────────────────────────────

// plotTolerance is the largest vertical distance, as a fraction of a
// chromosome track's range, between any marker and the plotted line: 0.1%,
// under half a pixel on a 500 px chart.
const plotTolerance = 0.001

// simplify reduces one chromosome's track (x sorted) to the [position, value]
// pairs of a line that passes within plotTolerance × the track's range of
// every marker (Ramer–Douglas–Peucker on vertical distance): a split keeps
// the marker farthest from the chord until every marker is within tolerance.
// Peaks, valleys and threshold crossings are therefore drawn to within the
// tolerance, while stretches the line already follows (the smoothed tracks
// are linear between smoothing-grid points where markers are dense) need no
// points. The first and last markers and the track's largest and smallest
// values are always kept, so the tooltip shows the exact peak. Values are
// rounded to decimals places; non-finite values are dropped.
func simplify(x []int64, y []float64, decimals int) [][2]float64 {
	xs, ys := make([]float64, 0, len(x)), make([]float64, 0, len(y))
	imin, imax := -1, -1
	for i := range x {
		if finite(y[i]) {
			if imin < 0 || y[i] < ys[imin] {
				imin = len(ys)
			}
			if imax < 0 || y[i] > ys[imax] {
				imax = len(ys)
			}
			xs, ys = append(xs, float64(x[i])), append(ys, y[i])
		}
	}
	n := len(xs)
	if n == 0 {
		return nil
	}
	eps := plotTolerance * (ys[imax] - ys[imin])
	keep := make([]bool, n)
	// Split first at the extremes, so they are exact; the segments between the
	// kept points are then refined.
	cuts := []int{0, n - 1, imin, imax}
	slices.Sort(cuts)
	var stack [][2]int
	for k := range cuts {
		keep[cuts[k]] = true
		if k > 0 && cuts[k] > cuts[k-1] {
			stack = append(stack, [2]int{cuts[k-1], cuts[k]})
		}
	}
	for len(stack) > 0 {
		a, b := stack[len(stack)-1][0], stack[len(stack)-1][1]
		stack = stack[:len(stack)-1]
		if b-a < 2 {
			continue
		}
		slope := 0.0
		if dx := xs[b] - xs[a]; dx > 0 {
			slope = (ys[b] - ys[a]) / dx
		}
		worst, wi := eps, -1
		for i := a + 1; i < b; i++ {
			if d := math.Abs(ys[i] - ys[a] - slope*(xs[i]-xs[a])); d > worst {
				worst, wi = d, i
			}
		}
		if wi >= 0 {
			keep[wi] = true
			stack = append(stack, [2]int{a, wi}, [2]int{wi, b})
		}
	}
	scale := math.Pow(10, float64(decimals))
	var out [][2]float64
	for i := range keep {
		if keep[i] {
			out = append(out, [2]float64{xs[i], math.Round(ys[i]*scale) / scale})
		}
	}
	return out
}

// ── shared chart options ─────────────────────────────────────────────────────

// jsNum formats v for JavaScript (null when not finite).
func jsNum(v float64) string {
	if !finite(v) {
		return "null"
	}
	return strconv.FormatFloat(v, 'g', -1, 64)
}

// yBounds returns the y-axis min and max: a non-negative chart starts at 0;
// otherwise each end is automatic unless a threshold line lies beyond the
// data, in which case the axis extends to it (plus 5% of the range), so the
// lines stay in view when zoomed into a stretch that does not reach them.
func yBounds(t thresholds, bidirectional bool) (yMin, yMax interface{}) {
	top, bottom := math.NaN(), math.NaN()
	for _, v := range []float64{t.p99, t.p95} {
		if finite(v) && !(v <= top) {
			top = v
		}
	}
	for _, v := range []float64{t.m99, t.m95} {
		if finite(v) && !(v >= bottom) {
			bottom = v
		}
	}
	lo := "Math.min(v.min, 0)"
	if bidirectional {
		lo = "v.min"
		if finite(bottom) {
			lo = "Math.min(v.min, " + jsNum(bottom) + ")"
		}
	}
	hi := "v.max"
	if finite(top) {
		hi = "Math.max(v.max, " + jsNum(top) + ")"
	}
	if !bidirectional {
		yMin = opts.Float(0.0)
	} else if finite(bottom) {
		yMin = opts.FuncOpts(fmt.Sprintf("function(v){var lo=%s,hi=%s;return lo-0.05*(hi-lo);}", lo, hi))
	}
	if finite(top) {
		yMax = opts.FuncOpts(fmt.Sprintf("function(v){var lo=%s,hi=%s;return hi+0.05*(hi-lo);}", lo, hi))
	}
	return yMin, yMax
}

// commonGlobalOpts returns the standard set of go-echarts global options used
// by every chart in GoBSAseq: a linear base-pair x axis spanning the data,
// and y bounds from yBounds.
func commonGlobalOpts(title, subtitle, yLabel, width, height string, t thresholds, bidirectional bool) []charts.GlobalOpts {
	yMin, yMax := yBounds(t, bidirectional)
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
			Max:          yMax,
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

// tooltipStyle is how a chart's tooltip reports the value at the pointer.
type tooltipStyle struct {
	decimals   int
	inclusive  bool   // a value on the line counts as beyond it
	prefix     string // before the position
	sig99      string // label of a value beyond a p99 line
	sig95      string
	thresholds bool // list the threshold values under the statistic
}

// tooltipJS returns the axis tooltip for a chart whose data series is main:
// the position and value of the nearest plotted point, ★/● when the value is
// beyond a p99/p95 line of either tail (upper lines: above; valley lines:
// below), and the threshold values.
func tooltipJS(main string, t thresholds, st tooltipStyle) string {
	gt, lt := ">", "<"
	if st.inclusive {
		gt, lt = ">=", "<="
	}
	list := ""
	if st.thresholds {
		list = fmt.Sprintf(`
		[['p99', T.u99, '%[1]s'], ['p95', T.u95, '%[2]s'], ['p99 valley', T.l99, '%[1]s'], ['p95 valley', T.l95, '%[2]s']].forEach(function(r) {
			if (r[1] !== null) out += '<span style="display:inline-block;margin-right:4px;border-radius:10px;width:10px;height:10px;background-color:' + r[2] + ';"></span>' + r[0] + ': ' + r[1].toFixed(D) + '<br/>';
		});`, colorP99, colorP95)
	}
	return fmt.Sprintf(`function(params) {
		var T = {u99: %[1]s, u95: %[2]s, l99: %[3]s, l95: %[4]s}, D = %[5]d;
		var main = null;
		params.forEach(function(p) { if (p.seriesName === '%[6]s') main = p; });
		if (!main) return '';
		var pos = main.value[0], val = main.value[1];
		var posFmt = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
		var sig = '';
		if ((T.u99 !== null && val %[7]s T.u99) || (T.l99 !== null && val %[8]s T.l99)) sig = ' <span style="color:%[9]s;font-weight:bold">%[10]s</span>';
		else if ((T.u95 !== null && val %[7]s T.u95) || (T.l95 !== null && val %[8]s T.l95)) sig = ' <span style="color:%[11]s">%[12]s</span>';
		var out = '<strong>%[13]s' + posFmt + '</strong><br/>' + main.marker + ' ' + main.seriesName + ': ' + val.toFixed(D) + sig + '<br/>';%[14]s
		return out;
	}`, jsNum(t.p99), jsNum(t.p95), jsNum(t.m99), jsNum(t.m95), st.decimals, main, gt, lt,
		colorP99, st.sig99, colorP95, st.sig95, st.prefix, list)
}

// addDataSeries adds a [position, value] line series.
func addDataSeries(line *charts.Line, name string, pts [][2]float64, options ...charts.SeriesOpts) {
	s := charts.SingleSeries{Name: name, Type: types.ChartLine, Data: pts}
	s.InitSeriesDefaultOpts(line.BaseConfiguration)
	s.ConfigureSeriesOpts(options...)
	line.MultiSeries = append(line.MultiSeries, s)
}

// addThreshold adds a horizontal threshold line at y as its own legend entry
// (toggling it hides the line): a mark line, which spans the plot at any zoom
// and takes no per-marker points, over a two-point series at the chromosome
// ends in the same style. No line when y is not finite.
func addThreshold(line *charts.Line, name string, y float64, span [2]float64, col string, width float32) {
	if !finite(y) {
		return
	}
	style := opts.LineStyle{Type: "dashed", Width: width, Color: col}
	addDataSeries(line, name, [][2]float64{{span[0], y}, {span[1], y}},
		charts.WithLineStyleOpts(style),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: col, Opacity: opts.Float(0)}),
		charts.WithMarkLineNameYAxisItemOpts(opts.MarkLineNameYAxisItem{Name: name, YAxis: y}),
		charts.WithMarkLineStyleOpts(opts.MarkLineStyle{
			Symbol:    []string{"none", "none"},
			Label:     &opts.Label{Show: opts.Bool(false)},
			LineStyle: &style,
		}),
	)
}

// addThresholds adds a chart's p99/p95 lines and, for two-tailed tracks, its
// valley lines.
func addThresholds(line *charts.Line, t thresholds, span [2]float64) {
	addThreshold(line, "p99", t.p99, span, colorP99, 1.8)
	addThreshold(line, "p95", t.p95, span, colorP95, 1.4)
	addThreshold(line, "p99 valley", t.m99, span, colorP99, 1.8)
	addThreshold(line, "p95 valley", t.m95, span, colorP95, 1.4)
}

// ── individual raw-value chart ───────────────────────────────────────────────

// siCharts draws the two bulks' SI. Without polarity (2b, --ref-parent none)
// SI is the REF fraction: it is labelled as such and drawn without
// thresholds, because a QTL's signal cancels across SNPs there.
func siCharts(chrom string, span [2]float64, hi, li [][2]float64, hiT, liT thresholds, signed bool, brm []stats.BRMBlock) []components.Charter {
	if !signed {
		return []components.Charter{
			createRawValueChart(chrom+" High bulk REF fraction (polarity unknown: not a QTL track)", span, hi, noThresholds(), false, brm),
			createRawValueChart(chrom+" Low bulk REF fraction (polarity unknown: not a QTL track)", span, li, noThresholds(), false, brm),
		}
	}
	return []components.Charter{
		createRawValueChart(chrom+" HighSI", span, hi, hiT, true, brm),
		createRawValueChart(chrom+" LowSI", span, li, liT, true, brm),
	}
}

// dsiChart draws ΔSI, or |ΔSI| (upper thresholds only) when polarity is unknown.
func dsiChart(chrom string, span [2]float64, pts [][2]float64, t thresholds, signed bool, brm []stats.BRMBlock) *charts.Line {
	if signed {
		return createRawValueChart(chrom+" ΔSI (DeltaSI)", span, pts, t, true, brm)
	}
	return createRawValueChart(chrom+" |ΔSI| (AbsDeltaSI, polarity unknown)", span, pts, upperOnly(t.p99, t.p95), false, brm)
}

// createRawValueChart builds one smoothed statistic's line chart with its
// genome-wide p99/p95 lines (valley lines for two-tailed tracks) and the BRM
// blocks shaded. bidirectional leaves the y axis free below 0.
func createRawValueChart(title string, span [2]float64, pts [][2]float64, t thresholds, bidirectional bool, brmBlocks []stats.BRMBlock) *charts.Line {
	subtitle := "no thresholds  |  shaded: BRM blocks"
	if finite(t.p99) {
		subtitle = fmt.Sprintf("p99 threshold: %.4f  |  p95 threshold: %.4f", t.p99, t.p95)
		if finite(t.m99) {
			subtitle += fmt.Sprintf("  |  valleys: %.4f, %.4f", t.m99, t.m95)
		}
		subtitle += "  |  shaded: BRM blocks"
	}

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Value", chartWidth, chartHeight, t, bidirectional)...)
	line.SetGlobalOptions(charts.WithTooltipOpts(opts.Tooltip{
		Show:        opts.Bool(true),
		Trigger:     "axis",
		AxisPointer: &opts.AxisPointer{Type: "cross"},
		Formatter: opts.FuncOpts(tooltipJS("Statistic", t, tooltipStyle{decimals: statDecimals, prefix: "Position: ",
			sig99: "★ p99", sig95: "● p95", thresholds: true})),
	}))

	statSeriesOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#1f77b4"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	statSeriesOpts = append(statSeriesOpts, brmMarkAreaOpts(brmBlocks)...)
	addDataSeries(line, "Statistic", pts, statSeriesOpts...)
	addThresholds(line, t, span)
	return line
}

// ── DsiGconsensusZ chart ──────────────────────────────────────────────────────

// createConsensusChart renders DsiGconsensusZ (the null-standardised ΔSI and
// G′ evidence combined by correlation-corrected Stouffer; dsigconsensus.go)
// with its genome-wide p99 / p95 lines from simulated no-QTL genomes, and the
// BRM blocks shaded.
func createConsensusChart(chrom string, span [2]float64, pts [][2]float64, brmBlocks []stats.BRMBlock, th stats.ConsensusThresholds, oneBulk bool) *charts.Line {
	t := upperOnly(th.P99, th.P95)
	parts := "ΔSI + G′"
	if oneBulk {
		parts = "AFDev + G1"
	}

	title := chrom + " — DsiGconsensusZ"
	subtitle := fmt.Sprintf(
		"%s, each standardised against simulated no-QTL genomes; Stouffer corrected for their null correlation r = %.2f. "+
			"Dashed lines: genome-wide p99 (%.3f) and p95 (%.3f). Shaded bands: BRM blocks.",
		parts, th.R, th.P99, th.P95,
	)

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "DsiGconsensusZ", chartWidth, chartHeight, t, true)...)
	line.SetGlobalOptions(
		charts.WithTooltipOpts(opts.Tooltip{
			Show:        opts.Bool(true),
			Trigger:     "axis",
			AxisPointer: &opts.AxisPointer{Type: "cross"},
			Formatter: opts.FuncOpts(tooltipJS("DsiGconsensusZ", t, tooltipStyle{decimals: consensusDecimals, inclusive: true,
				sig99: "★ p99 significant", sig95: "● p95 suggestive", thresholds: true})),
		}),
	)

	seriesOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#2ca02c"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	seriesOpts = append(seriesOpts, brmMarkAreaOpts(brmBlocks)...)
	addDataSeries(line, "DsiGconsensusZ", pts, seriesOpts...)
	addThreshold(line, fmt.Sprintf("p95 (%.3f)", th.P95), th.P95, span, colorP95, 1.4)
	addThreshold(line, fmt.Sprintf("p99 (%.3f)", th.P99), th.P99, span, colorP99, 1.8)
	return line
}

// ── BRM mark-area helpers ─────────────────────────────────────────────────────

// brmMarkAreaOpts returns go-echarts SeriesOpts that shade the BRM blocks,
// [Start, Stop] in bp on the position axis, in the same amber colour used by
// the legacy code. It is safe to call with nil or empty blocks.
func brmMarkAreaOpts(blocks []stats.BRMBlock) []charts.SeriesOpts {
	if len(blocks) == 0 {
		return nil
	}
	areas := make([][]opts.MarkAreaData, 0, len(blocks))
	for _, b := range blocks {
		areas = append(areas, []opts.MarkAreaData{
			{Name: fmt.Sprintf("BRM block %.4f", b.Peak), XAxis: float64(b.Start)},
			{XAxis: float64(b.Stop)},
		})
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
