package plots

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
)

const (
	chartTheme   = types.ThemeWesteros
	chartWidth   = "900px"
	chartHeight  = "500px"
	zSig         = 3.0
	zSugg        = 2.0
	posFormatter = `function(value) {
		if (value >= 1000000) { return (value / 1000000).toFixed(2) + ' Mb'; }
		if (value >= 1000)    { return (value / 1000).toFixed(1)    + ' kb'; }
		return value;
	}`
)

// CreatePlots generates all HTML plots for the analysis
func CreatePlots(cfg utils.AnalysisConfig, bsaType string, sm []stats.SmoothedStats, brmBlocks []stats.BRMBlock) error {
	if len(sm) == 0 {
		return fmt.Errorf("no smoothed stats for plotting")
	}

	outDir := filepath.Join(cfg.OutputDir, "plots")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return err
	}

	hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk := bulkFlags(bsaType)

	// Group by chromosome
	byChrom := make(map[string][]stats.SmoothedStats)
	for _, s := range sm {
		byChrom[s.CHROM] = append(byChrom[s.CHROM], s)
	}

	chroms := make([]string, 0, len(byChrom))
	for c := range byChrom {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	// Convert BRM blocks for plotting
	plotBlocks := convertBRMBlocks(brmBlocks)

	// Create pages
	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics"

	robustZPage := components.NewPage()
	robustZPage.SetLayout(components.PageFlexLayout)
	robustZPage.PageTitle = "GoBSAseq — Robust Z-score Overlay"

	compositePage := components.NewPage()
	compositePage.SetLayout(components.PageFlexLayout)
	compositePage.PageTitle = "GoBSAseq — Composite Signal"

	for _, chrom := range chroms {
		stats := byChrom[chrom]
		if len(stats) == 0 {
			continue
		}

		// Extract data arrays
		n := len(stats)
		x := make([]int64, n)
		hi := make([]float64, n)
		li := make([]float64, n)
		dsi := make([]float64, n)
		gs := make([]float64, n)
		ed := make([]float64, n)
		lod := make([]float64, n)
		bbl := make([]float64, n)

		for i, s := range stats {
			x[i] = s.POS
			if hasHighBulk {
				hi[i] = s.SmHighSI
			}
			if hasLowBulk {
				li[i] = s.SmLowSI
			}
			if hasBothBulks {
				dsi[i] = s.SmDeltaSI
				gs[i] = s.SmGstat
				ed[i] = s.SmED
				lod[i] = s.SmLOD
				bbl[i] = s.SmBBLogBF
			}
			if hasOneBulk {
				hi[i] = s.SmAFDev
			}
		}

		// Compute Z-scores for overlay if using both bulks
		if hasBothBulks {
			// Use composite and maxAbsZ for plotting
			individualPage.AddCharts(
				createInteractiveLineChart(chrom+" DeltaSI", x, dsi, plotBlocks, true),
				createInteractiveLineChart(chrom+" Gstat", x, gs, plotBlocks, false),
				createInteractiveLineChart(chrom+" ED4", x, ed, plotBlocks, false),
				createInteractiveLineChart(chrom+" LOD", x, lod, plotBlocks, false),
				createInteractiveLineChart(chrom+" BBLogBF", x, bbl, plotBlocks, false),
			)

			// Robust Z overlay
			robustZPage.AddCharts(createRobustZOverlayChart(chrom, x,
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZHighSI }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZLowSI }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZDeltaSI }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZGstat }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZED }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZLOD }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZBBLogBF }),
				plotBlocks))

			// Composite signal
			compositePage.AddCharts(createCompositeSignalChart(chrom, x,
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZHighSI }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZLowSI }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZDeltaSI }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZGstat }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZED }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZLOD }),
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZBBLogBF }),
				plotBlocks))
		}

		if hasOneBulk {
			individualPage.AddCharts(
				createInteractiveLineChart(chrom+" AFDev", x, hi, plotBlocks, true),
				createInteractiveLineChart(chrom+" OneBulkG", x, statsToFloat(stats, func(s stats.SmoothedStats) float64 { return s.SmOneBulkG }), plotBlocks, false),
				createInteractiveLineChart(chrom+" OneBulkLOD", x, statsToFloat(stats, func(s stats.SmoothedStats) float64 { return s.SmOneBulkLOD }), plotBlocks, false),
				createInteractiveLineChart(chrom+" OneBulkBBLogBF", x, statsToFloat(stats, func(s stats.SmoothedStats) float64 { return s.SmOneBulkBBLogBF }), plotBlocks, false),
			)

			robustZPage.AddCharts(createRobustZOverlayChart(chrom, x,
				statsToZ(stats, func(s stats.SmoothedStats) float64 { return s.ZAFDev }),
				nil, nil, nil, nil, nil, nil, plotBlocks))
		}
	}

	// Write HTML files
	if err := writeHTMLPage(individualPage, filepath.Join(outDir, "GoBSAseq_IndividualPlots.html")); err != nil {
		return err
	}
	if err := writeHTMLPage(robustZPage, filepath.Join(outDir, "GoBSAseq_RobustZScore.html")); err != nil {
		return err
	}
	if err := writeHTMLPage(compositePage, filepath.Join(outDir, "GoBSAseq_CompositeSignal.html")); err != nil {
		return err
	}

	return nil
}

// PlotBlock is a simplified block for plotting
type PlotBlock struct {
	Chrom   string
	Start   int64
	Stop    int64
	PeakPos int64
	Peak    float64
}

func convertBRMBlocks(blocks []stats.BRMBlock) []PlotBlock {
	result := make([]PlotBlock, len(blocks))
	for i, b := range blocks {
		result[i] = PlotBlock{
			Chrom:   b.Chrom,
			Start:   b.Start,
			Stop:    b.Stop,
			PeakPos: b.PeakPos,
			Peak:    b.Peak,
		}
	}
	return result
}

func statsToFloat(stats []stats.SmoothedStats, fn func(stats.SmoothedStats) float64) []float64 {
	result := make([]float64, len(stats))
	for i, s := range stats {
		result[i] = fn(s)
	}
	return result
}

func statsToZ(stats []stats.SmoothedStats, fn func(stats.SmoothedStats) float64) []float64 {
	return statsToFloat(stats, fn)
}

func createInteractiveLineChart(title string, x []int64, y []float64, blocks []PlotBlock, hasNegativeThresh bool) *charts.Line {
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, "Statistic value", chartWidth, chartHeight, hasNegativeThresh)...)

	line.SetGlobalOptions(chartsWithThresholdTooltip())

	n := len(y)
	yData := make([]opts.LineData, n)
	for i, v := range y {
		yData[i] = opts.LineData{Value: v}
	}

	line.SetXAxis(positionLabels(x)).
		AddSeries("Statistic", yData,
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#1f77b4"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	if len(blocks) > 0 {
		line = line.AddSeries("BRM Blocks", make([]opts.LineData, n),
			charts.WithMarkAreaData(brmBlockMarkAreas(blocks, x)...),
			charts.WithMarkAreaStyleOpts(opts.MarkAreaStyle{
				Label:     &opts.Label{Show: opts.Bool(false)},
				ItemStyle: &opts.ItemStyle{Color: "rgba(243, 156, 18, 0.22)"},
			}),
		)
	}

	return line
}

func createRobustZOverlayChart(chrom string, x []int64,
	hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ []float64,
	blocks []PlotBlock) *charts.Line {

	title := chrom + " — Robust Z-score Overlay"
	subtitle := "Genome-wide robust Z-score. z = ±2 suggestive · z = ±3 significant"

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
		chartsWithZTooltip(),
	)

	n := len(x)
	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}

	line.SetXAxis(positionLabels(x)).
		AddSeries("z=0", mkRef(0),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#bdc3c7"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
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

	series := []struct {
		name  string
		data  []float64
		color string
		width float32
	}{
		{"HighSI", hiZ, "#1f77b4", 2.0},
		{"LowSI", liZ, "#ff7f0e", 2.0},
		{"DeltaSI", dsiZ, "#2ca02c", 3.0},
		{"Gstat", gsZ, "#17becf", 2.0},
		{"ED", edZ, "#d62728", 2.0},
		{"LOD", lodZ, "#9467bd", 2.0},
		{"BBLogBF", bblZ, "#8c564b", 2.0},
	}

	for _, s := range series {
		if s.data == nil {
			continue
		}
		ld := floatSliceToLineData(s.data)
		line.AddSeries(s.name, ld,
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: s.width, Color: s.color}),
			charts.WithItemStyleOpts(opts.ItemStyle{Color: s.color, Opacity: opts.Float(0)}),
		)
	}

	if len(blocks) > 0 {
		line = line.AddSeries("BRM Blocks", make([]opts.LineData, n),
			charts.WithMarkAreaData(brmBlockMarkAreas(blocks, x)...),
			charts.WithMarkAreaStyleOpts(opts.MarkAreaStyle{
				Label:     &opts.Label{Show: opts.Bool(false)},
				ItemStyle: &opts.ItemStyle{Color: "rgba(243, 156, 18, 0.22)"},
			}),
		)
	}

	return line
}

func createCompositeSignalChart(chrom string, x []int64,
	hiZ, liZ, dsiZ, gsZ, edZ, lodZ, bblZ []float64,
	blocks []PlotBlock) *charts.Line {

	n := len(x)
	composite := make([]float64, n)
	for i := range composite {
		composite[i] = maxAbs(hiZ[i], liZ[i], dsiZ[i], gsZ[i], edZ[i], lodZ[i], bblZ[i])
	}

	title := chrom + " — Composite Signal (max |Z|)"
	subtitle := "Max absolute robust Z-score across all statistics. z=2 suggestive · z=3 significant"

	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "max |Z-score|", chartWidth, chartHeight, false)...)
	line.SetGlobalOptions(chartsWithZTooltip())

	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}

	compositeData := floatSliceToLineData(composite)

	line.SetXAxis(positionLabels(x)).
		AddSeries("Composite", compositeData,
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#2ca02c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=2 (sugg.)", mkRef(zSugg),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		).
		AddSeries("z=3 (sig.)", mkRef(zSig),
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)

	if len(blocks) > 0 {
		line = line.AddSeries("BRM Blocks", make([]opts.LineData, n),
			charts.WithMarkAreaData(brmBlockMarkAreas(blocks, x)...),
			charts.WithMarkAreaStyleOpts(opts.MarkAreaStyle{
				Label:     &opts.Label{Show: opts.Bool(false)},
				ItemStyle: &opts.ItemStyle{Color: "rgba(243, 156, 18, 0.22)"},
			}),
		)
	}

	return line
}

// Helper functions

func maxAbs(vals ...float64) float64 {
	max := 0.0
	for _, v := range vals {
		if a := abs(v); a > max {
			max = a
		}
	}
	return max
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

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

func chartsWithThresholdTooltip() charts.GlobalOpts {
	return charts.WithTooltipOpts(opts.Tooltip{
		Show:        opts.Bool(true),
		Trigger:     "axis",
		AxisPointer: &opts.AxisPointer{Type: "cross"},
	})
}

func chartsWithZTooltip() charts.GlobalOpts {
	return charts.WithTooltipOpts(opts.Tooltip{
		Show:        opts.Bool(true),
		Trigger:     "axis",
		AxisPointer: &opts.AxisPointer{Type: "cross"},
		Formatter: opts.FuncOpts(`function(params) {
			let pos = params[0].axisValue;
			let posStr = pos >= 1e6 ? (pos/1e6).toFixed(3)+' Mb' : pos >= 1000 ? (pos/1000).toFixed(2)+' kb' : pos+' bp';
			let result = '<strong>' + posStr + '</strong><br/>';
			params.forEach(function(item) {
				let val = parseFloat(item.value);
				if (isNaN(val)) return;
				let sig = '';
				if (Math.abs(val) >= 3.0)      sig = ' <span style="color:#e74c3c;font-weight:bold">★ significant</span>';
				else if (Math.abs(val) >= 2.0)  sig = ' <span style="color:#f39c12">● suggestive</span>';
				result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
			});
			return result;
		}`),
	})
}

func brmBlockMarkAreas(blocks []PlotBlock, x []int64) [][]opts.MarkAreaData {
	if len(blocks) == 0 || len(x) == 0 {
		return nil
	}

	xLabels := positionLabels(x)
	areas := make([][]opts.MarkAreaData, 0, len(blocks))
	for _, b := range blocks {
		startIdx := sort.Search(len(x), func(i int) bool { return x[i] >= b.Start })
		stopIdx := sort.Search(len(x), func(i int) bool { return x[i] > b.Stop }) - 1
		if startIdx >= len(x) || stopIdx < 0 {
			continue
		}
		if stopIdx < startIdx {
			stopIdx = startIdx
		}

		areas = append(areas, []opts.MarkAreaData{
			{
				Name:  fmt.Sprintf("BRM block %.4f", b.Peak),
				XAxis: xLabels[startIdx],
			},
			{XAxis: xLabels[stopIdx]},
		})
	}
	return areas
}

func positionLabels(x []int64) []string {
	labels := make([]string, len(x))
	for i, v := range x {
		labels[i] = fmt.Sprintf("%d", v)
	}
	return labels
}

func floatSliceToLineData(vals []float64) []opts.LineData {
	ld := make([]opts.LineData, len(vals))
	for i, v := range vals {
		ld[i] = opts.LineData{Value: v}
	}
	return ld
}

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

func bulkFlags(bsaType string) (hasHighBulk, hasLowBulk, hasBothBulks, hasOneBulk bool) {
	hasHighBulk = strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLowBulk = strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	hasBothBulks = hasHighBulk && hasLowBulk
	hasOneBulk = hasHighBulk != hasLowBulk
	return
}
