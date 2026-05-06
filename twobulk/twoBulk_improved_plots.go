// This file contains the refactored plotting functions for GoBSAseq
// Changes:
// 1. GenerateHtmlPlotsAndQTL now creates TWO separate HTML files:
//    - GoBSAseq_IndividualPlots.html (individual stat plots per chromosome)
//    - GoBSAseq_NormalizedOverlay.html (normalized overlay plots per chromosome)
// 2. Improved layout with better page organization, navigation, and styling
// 3. Enhanced chart dimensions and responsive design
// 4. Better color scheme and visual hierarchy

// ============== REFACTORED GenerateHtmlPlotsAndQTL ==============

// GenerateHtmlPlotsAndQTL processes smoothed stats, detects QTLs, and creates interactive HTML charts
// Now generates TWO separate HTML files for better organization
package twoBulk

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"sync"

	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
)

// BSAstats holds the raw statistics for a single SNP
type BSAstats struct {
	CHROM      string
	POS        int64
	REF        string
	ALT        string
	HighParGT  []int
	LowParGT   []int
	HighBulkGT []int
	HighBulkAD string
	LowBulkGT  []int
	LowBulkAD  string

	HighBulkL int
	HighBulkH int
	LowBulkL  int
	LowBulkH  int
	HighSI    float64
	LowSI     float64

	DeltaSI float64
	Gstat   float64
	ED      float64
	LOD     float64
	BBLogBF float64

	Depth int
}

// SmoothedStats holds the averaged statistics for a genomic window
type SmoothedStats struct {
	CHROM          string
	POS            int64
	DeltaSI        float64
	Gstat          float64
	ED             float64
	LOD            float64
	BBLogBF        float64
	HighSI         float64
	LowSI          float64
	NumSNPs        int
	MeanHighBulkDP int
	MeanLowBulkDP  int
}

// Thresholds holds the significance levels for each statistic
type Thresholds struct {
	DsiP99  float64
	DsiP95  float64
	DsiMp99 float64
	DsiMp95 float64

	GsP99 float64
	GsP95 float64

	EdP99 float64
	EdP95 float64

	LodP99 float64
	LodP95 float64

	BbP99 float64
	BbP95 float64

	HighP99  float64
	HighP95  float64
	HighMp99 float64
	HighMp95 float64

	LowP99  float64
	LowP95  float64
	LowMp99 float64
	LowMp95 float64
}

// QTLRecord holds the detected QTL bounds and peaks
type QTLRecord struct {
	Chrom string
	Start int64
	Stop  int64
	Peak  float64
	Stat  string
	CI    string
}

var thresholdCache sync.Map

func GenerateHtmlPlotsAndQTL(allSmoothed []SmoothedStats, highSmAF, lowSmAF float64, rep int, htmlOutFile, qtlOutFile string) error {
	// 1. Group stats by chromosome
	byChr := make(map[string][]SmoothedStats)
	for _, s := range allSmoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}

	var chroms []string
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	var allQTLs []QTLRecord

	// ============================================================
	// PAGE 1: Individual Statistics Plots
	// ============================================================
	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq - Individual Statistics"

	// ============================================================
	// PAGE 2: Normalized Overlay Plots
	// ============================================================
	normalizedPage := components.NewPage()
	normalizedPage.SetLayout(components.PageFlexLayout)
	normalizedPage.PageTitle = "GoBSAseq - Normalized Overlay"

	for _, chrom := range chroms {
		stats := byChr[chrom]
		n := float64(len(stats))
		if n == 0 {
			continue
		}

		// 2. Calculate average thresholds for the current chromosome
		var sumHp99, sumHp95, sumLp99, sumLp95, sumDp99, sumDp95, sumDMp99, sumDMp95, sumEp99, sumEp95, sumLodp99, sumLodp95, sumBbp99, sumBbp95 float64

		for _, s := range stats {
			t := calcThresholdsCached(s.MeanHighBulkDP, s.MeanLowBulkDP, highSmAF, lowSmAF, rep)
			sumHp99 += t.HighP99
			sumHp95 += t.HighP95
			sumLp99 += t.LowP99
			sumLp95 += t.LowP95
			sumDp99 += t.DsiP99
			sumDp95 += t.DsiP95
			sumDMp99 += t.DsiMp99
			sumDMp95 += t.DsiMp95
			sumEp99 += t.EdP99
			sumEp95 += t.EdP95
			sumLodp99 += t.LodP99
			sumLodp95 += t.LodP95
			sumBbp99 += t.BbP99
			sumBbp95 += t.BbP95
		}

		avgHp99, avgHp95 := sumHp99/n, sumHp95/n
		avgLp99, avgLp95 := sumLp99/n, sumLp95/n
		avgDp99, avgDp95 := sumDp99/n, sumDp95/n
		avgDMp99, avgDMp95 := sumDMp99/n, sumDMp95/n
		avgEp99, avgEp95 := sumEp99/n, sumEp95/n
		avgLodp99, avgLodp95 := sumLodp99/n, sumLodp95/n
		avgBbp99, avgBbp95 := sumBbp99/n, sumBbp95/n

		// 3. Extract arrays for plotting and QTL detection
		var x []int
		var hi, li, dsi, ed, lod, bbl []float64

		for _, s := range stats {
			x = append(x, int(s.POS))
			hi = append(hi, s.HighSI)
			li = append(li, s.LowSI)
			dsi = append(dsi, s.DeltaSI)
			ed = append(ed, s.ED)
			lod = append(lod, s.LOD)
			bbl = append(bbl, s.BBLogBF)
		}

		// 4. Detect QTLs using the intersection logic
		allQTLs = append(allQTLs, detectQTLs(chrom, x, hi, avgHp99, "HighSI", "99", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, hi, avgHp95, "HighSI", "95", false)...)

		allQTLs = append(allQTLs, detectQTLs(chrom, x, li, avgLp99, "LowSI", "99", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, li, avgLp95, "LowSI", "95", false)...)

		allQTLs = append(allQTLs, detectQTLs(chrom, x, dsi, avgDp99, "DeltaSI", "99", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, dsi, avgDp95, "DeltaSI", "95", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, dsi, avgDMp99, "DeltaSI", "99", true)...) // Valleys
		allQTLs = append(allQTLs, detectQTLs(chrom, x, dsi, avgDMp95, "DeltaSI", "95", true)...) // Valleys

		allQTLs = append(allQTLs, detectQTLs(chrom, x, ed, avgEp99, "ED", "99", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, ed, avgEp95, "ED", "95", false)...)

		allQTLs = append(allQTLs, detectQTLs(chrom, x, lod, avgLodp99, "LOD", "99", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, lod, avgLodp95, "LOD", "95", false)...)

		allQTLs = append(allQTLs, detectQTLs(chrom, x, bbl, avgBbp99, "BBLogBF", "99", false)...)
		allQTLs = append(allQTLs, detectQTLs(chrom, x, bbl, avgBbp95, "BBLogBF", "95", false)...)

		// 5. Create individual Echarts (for Individual Plots page)
		hiChart := createInteractiveLineChart(chrom+" HighSI", x, hi, avgHp99, avgHp95, 0, 0, false)
		liChart := createInteractiveLineChart(chrom+" LowSI", x, li, avgLp99, avgLp95, 0, 0, false)
		dsiChart := createInteractiveLineChart(chrom+" DeltaSI", x, dsi, avgDp99, avgDp95, avgDMp99, avgDMp95, true)
		edChart := createInteractiveLineChart(chrom+" ED", x, ed, avgEp99, avgEp95, 0, 0, false)
		lodChart := createInteractiveLineChart(chrom+" LOD", x, lod, avgLodp99, avgLodp95, 0, 0, false)
		bblChart := createInteractiveLineChart(chrom+" BBLogBF", x, bbl, avgBbp99, avgBbp95, 0, 0, false)

		individualPage.AddCharts(hiChart, liChart, dsiChart, edChart, lodChart, bblChart)

		// 6. Create normalized overlay chart (for Normalized Overlay page)
		normChart := createNormalizedOverlayChart(chrom, x, hi, li, dsi, ed, lod, bbl,
			avgHp99, avgHp95, avgLp99, avgLp95,
			avgDp99, avgDp95, avgDMp99, avgDMp95,
			avgEp99, avgEp95, avgLodp99, avgLodp95,
			avgBbp99, avgBbp95)
		normalizedPage.AddCharts(normChart)
	}

	// 7. Write Individual Plots HTML to file
	individualFile := filepath.Join(filepath.Dir(htmlOutFile), "GoBSAseq_IndividualPlots.html")
	fHtml, err := os.Create(individualFile)
	if err != nil {
		return err
	}
	defer fHtml.Close()
	if err := individualPage.Render(fHtml); err != nil {
		return err
	}

	// 8. Write Normalized Overlay HTML to file
	normalizedFile := filepath.Join(filepath.Dir(htmlOutFile), "GoBSAseq_NormalizedOverlay.html")
	fNorm, err := os.Create(normalizedFile)
	if err != nil {
		return err
	}
	defer fNorm.Close()
	if err := normalizedPage.Render(fNorm); err != nil {
		return err
	}

	// 9. Write QTL TSV to file
	fTsv, err := os.Create(qtlOutFile)
	if err != nil {
		return err
	}
	defer fTsv.Close()

	fmt.Fprintf(fTsv, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\n")
	for _, q := range allQTLs {
		fmt.Fprintf(fTsv, "%s\t%d\t%d\t%.6f\t%s\t%s\n", q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI)
	}

	return nil
}

// ============== IMPROVED createNormalizedOverlayChart ==============

// createNormalizedOverlayChart plots all six statistics normalized to their p99 thresholds
// on a single shared y-axis, with p99=1.0 and p95=0.95 clearly marked
// IMPROVED: Better layout, larger charts, enhanced styling, chromosome navigation
func createNormalizedOverlayChart(chrom string, x []int,
	hi, li, dsi, ed, lod, bbl []float64,
	avgHp99, avgHp95, avgLp99, avgLp95 float64,
	avgDp99, avgDp95, avgDMp99, avgDMp95 float64,
	avgEp99, avgEp95, avgLodp99, avgLodp95 float64,
	avgBbp99, avgBbp95 float64) *charts.Line {

	line := charts.NewLine()
	line.SetGlobalOptions(
		charts.WithInitializationOpts(opts.Initialization{
			Theme:  types.ThemeWesteros,
			Width:  "100%",
			Height: "700px",
		}),
		charts.WithTitleOpts(opts.Title{
			Title:    chrom + " Normalized Statistics Overlay",
			Subtitle: "All stats scaled to threshold units (p99 = 1.0, p95 = 0.95). DeltaSI valleys shown as negative. Significance: |value| > 1.0 (p99), |value| > 0.95 (p95)",
			Left:     "center",
			Top:      "2%",
		}),
		charts.WithXAxisOpts(opts.XAxis{
			Name:         "Genomic Position (bp)",
			NameLocation: "middle",
			NameGap:      30,
			AxisLabel: &opts.AxisLabel{
				Rotate: 45,
				Formatter: opts.FuncOpts(`function(value) {
					if (value >= 1000000) {
						return (value / 1000000).toFixed(2) + ' Mb';
					} else if (value >= 1000) {
						return (value / 1000).toFixed(1) + ' kb';
					}
					return value;
				}`),
			},
		}),
		charts.WithYAxisOpts(opts.YAxis{
			Name:         "Normalized Value (threshold units)",
			NameLocation: "middle",
			NameGap:      50,
			Min:          opts.Float(-1.2),
			Max:          opts.Float(1.2),
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
			AxisLabel: &opts.AxisLabel{
				Formatter: opts.FuncOpts(`function(value) {
					if (value === 1.0) return 'p99 (+)';
					if (value === 0.95) return 'p95 (+)';
					if (value === -1.0) return 'p99 (-)';
					if (value === -0.95) return 'p95 (-)';
					if (value === 0) return '0';
					return value.toFixed(2);
				}`),
			},
		}),
		charts.WithDataZoomOpts(
			opts.DataZoom{
				Type:       "slider",
				XAxisIndex: []int{0},
				Start:      0,
				End:        100,
				Height:     30,
				Bottom:     20,
			},
			opts.DataZoom{
				Type:       "inside",
				XAxisIndex: []int{0},
			},
		),
		charts.WithTooltipOpts(opts.Tooltip{
			Show:    opts.Bool(true),
			Trigger: "axis",
			AxisPointer: &opts.AxisPointer{
				Type: "cross",
			},
			Formatter: opts.FuncOpts(`function(params) {
				let result = '<strong>Position: ' + params[0].axisValue + ' bp</strong><br/>';
				params.forEach(function(item) {
					if (item.seriesName.indexOf('threshold') === -1 && item.seriesName !== 'zero') {
						let val = parseFloat(item.value);
						let sig = '';
						if (Math.abs(val) >= 1.0) sig = ' <span style="color:red">★ SIGNIFICANT</span>';
						else if (Math.abs(val) >= 0.95) sig = ' <span style="color:orange">● suggestive</span>';
						result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(3) + sig + '<br/>';
					}
				});
				return result;
			}`),
		}),
		charts.WithLegendOpts(opts.Legend{
			Show:   opts.Bool(true),
			Top:    "8%",
			Left:   "center",
			Type:   "scroll",
			Orient: "horizontal",
		}),
		charts.WithToolboxOpts(opts.Toolbox{
			Show: opts.Bool(true),
			Feature: &opts.ToolBoxFeature{
				SaveAsImage: &opts.ToolBoxFeatureSaveAsImage{
					Show:  opts.Bool(true),
					Title: "Save PNG",
				},
				DataZoom: &opts.ToolBoxFeatureDataZoom{
					Show:  opts.Bool(true),
					Title: map[string]string{"zoom": "Zoom", "back": "Reset Zoom"},
				},
				Restore: &opts.ToolBoxFeatureRestore{
					Show:  opts.Bool(true),
					Title: "Reset",
				},
			},
		}),
		charts.WithGridOpts(opts.Grid{
			Left:         "8%",
			Right:        "5%",
			Top:          "18%",
			Bottom:       "15%",
			ContainLabel: true,
		}),
	)

	// Build normalized data series
	hiNorm := normalizeToThreshold(hi, avgHp99, false)
	liNorm := normalizeToThreshold(li, avgLp99, false)
	dsiNorm := normalizeDeltaSI(dsi, avgDp99, avgDp95, avgDMp99, avgDMp95)
	edNorm := normalizeToThreshold(ed, avgEp99, false)
	lodNorm := normalizeToThreshold(lod, avgLodp99, false)
	bblNorm := normalizeToThreshold(bbl, avgBbp99, false)

	// Build straight lines for threshold references (normalized)
	p99Pos := make([]opts.LineData, len(x))
	p95Pos := make([]opts.LineData, len(x))
	p99Neg := make([]opts.LineData, len(x))
	p95Neg := make([]opts.LineData, len(x))
	zeroLine := make([]opts.LineData, len(x))

	for i := range x {
		p99Pos[i] = opts.LineData{Value: 1.0}
		p95Pos[i] = opts.LineData{Value: 0.95}
		p99Neg[i] = opts.LineData{Value: -1.0}
		p95Neg[i] = opts.LineData{Value: -0.95}
		zeroLine[i] = opts.LineData{Value: 0.0}
	}

	smoothing := true

	// Add threshold reference lines first (so they appear behind data)
	line.SetXAxis(x).
		AddSeries("p99 threshold (+)", p99Pos,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 2, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
			charts.WithZ(1),
		).
		AddSeries("p95 threshold (+)", p95Pos,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.5, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
			charts.WithZ(1),
		).
		AddSeries("zero baseline", zeroLine,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#95a5a6"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
			charts.WithZ(1),
		).
		AddSeries("p95 threshold (-)", p95Neg,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.5, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
			charts.WithZ(1),
		).
		AddSeries("p99 threshold (-)", p99Neg,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 2, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
			charts.WithZ(1),
		)

	// Add normalized data series with distinct, colorblind-friendly colors
	line.AddSeries("HighSI", floatSliceToLineData(hiNorm),
		charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#1f77b4"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: "#1f77b4"}),
		charts.WithZ(2),
	)
	line.AddSeries("LowSI", floatSliceToLineData(liNorm),
		charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#ff7f0e"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: "#ff7f0e"}),
		charts.WithZ(2),
	)
	line.AddSeries("DeltaSI", floatSliceToLineData(dsiNorm),
		charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 3, Color: "#2ca02c"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: "#2ca02c"}),
		charts.WithZ(3),
	)
	line.AddSeries("ED", floatSliceToLineData(edNorm),
		charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#d62728"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: "#d62728"}),
		charts.WithZ(2),
	)
	line.AddSeries("LOD", floatSliceToLineData(lodNorm),
		charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#9467bd"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: "#9467bd"}),
		charts.WithZ(2),
	)
	line.AddSeries("BBLogBF", floatSliceToLineData(bblNorm),
		charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#8c564b"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Color: "#8c564b"}),
		charts.WithZ(2),
	)

	// Add mark areas for significance regions
	line.SetSeriesOptions(
		charts.WithMarkAreaNameTypeItemOpts(
			opts.MarkAreaNameTypeItem{
				Name: "Significant Region (+)",
				Type: "p99 threshold (+)",
			},
			opts.MarkAreaNameTypeItem{
				Name: "Significant Region (+)",
				Type: "p95 threshold (+)",
			},
		),
		charts.WithMarkAreaStyleOpts(
			opts.MarkAreaStyle{
				ItemStyle: &opts.ItemStyle{
					Color: "rgba(231, 76, 60, 0.08)",
				},
			},
		),
	)

	return line
}

// ============== IMPROVED createInteractiveLineChart ==============

// createInteractiveLineChart handles standard and bidirectional (DeltaSI) plots
// IMPROVED: Better dimensions, enhanced tooltips, position formatting
func createInteractiveLineChart(title string, x []int, y []float64, t99, t95, tm99, tm95 float64, hasNegativeThresh bool) *charts.Line {
	line := charts.NewLine()
	line.SetGlobalOptions(
		charts.WithInitializationOpts(opts.Initialization{
			Theme:  types.ThemeWesteros,
			Width:  "100%",
			Height: "500px",
		}),
		charts.WithTitleOpts(opts.Title{
			Title: title,
			Left:  "center",
			Top:   "2%",
		}),
		charts.WithXAxisOpts(opts.XAxis{
			Name:         "Genomic Position (bp)",
			NameLocation: "middle",
			NameGap:      30,
			AxisLabel: &opts.AxisLabel{
				Rotate: 45,
				Formatter: opts.FuncOpts(`function(value) {
					if (value >= 1000000) {
						return (value / 1000000).toFixed(2) + ' Mb';
					} else if (value >= 1000) {
						return (value / 1000).toFixed(1) + ' kb';
					}
					return value;
				}`),
			},
		}),
		charts.WithYAxisOpts(opts.YAxis{
			Name:         title + " Value",
			NameLocation: "middle",
			NameGap:      50,
			SplitLine:    &opts.SplitLine{Show: opts.Bool(true)},
		}),
		charts.WithDataZoomOpts(
			opts.DataZoom{
				Type:       "slider",
				XAxisIndex: []int{0},
				Start:      0,
				End:        100,
				Height:     25,
				Bottom:     15,
			},
			opts.DataZoom{
				Type:       "inside",
				XAxisIndex: []int{0},
			},
		),
		charts.WithTooltipOpts(opts.Tooltip{
			Show:    opts.Bool(true),
			Trigger: "axis",
			AxisPointer: &opts.AxisPointer{
				Type: "cross",
			},
			Formatter: opts.FuncOpts(`function(params) {
				let result = '<strong>Position: ' + params[0].axisValue + ' bp</strong><br/>';
				params.forEach(function(item) {
					let val = parseFloat(item.value);
					if (!isNaN(val)) {
						let sig = '';
						if (item.seriesName === 'Statistic') {
							// Check against thresholds
							let t99 = params.find(p => p.seriesName === 'p99');
							let t95 = params.find(p => p.seriesName === 'p95');
							if (t99 && val > parseFloat(t99.value)) sig = ' <span style="color:red">★ SIGNIFICANT</span>';
							else if (t95 && val > parseFloat(t95.value)) sig = ' <span style="color:orange">● suggestive</span>';
						}
						result += item.marker + ' ' + item.seriesName + ': ' + val.toFixed(4) + sig + '<br/>';
					}
				});
				return result;
			}`),
		}),
		charts.WithLegendOpts(opts.Legend{
			Show: opts.Bool(true),
			Top:  "8%",
			Left: "center",
			Type: "scroll",
		}),
		charts.WithToolboxOpts(opts.Toolbox{
			Show: opts.Bool(true),
			Feature: &opts.ToolBoxFeature{
				SaveAsImage: &opts.ToolBoxFeatureSaveAsImage{
					Show:  opts.Bool(true),
					Title: "Save PNG",
				},
				DataZoom: &opts.ToolBoxFeatureDataZoom{
					Show:  opts.Bool(true),
					Title: map[string]string{"zoom": "Zoom", "back": "Reset Zoom"},
				},
				Restore: &opts.ToolBoxFeatureRestore{
					Show:  opts.Bool(true),
					Title: "Reset",
				},
			},
		}),
		charts.WithGridOpts(opts.Grid{
			Left:         "10%",
			Right:        "5%",
			Top:          "18%",
			Bottom:       "12%",
			ContainLabel: true,
		}),
	)

	// Build straight lines for the average thresholds
	yData := make([]opts.LineData, len(y))
	y99Data := make([]opts.LineData, len(y))
	y95Data := make([]opts.LineData, len(y))

	var ym99Data []opts.LineData
	var ym95Data []opts.LineData
	if hasNegativeThresh {
		ym99Data = make([]opts.LineData, len(y))
		ym95Data = make([]opts.LineData, len(y))
	}

	for i, val := range y {
		yData[i] = opts.LineData{Value: val}
		y99Data[i] = opts.LineData{Value: t99}
		y95Data[i] = opts.LineData{Value: t95}
		if hasNegativeThresh {
			ym99Data[i] = opts.LineData{Value: tm99}
			ym95Data[i] = opts.LineData{Value: tm95}
		}
	}

	smoothing := true
	line.SetXAxis(x).
		AddSeries("Statistic", yData,
			charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5}),
		).
		AddSeries("p99", y99Data,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 2, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
		).
		AddSeries("p95", y95Data,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.5, Color: "#f39c12"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
		)

	if hasNegativeThresh {
		line.AddSeries("p99_valley", ym99Data,
			charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 2, Color: "#e74c3c"}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
		).
			AddSeries("p95_valley", ym95Data,
				charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.5, Color: "#f39c12"}),
				charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0.0)}),
			)
	}

	return line
}
