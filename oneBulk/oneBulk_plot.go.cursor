package oneBulk

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"

	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
	"gonum.org/v1/gonum/stat/distuv"
)

const (
	chartTheme  = types.ThemeWesteros
	chartWidth  = "900px"
	chartHeight = "500px"
	zSig        = 3.0
	zSugg       = 2.0
)

func robustBackground(vals []float64, trimFrac float64) (median, mad float64) {
	if len(vals) == 0 {
		return 0, 0
	}
	sorted := make([]float64, len(vals))
	copy(sorted, vals)
	sort.Float64s(sorted)

	cutIdx := int(math.Round(float64(len(sorted)) * (1.0 - trimFrac)))
	if cutIdx < 1 {
		cutIdx = 1
	}
	trimmed := sorted[:cutIdx]
	median = statQuantile(trimmed, 0.5)

	devs := make([]float64, len(trimmed))
	for i, v := range trimmed {
		devs[i] = math.Abs(v - median)
	}
	sort.Float64s(devs)
	mad = statQuantile(devs, 0.5)
	return median, mad
}

func statQuantile(sorted []float64, p float64) float64 {
	n := len(sorted)
	if n == 0 {
		return 0
	}
	if p <= 0 {
		return sorted[0]
	}
	if p >= 1 {
		return sorted[n-1]
	}
	pos := p * float64(n-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	if lo == hi {
		return sorted[lo]
	}
	frac := pos - float64(lo)
	return sorted[lo]*(1-frac) + sorted[hi]*frac
}

func robustZScore(vals []float64, median, mad float64) []float64 {
	out := make([]float64, len(vals))
	scale := mad * 1.4826
	if scale == 0 {
		return out
	}
	for i, v := range vals {
		out[i] = (v - median) / scale
	}
	return out
}

func commonGlobalOpts(title, subtitle, yLabel string, hasNegativeThresh bool) []charts.GlobalOpts {
	yMin := opts.Float(0)
	if hasNegativeThresh {
		yMin = nil
	}
	posFormatter := `function(v) {
		if (v >= 1e9) return (v/1e9).toFixed(2)+' Gb';
		if (v >= 1e6) return (v/1e6).toFixed(2)+' Mb';
		if (v >= 1e3) return (v/1e3).toFixed(1)+' kb';
		return v;
	}`
	return []charts.GlobalOpts{
		charts.WithInitializationOpts(opts.Initialization{
			Theme:  chartTheme,
			Width:  chartWidth,
			Height: chartHeight,
		}),
		charts.WithTitleOpts(opts.Title{
			Title:    title,
			Subtitle: subtitle,
			Left:     "center",
			Top:      "2%",
		}),
		charts.WithXAxisOpts(opts.XAxis{
			Name:         "Genomic Position",
			NameLocation: "middle",
			NameGap:      35,
			AxisLabel: &opts.AxisLabel{Rotate: 30, Formatter: opts.FuncOpts(posFormatter)},
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
		charts.WithLegendOpts(opts.Legend{Show: opts.Bool(true), Top: "9%", Left: "center", Type: "scroll"}),
		charts.WithToolboxOpts(opts.Toolbox{
			Show:  opts.Bool(true),
			Right: "2%",
			Feature: &opts.ToolBoxFeature{
				SaveAsImage: &opts.ToolBoxFeatureSaveAsImage{Show: opts.Bool(true)},
				DataZoom:    &opts.ToolBoxFeatureDataZoom{Show: opts.Bool(true)},
				Restore:     &opts.ToolBoxFeatureRestore{Show: opts.Bool(true)},
			},
		}),
		charts.WithGridOpts(opts.Grid{Left: "8%", Right: "4%", Top: "20%", Bottom: "14%", ContainLabel: opts.Bool(true)}),
	}
}

func createInteractiveLineChart(title string, x []int64, y []float64, t99, t95, tm99, tm95 float64, hasNegativeThresh bool, brmBlocks []BRMBlock) *charts.Line {
	subtitle := fmt.Sprintf("p99 threshold: %.4f  |  p95 threshold: %.4f  |  shaded: BRM blocks", t99, t95)
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, subtitle, "Value", hasNegativeThresh)...)

	n := len(y)
	yData := make([]opts.LineData, n)
	y99 := make([]opts.LineData, n)
	y95 := make([]opts.LineData, n)
	var ym99, ym95 []opts.LineData
	if hasNegativeThresh {
		ym99 = make([]opts.LineData, n)
		ym95 = make([]opts.LineData, n)
	}
	for i, v := range y {
		yData[i] = opts.LineData{Value: v}
		y99[i] = opts.LineData{Value: t99}
		y95[i] = opts.LineData{Value: t95}
		if hasNegativeThresh {
			ym99[i] = opts.LineData{Value: tm99}
			ym95[i] = opts.LineData{Value: tm95}
		}
	}

	statOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#1f77b4"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	statOpts = append(statOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	line.SetXAxis(positionLabels(x)).
		AddSeries("Statistic", yData, statOpts...).
		AddSeries("p99", y99, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"})).
		AddSeries("p95", y95, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}))

	if hasNegativeThresh {
		line.AddSeries("p99 valley", ym99, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.8, Color: "#e74c3c"})).
			AddSeries("p95 valley", ym95, charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Width: 1.4, Color: "#f39c12"}))
	}
	return line
}

func createRobustZOverlayChartOneBulk(chrom string, x []int64, siZ, absZ, gsZ, edZ, lodZ, bblZ []float64, brmBlocks []BRMBlock) *charts.Line {
	title := chrom + " — Robust Z-score Overlay"
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(title, "Single-bulk robust Z-scores. z=±2 suggestive, z=±3 significant.", "Robust Z-score", true)...)

	n := len(x)
	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}
	zeroOpts := []charts.SeriesOpts{
		charts.WithLineStyleOpts(opts.LineStyle{Type: "solid", Width: 1, Color: "#bdc3c7"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	zeroOpts = append(zeroOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	line.SetXAxis(positionLabels(x)).
		AddSeries("z=0", mkRef(0), zeroOpts...).
		AddSeries("z=+3", mkRef(zSig), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#e74c3c"})).
		AddSeries("z=-3", mkRef(-zSig), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#e74c3c"}))

	type seriesDef struct {
		name string
		data []float64
		col  string
	}
	for _, s := range []seriesDef{
		{"SI", siZ, "#1f77b4"},
		{"AbsSI", absZ, "#ff7f0e"},
		{"Gstat", gsZ, "#17becf"},
		{"ED4", edZ, "#d62728"},
		{"LOD", lodZ, "#9467bd"},
		{"BBLogBF", bblZ, "#8c564b"},
	} {
		line.AddSeries(s.name, floatSliceToLineData(s.data),
			charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
			charts.WithLineStyleOpts(opts.LineStyle{Width: 2, Color: s.col}),
			charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
		)
	}
	return line
}

func createCompositeSignalChartOneBulk(chrom string, x []int64, siZ, absZ, gsZ, edZ, lodZ, bblZ []float64, brmBlocks []BRMBlock) *charts.Line {
	n := len(x)
	composite := make([]float64, n)
	for i := range composite {
		composite[i] = math.Max(math.Abs(siZ[i]),
			math.Max(math.Abs(absZ[i]),
				math.Max(math.Abs(gsZ[i]),
					math.Max(math.Abs(edZ[i]),
						math.Max(math.Abs(lodZ[i]), math.Abs(bblZ[i]))))))
	}
	line := charts.NewLine()
	line.SetGlobalOptions(commonGlobalOpts(chrom+" — Composite Signal (max |Z|)", "Max |Z| across all statistics.", "max |Z-score|", false)...)

	compositeData := floatSliceToLineData(composite)
	compositeOpts := []charts.SeriesOpts{
		charts.WithLineChartOpts(opts.LineChart{Smooth: opts.Bool(true)}),
		charts.WithLineStyleOpts(opts.LineStyle{Width: 2.5, Color: "#2ca02c"}),
		charts.WithItemStyleOpts(opts.ItemStyle{Opacity: opts.Float(0)}),
	}
	compositeOpts = append(compositeOpts, brmBlockMarkAreaOpts(brmBlocks, x)...)

	mkRef := func(val float64) []opts.LineData {
		d := make([]opts.LineData, n)
		for i := range d {
			d[i] = opts.LineData{Value: val}
		}
		return d
	}
	line.SetXAxis(positionLabels(x)).
		AddSeries("Composite", compositeData, compositeOpts...).
		AddSeries("z=2", mkRef(zSugg), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#f39c12"})).
		AddSeries("z=3", mkRef(zSig), charts.WithLineStyleOpts(opts.LineStyle{Type: "dashed", Color: "#e74c3c"}))
	return line
}

func brmBlockMarkAreaOpts(blocks []BRMBlock, x []int64) []charts.SeriesOpts {
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
		areas = append(areas, []opts.MarkAreaData{{XAxis: xLabels[startIdx]}, {XAxis: xLabels[stopIdx]}})
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

func floatSliceToLineData(vals []float64) []opts.LineData {
	ld := make([]opts.LineData, len(vals))
	for i, v := range vals {
		ld[i] = opts.LineData{Value: v}
	}
	return ld
}

func positionLabels(x []int64) []string {
	labels := make([]string, len(x))
	for i, v := range x {
		labels[i] = fmt.Sprintf("%d", v)
	}
	return labels
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
	return f.Close()
}

// GenerateHtmlPlotsAndQTLOneBulk creates HTML plots and QTL tables for single-bulk analysis.
func GenerateHtmlPlotsAndQTLOneBulk(allSmoothed []SmoothedStats, bulkSize int, population string, alphas []float64, htmlOutFile, qtlOutFile string) ([]QTLRecord, error) {
	outDir := filepath.Dir(htmlOutFile)
	const trim = 0.01

	collectStat := func(fn func(SmoothedStats) float64) []float64 {
		v := make([]float64, len(allSmoothed))
		for i, s := range allSmoothed {
			v[i] = fn(s)
		}
		return v
	}
	siMed, siMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.SI }), trim)
	absMed, absMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.AbsSI }), trim)
	gsMed, gsMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.Gstat }), trim)
	edMed, edMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.ED }), trim)
	lodMed, lodMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.LOD }), trim)
	bblMed, bblMAD := robustBackground(collectStat(func(s SmoothedStats) float64 { return s.BBLogBF }), trim)

	byChr := make(map[string][]SmoothedStats)
	for _, s := range allSmoothed {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}
	chroms := make([]string, 0, len(byChr))
	for c := range byChr {
		chroms = append(chroms, c)
	}
	sort.Strings(chroms)

	var allQTLs, allConsensusQTLs, allMaxZQTLs []QTLRecord
	var allBRMBlocks []BRMBlock

	popLevel := 0
	if population == "F2" {
		popLevel = 1
	}
	brmAlpha := defaultBRMAlpha
	for _, alpha := range alphas {
		if alpha > brmAlpha && alpha < 1 {
			brmAlpha = alpha
		}
	}
	brmUAlpha := distuv.UnitNormal.Quantile(1 - brmAlpha/2)

	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Statistics (Single Bulk)"

	robustZPage := components.NewPage()
	robustZPage.SetLayout(components.PageFlexLayout)
	robustZPage.PageTitle = "GoBSAseq — Robust Z-score Overlay (Single Bulk)"

	compositePage := components.NewPage()
	compositePage.SetLayout(components.PageFlexLayout)
	compositePage.PageTitle = "GoBSAseq — Composite Signal (Single Bulk)"

	for _, chrom := range chroms {
		stats := byChr[chrom]
		if len(stats) == 0 {
			continue
		}

		nf := float64(len(stats))
		var sumSp99, sumSp95, sumSMp99, sumSMp95, sumAp99, sumAp95, sumGs99, sumGs95, sumEp99, sumEp95, sumLod99, sumLod95, sumBb99, sumBb95 float64
		for _, s := range stats {
			t := s.thresholds
			sumSp99 += t.SiP99
			sumSp95 += t.SiP95
			sumSMp99 += t.SiMp99
			sumSMp95 += t.SiMp95
			sumAp99 += t.AbsP99
			sumAp95 += t.AbsP95
			sumGs99 += t.GsP99
			sumGs95 += t.GsP95
			sumEp99 += t.EdP99
			sumEp95 += t.EdP95
			sumLod99 += t.LodP99
			sumLod95 += t.LodP95
			sumBb99 += t.BbP99
			sumBb95 += t.BbP95
		}

		n := len(stats)
		x := make([]int64, n)
		si := make([]float64, n)
		abs := make([]float64, n)
		gs := make([]float64, n)
		ed := make([]float64, n)
		lod := make([]float64, n)
		bbl := make([]float64, n)
		siT99, siTM99 := make([]float64, n), make([]float64, n)
		siT95, siTM95 := make([]float64, n), make([]float64, n)
		absT99, absT95 := make([]float64, n), make([]float64, n)
		gsT99, gsT95 := make([]float64, n), make([]float64, n)
		edT99, edT95 := make([]float64, n), make([]float64, n)
		lodT99, lodT95 := make([]float64, n), make([]float64, n)
		bblT99, bblT95 := make([]float64, n), make([]float64, n)

		for i, s := range stats {
			x[i] = s.POS
			si[i], abs[i], gs[i], ed[i], lod[i], bbl[i] = s.SI, s.AbsSI, s.Gstat, s.ED, s.LOD, s.BBLogBF
			t := s.thresholds
			siT99[i], siTM99[i] = t.SiP99, t.SiMp99
			siT95[i], siTM95[i] = t.SiP95, t.SiMp95
			absT99[i], absT95[i] = t.AbsP99, t.AbsP95
			gsT99[i], gsT95[i] = t.GsP99, t.GsP95
			edT99[i], edT95[i] = t.EdP99, t.EdP95
			lodT99[i], lodT95[i] = t.LodP99, t.LodP95
			bblT99[i], bblT95[i] = t.BbP99, t.BbP95
		}

		var chromQTLs []QTLRecord
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siT99, "SI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siTM99, "SI", "99", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siT95, "SI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, si, siTM95, "SI", "95", true, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, abs, absT99, "AbsSI", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, abs, absT95, "AbsSI", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, gs, gsT99, "Gstat", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, gs, gsT95, "Gstat", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, ed, edT99, "ED4", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, ed, edT95, "ED4", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, lod, lodT99, "LOD", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, lod, lodT95, "LOD", "95", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, bbl, bblT99, "BBLogBF", "99", false, "Permutation")...)
		chromQTLs = append(chromQTLs, detectQTLsAdaptive(chrom, x, bbl, bblT95, "BBLogBF", "95", false, "Permutation")...)

		siZ := robustZScore(si, siMed, siMAD)
		absZ := robustZScore(abs, absMed, absMAD)
		gsZ := robustZScore(gs, gsMed, gsMAD)
		edZ := robustZScore(ed, edMed, edMAD)
		lodZ := robustZScore(lod, lodMed, lodMAD)
		bblZ := robustZScore(bbl, bblMed, bblMAD)

		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, siZ, zSig, "SI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, absZ, zSig, "AbsSI_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, gsZ, zSig, "Gstat_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, edZ, zSig, "ED4_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, lodZ, zSig, "LOD_Z", "z3", false, "ZScore")...)
		chromQTLs = append(chromQTLs, detectQTLs(chrom, x, bblZ, zSig, "BBLogBF_Z", "z3", false, "ZScore")...)

		cQTLs := detectConsensusQTLsOneBulk(chrom, stats)
		allConsensusQTLs = append(allConsensusQTLs, cQTLs...)

		chromBRMBlocks := calculateBRMBlocksOneBulk(chrom, stats, bulkSize, popLevel, brmUAlpha)
		allBRMBlocks = append(allBRMBlocks, chromBRMBlocks...)

		composite := make([]float64, n)
		for i := range composite {
			composite[i] = math.Max(math.Abs(siZ[i]),
				math.Max(math.Abs(absZ[i]),
					math.Max(math.Abs(gsZ[i]),
						math.Max(math.Abs(edZ[i]),
							math.Max(math.Abs(lodZ[i]), math.Abs(bblZ[i]))))))
		}
		maxZQTLs := detectQTLs(chrom, x, composite, zSig, "Composite_Z", "z3", false, "MaxZ")
		allMaxZQTLs = append(allMaxZQTLs, maxZQTLs...)
		allMaxZQTLs = append(allMaxZQTLs, intersectQTLsWithBRM(maxZQTLs, chromBRMBlocks, "CompositeHighConfidence")...)

		chromQTLs = append(chromQTLs, intersectQTLsWithBRM(chromQTLs, chromBRMBlocks, "HighConfidence")...)
		allQTLs = append(allQTLs, chromQTLs...)

		robustZPage.AddCharts(createRobustZOverlayChartOneBulk(chrom, x, siZ, absZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))
		individualPage.AddCharts(
			createInteractiveLineChart(chrom+" SI", x, si, sumSp99/nf, sumSp95/nf, sumSMp99/nf, sumSMp95/nf, true, chromBRMBlocks),
			createInteractiveLineChart(chrom+" AbsSI", x, abs, sumAp99/nf, sumAp95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" Gstat", x, gs, sumGs99/nf, sumGs95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" ED4", x, ed, sumEp99/nf, sumEp95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" LOD", x, lod, sumLod99/nf, sumLod95/nf, 0, 0, false, chromBRMBlocks),
			createInteractiveLineChart(chrom+" BBLogBF", x, bbl, sumBb99/nf, sumBb95/nf, 0, 0, false, chromBRMBlocks),
		)
		compositePage.AddCharts(createCompositeSignalChartOneBulk(chrom, x, siZ, absZ, gsZ, edZ, lodZ, bblZ, chromBRMBlocks))
	}

	if err := writeHTMLPage(individualPage, filepath.Join(outDir, "GoBSAseq_IndividualPlots.html")); err != nil {
		return nil, err
	}
	if err := writeHTMLPage(robustZPage, htmlOutFile); err != nil {
		return nil, err
	}
	if err := writeHTMLPage(compositePage, filepath.Join(outDir, "GoBSAseq_CompositeSignal.html")); err != nil {
		return nil, err
	}

	fTsv, err := os.Create(qtlOutFile)
	if err != nil {
		return nil, fmt.Errorf("create qtl file: %w", err)
	}
	fmt.Fprintf(fTsv, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\tSOURCE\n")
	for _, q := range allQTLs {
		fmt.Fprintf(fTsv, "%s\t%d\t%d\t%.6f\t%s\t%s\t%s\n", q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI, q.Source)
	}
	_ = fTsv.Close()

	fCons, err := os.Create(filepath.Join(outDir, "GoBSAseq_QTL_CONSENSUS.tsv"))
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(fCons, "CHROM\tSTART\tSTOP\t#STATS\tSTAT\tCI\tSTATS\n")
	for _, q := range allConsensusQTLs {
		fmt.Fprintf(fCons, "%s\t%d\t%d\t%d\t%s\t%s\t%s\n", q.Chrom, q.Start, q.Stop, int(q.Peak), q.Stat, q.CI, q.Source)
	}
	_ = fCons.Close()

	fMaxZ, err := os.Create(filepath.Join(outDir, "GoBSAseq_QTL_MAX_Z.tsv"))
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(fMaxZ, "CHROM\tSTART\tSTOP\tPEAK\tSTAT\tCI\tSOURCE\n")
	for _, q := range allMaxZQTLs {
		fmt.Fprintf(fMaxZ, "%s\t%d\t%d\t%.6f\t%s\t%s\t%s\n", q.Chrom, q.Start, q.Stop, q.Peak, q.Stat, q.CI, q.Source)
	}
	_ = fMaxZ.Close()

	fBRM, err := os.Create(filepath.Join(outDir, "GoBSAseq_BRMBlocks.tsv"))
	if err != nil {
		return nil, err
	}
	fmt.Fprintf(fBRM, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK_ABS_SI\tBRM_THRESHOLD\n")
	for _, b := range allBRMBlocks {
		fmt.Fprintf(fBRM, "%s\t%d\t%d\t%d\t%.6f\t%.6f\n", b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.Threshold)
	}
	_ = fBRM.Close()

	return allMaxZQTLs, nil
}
