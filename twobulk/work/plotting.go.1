package twobulk

import (
	"fmt"
	"os"
	"sort"

	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
)

// QTLRecord holds the detected QTL bounds and peaks
type QTLRecord struct {
	Chrom string
	Start int64
	Stop  int64
	Peak  float64
	Stat  string
	CI    string
}

// GenerateHtmlPlotsAndQTL processes smoothed stats, detects QTLs, and creates interactive HTML charts
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
	page := components.NewPage()
	page.SetLayout(components.PageFlexLayout)

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

		// 5. Create Echarts
		hiChart := createInteractiveLineChart(chrom+" HighSI", x, hi, avgHp99, avgHp95, 0, 0, false)
		liChart := createInteractiveLineChart(chrom+" LowSI", x, li, avgLp99, avgLp95, 0, 0, false)
		dsiChart := createInteractiveLineChart(chrom+" DeltaSI", x, dsi, avgDp99, avgDp95, avgDMp99, avgDMp95, true)
		edChart := createInteractiveLineChart(chrom+" ED", x, ed, avgEp99, avgEp95, 0, 0, false)
		lodChart := createInteractiveLineChart(chrom+" LOD", x, lod, avgLodp99, avgLodp95, 0, 0, false)
		bblChart := createInteractiveLineChart(chrom+" BBLogBF", x, bbl, avgBbp99, avgBbp95, 0, 0, false)

		page.AddCharts(hiChart, liChart, dsiChart, edChart, lodChart, bblChart)
	}

	// 6. Write HTML to file
	fHtml, err := os.Create(htmlOutFile)
	if err != nil {
		return err
	}
	defer fHtml.Close()
	if err := page.Render(fHtml); err != nil {
		return err
	}

	// 7. Write QTL TSV to file
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

// detectQTLs identifies consecutive points crossing the threshold and returns the START, STOP, and PEAK
// detectQTLs identifies consecutive points crossing the threshold and returns the START, STOP, and PEAK
func detectQTLs(chrom string, x []int, y []float64, threshold float64, statName, ci string, isValley bool) []QTLRecord {
	var qtls []QTLRecord
	inQTL := false
	var start, stop int64
	var peak float64

	for i, val := range y {
		condition := val > threshold
		if isValley {
			condition = val < threshold
		}

		if condition {
			if !inQTL {
				inQTL = true
				start = int64(x[i])
				peak = val
			} else {
				if (isValley && val < peak) || (!isValley && val > peak) {
					peak = val
				}
			}
		} else {
			if inQTL {
				stop = int64(x[i-1]) // Record last point above threshold
				qtls = append(qtls, QTLRecord{
					Chrom: chrom,
					Start: start,
					Stop:  stop,
					Peak:  peak,
					Stat:  statName,
					CI:    ci,
				})
				inQTL = false
			}
		}
	}

	// Handle QTL that extends to the very last window of the chromosome
	if inQTL {
		stop = int64(x[len(x)-1])
		qtls = append(qtls, QTLRecord{
			Chrom: chrom,
			Start: start,
			Stop:  stop,
			Peak:  peak,
			Stat:  statName,
			CI:    ci,
		})
	}
	return qtls
}

// createInteractiveLineChart handles standard and bidirectional (DeltaSI) plots[cite: 3]
func createInteractiveLineChart(title string, x []int, y []float64, t99, t95, tm99, tm95 float64, hasNegativeThresh bool) *charts.Line {
	line := charts.NewLine()
	line.SetGlobalOptions(
		charts.WithInitializationOpts(opts.Initialization{Theme: types.ThemeWesteros}),
		charts.WithTitleOpts(opts.Title{Title: title}),
		charts.WithXAxisOpts(opts.XAxis{Name: "Position (bp)"}),
		charts.WithDataZoomOpts(opts.DataZoom{Type: "slider", XAxisIndex: []int{0}}),
		charts.WithTooltipOpts(opts.Tooltip{Show: opts.Bool(true), Trigger: "axis"}),
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
		AddSeries("Statistic", yData, charts.WithLineChartOpts(opts.LineChart{Smooth: &smoothing})).
		AddSeries("p99", y99Data).
		AddSeries("p95", y95Data)

	if hasNegativeThresh {
		line.AddSeries("p99_valley", ym99Data)
		line.AddSeries("p95_valley", ym95Data)
	}

	return line
}
