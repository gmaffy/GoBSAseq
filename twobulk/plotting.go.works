package twobulk

//
//import (
//	"fmt"
//	"image/color"
//	"math"
//	"sort"
//	//"github.com/fatih/color"
//	"gonum.org/v1/plot"
//	"gonum.org/v1/plot/plotter"
//	"gonum.org/v1/plot/vg"
//)
//
//// generatePlots creates publication-quality plots for each chromosome and statistic
//func generatePlots(byChr map[string][]BSAstats, thresholds map[string]Thresholds, config AnalysisConfig) error {
//	// Statistics to plot
//	statsToPlot := []struct {
//		name      string
//		getValue  func(BSAstats) float64
//		getThresh func(Thresholds) float64
//		color     color.Color
//	}{
//		{"DeltaSI",
//			func(s BSAstats) float64 { return s.DeltaSIK },
//			func(t Thresholds) float64 { return t.Delta },
//			color.RGBA{0, 0, 255, 255}}, // Blue
//		{"ED",
//			func(s BSAstats) float64 { return s.EDK },
//			func(t Thresholds) float64 { return t.ED },
//			color.RGBA{255, 0, 0, 255}}, // Red
//		{"LOD",
//			func(s BSAstats) float64 { return s.LODK },
//			func(t Thresholds) float64 { return t.LOD },
//			color.RGBA{0, 255, 0, 255}}, // Green
//		{"Gstat",
//			func(s BSAstats) float64 { return s.Gprime },
//			func(t Thresholds) float64 { return t.G },
//			color.RGBA{255, 255, 0, 255}}, // Yellow
//		{"BBLogBF",
//			func(s BSAstats) float64 { return s.BBLogBFK },
//			func(t Thresholds) float64 { return t.BB },
//			color.RGBA{255, 0, 255, 255}}, // Purple
//	}
//
//	// Create a multi-page plot or individual plots per chromosome
//	for chr, stats := range byChr {
//		if len(stats) == 0 {
//			continue
//		}
//
//		// Create a new plot for each statistic for this chromosome
//		for _, stat := range statsToPlot {
//			p := plot.New()
//			p.Title.Text = fmt.Sprintf("Chromosome %s - %s", chr, stat.name)
//			p.X.Label.Text = "Position (bp)"
//			p.Y.Label.Text = stat.name
//
//			// Create line plot for smoothed statistic
//			pts := make(plotter.XYs, len(stats))
//			for i, s := range stats {
//				pts[i].X = float64(s.POS)
//				pts[i].Y = stat.getValue(s)
//			}
//
//			line, err := plotter.NewLine(pts)
//			if err != nil {
//				return err
//			}
//			line.Color = stat.color
//			line.Width = vg.Points(1.5)
//			p.Add(line)
//
//			// Add threshold line
//			thresh := stat.getThresh(thresholds[chr])
//			if thresh > 0 {
//				threshLine := plotter.NewFunction(func(x float64) float64 { return thresh })
//				threshLine.Color = color.RGBA{255, 0, 0, 255} // Red threshold
//				threshLine.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
//				threshLine.Width = vg.Points(1)
//				p.Add(threshLine)
//
//				// Add negative threshold for DeltaSI if needed
//				if stat.name == "DeltaSI" {
//					negThreshLine := plotter.NewFunction(func(x float64) float64 { return -thresh })
//					negThreshLine.Color = color.RGBA{255, 0, 0, 255}
//					negThreshLine.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
//					negThreshLine.Width = vg.Points(1)
//					p.Add(negThreshLine)
//				}
//			}
//
//			// Customize plot
//			p.Legend.Add(stat.name, line)
//			if thresh > 0 {
//				p.Legend.Add(fmt.Sprintf("Threshold (α=%.3f)", config.Alpha),
//					plotter.NewFunction(func(x float64) float64 { return thresh }))
//			}
//
//			// Save plot
//			filename := fmt.Sprintf("%s_%s_%s.png", config.OutputFile, chr, stat.name)
//			if err := p.Save(12*vg.Inch, 6*vg.Inch, filename); err != nil {
//				fmt.Printf("Error saving plot %s: %v\n", filename, err)
//
//			} else {
//				fmt.Printf("Saved plot: %s\n", filename)
//			}
//		}
//
//		// Create combined plot with all statistics (normalized)
//		combinedPlot := plot.New()
//		combinedPlot.Title.Text = fmt.Sprintf("Chromosome %s - All Statistics (Normalized)", chr)
//		combinedPlot.X.Label.Text = "Position (bp)"
//		combinedPlot.Y.Label.Text = "Normalized Value"
//
//		colors := []color.Color{
//			color.RGBA{0, 0, 255, 255},
//			color.RGBA{255, 0, 0, 255},
//			color.RGBA{0, 255, 0, 255},
//			color.RGBA{255, 255, 0, 255},
//			color.RGBA{255, 0, 255, 255},
//		}
//
//		// Normalize and plot each statistic
//		for i, stat := range statsToPlot {
//			// Find max absolute value for normalization
//			maxVal := 0.0
//			for _, s := range stats {
//				val := stat.getValue(s)
//				if math.Abs(val) > maxVal {
//					maxVal = math.Abs(val)
//				}
//			}
//
//			if maxVal == 0 {
//				continue
//			}
//
//			pts := make(plotter.XYs, len(stats))
//			for j, s := range stats {
//				pts[j].X = float64(s.POS)
//				normalized := stat.getValue(s) / maxVal
//				if math.IsNaN(normalized) {
//					normalized = 0
//				}
//				pts[j].Y = normalized
//			}
//
//			line, err := plotter.NewLine(pts)
//			if err != nil {
//				return err
//			}
//			line.Color = colors[i%len(colors)]
//			line.Width = vg.Points(1.5)
//			combinedPlot.Add(line)
//			combinedPlot.Legend.Add(stat.name, line)
//		}
//
//		combinedFilename := fmt.Sprintf("%s_%s_combined.png", config.OutputFile, chr)
//		if err := combinedPlot.Save(12*vg.Inch, 6*vg.Inch, combinedFilename); err != nil {
//			fmt.Printf("Error saving combined plot: %v\n", err)
//		} else {
//			fmt.Printf("Saved combined plot: %s\n", combinedFilename)
//		}
//	}
//
//	// Create genome-wide Manhattan-style plot
//	if err := generateManhattanPlot(byChr, thresholds, config); err != nil {
//		fmt.Printf("Error generating Manhattan plot: %v\n", err)
//	}
//
//	return nil
//}
//
//// generateManhattanPlot creates a genome-wide Manhattan plot
//func generateManhattanPlot(byChr map[string][]BSAstats, thresholds map[string]Thresholds, config AnalysisConfig) error {
//	// Collect all chromosomes and sort them
//	chroms := make([]string, 0, len(byChr))
//	for chr := range byChr {
//		chroms = append(chroms, chr)
//	}
//	sort.Strings(chroms)
//
//	// Calculate cumulative positions
//	cumulativePos := make(map[string]float64)
//	currentPos := 0.0
//	for _, chr := range chroms {
//		stats := byChr[chr]
//		if len(stats) == 0 {
//			continue
//		}
//		cumulativePos[chr] = currentPos
//		maxPos := stats[len(stats)-1].POS
//		currentPos += float64(maxPos) + 1e6 // Add 1Mb gap between chromosomes
//	}
//
//	// Plot LOD scores (most commonly used in BSA-seq)
//	p := plot.New()
//	p.Title.Text = "Genome-wide Manhattan Plot - LOD Score"
//	p.X.Label.Text = "Chromosome"
//	p.Y.Label.Text = "LOD Score"
//
//	// Prepare data points
//	points := make(plotter.XYs, 0)
//	for _, chr := range chroms {
//		stats := byChr[chr]
//		offset := cumulativePos[chr]
//		for _, s := range stats {
//			points = append(points, plotter.XY{
//				X: offset + float64(s.POS),
//				Y: s.LODK,
//			})
//		}
//	}
//
//	scatter, err := plotter.NewScatter(points)
//	if err != nil {
//		return err
//	}
//	scatter.GlyphStyle.Color = color.RGBA{0, 0, 255, 100}
//	scatter.GlyphStyle.Radius = vg.Points(2)
//	p.Add(scatter)
//
//	// Add threshold lines (use average threshold across chromosomes)
//	avgThresh := 0.0
//	for _, chr := range chroms {
//		avgThresh += thresholds[chr].LOD
//	}
//	avgThresh /= float64(len(chroms))
//
//	threshLine := plotter.NewFunction(func(x float64) float64 { return avgThresh })
//	threshLine.Color = color.RGBA{255, 0, 0, 255}
//	threshLine.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
//	threshLine.Width = vg.Points(1.5)
//	p.Add(threshLine)
//
//	// Customize x-axis to show chromosome boundaries
//	p.X.Tick.Marker = &chromosomeTicker{chroms: chroms, positions: cumulativePos}
//
//	filename := fmt.Sprintf("%s_manhattan.png", config.OutputFile)
//	if err := p.Save(16*vg.Inch, 8*vg.Inch, filename); err != nil {
//		return err
//	}
//
//	fmt.Printf("Saved Manhattan plot: %s\n", filename)
//	return nil
//}
//
//// chromosomeTicker customizes x-axis ticks to show chromosome names
//type chromosomeTicker struct {
//	chroms    []string
//	positions map[string]float64
//}
//
//func (t *chromosomeTicker) Ticks(min, max float64) []plot.Tick {
//	ticks := make([]plot.Tick, 0, len(t.chroms))
//	for _, chr := range t.chroms {
//		pos := t.positions[chr]
//		if pos >= min && pos <= max {
//			ticks = append(ticks, plot.Tick{
//				Value: pos,
//				Label: chr,
//			})
//		}
//	}
//	return ticks
//}
