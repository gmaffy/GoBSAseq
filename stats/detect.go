package stats

import (
	"fmt"
	"sort"

	"github.com/xuri/excelize/v2"
)

type PeakIntersection struct {
	Stat       string
	Chrom      string
	Start      float64
	End        float64
	PeakPos    int64
	PeakValue  float64
	Threshold  float64
	StartIndex int
	EndIndex   int
	PeakIndex  int
}

type Peak struct {
	Pos   int64
	Value float64
}

type ConsolidatedQTL struct {
	Chrom string

	Start float64
	End   float64

	// Highest threshold supporting this QTL (optional)
	Threshold string

	// Two-bulk
	HighSI  *Peak
	LowSI   *Peak
	DeltaSI *Peak
	Gstat   *Peak
	ED      *Peak
	LOD     *Peak
	BBLogBF *Peak

	// One-bulk
	AFDev          *Peak
	OneBulkG       *Peak
	OneBulkLOD     *Peak
	OneBulkBBLogBF *Peak
}
type Tail int

const (
	UpperTail Tail = iota
	LowerTail
)

func FindPeakIntersections(statName string, smoothed []SmoothedStats, thresholds []Thresholds, valueFn func(SmoothedStats) float64, threshFn func(Thresholds) float64, tail Tail) []PeakIntersection {

	if len(smoothed) != len(thresholds) || len(smoothed) < 2 {
		return nil
	}

	// Average threshold across chromosome.
	var threshSum float64
	for _, t := range thresholds {
		threshSum += threshFn(t)
	}
	avgThresh := threshSum / float64(len(thresholds))

	var peaks []PeakIntersection

	inPeak := false

	var (
		startPos   float64
		startIndex int

		peakPos    int64
		peakValue  float64
		peakThresh float64
		peakIndex  int
	)

	for i := 0; i < len(smoothed)-1; i++ {

		y1 := valueFn(smoothed[i])
		y2 := valueFn(smoothed[i+1])

		t1 := avgThresh
		t2 := avgThresh

		var d1, d2 float64

		switch tail {
		case UpperTail:
			d1 = y1 - t1
			d2 = y2 - t2

		case LowerTail:
			d1 = t1 - y1
			d2 = t2 - y2
		}

		//--------------------------------------------------------
		// Enter region
		//--------------------------------------------------------
		if !inPeak && d1 <= 0 && d2 > 0 {

			f := d1 / (d1 - d2)

			startPos = float64(smoothed[i].POS) +
				f*float64(smoothed[i+1].POS-smoothed[i].POS)

			startIndex = i + 1

			inPeak = true

			// Always initialise from the first point inside the region.
			peakPos = smoothed[i+1].POS
			peakValue = y2
			peakThresh = t2
			peakIndex = i + 1

			continue
		}

		if !inPeak {
			continue
		}

		//--------------------------------------------------------
		// Update best point
		//--------------------------------------------------------
		better := false

		switch tail {
		case UpperTail:
			better = y2 > peakValue

		case LowerTail:
			better = y2 < peakValue
		}

		if better {
			peakValue = y2
			peakPos = smoothed[i+1].POS
			peakThresh = t2
			peakIndex = i + 1
		}

		//--------------------------------------------------------
		// Leave region
		//--------------------------------------------------------
		if d1 > 0 && d2 <= 0 {

			f := d1 / (d1 - d2)

			endPos := float64(smoothed[i].POS) +
				f*float64(smoothed[i+1].POS-smoothed[i].POS)

			peaks = append(peaks, PeakIntersection{
				Stat:  statName,
				Chrom: smoothed[i].CHROM,

				Start: startPos,
				End:   endPos,

				PeakPos:   peakPos,
				PeakValue: peakValue,
				Threshold: peakThresh,

				StartIndex: startIndex,
				EndIndex:   i,
				PeakIndex:  peakIndex,
			})

			inPeak = false
		}
	}

	//--------------------------------------------------------
	// Region continues to chromosome end
	//--------------------------------------------------------
	if inPeak {

		last := len(smoothed) - 1

		peaks = append(peaks, PeakIntersection{
			Stat:  statName,
			Chrom: smoothed[last].CHROM,

			Start: startPos,
			End:   float64(smoothed[last].POS),

			PeakPos:   peakPos,
			PeakValue: peakValue,
			Threshold: peakThresh,

			StartIndex: startIndex,
			EndIndex:   last,
			PeakIndex:  peakIndex,
		})
	}

	return peaks
}

// findPeaksWithFallback tries P99 threshold first, falls back to P95 if no peaks found
func findPeaksWithFallback(statName string, sm []SmoothedStats, th []Thresholds, valueFn func(SmoothedStats) float64, threshP99 func(Thresholds) float64, threshP95 func(Thresholds) float64, tail Tail, label string) []PeakIntersection {
	p99Peaks := FindPeakIntersections(statName, sm, th, valueFn, threshP99, tail)
	if len(p99Peaks) > 0 {
		fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, label, len(p99Peaks))
		fmt.Println(p99Peaks)
		return p99Peaks
	}
	
	p95Peaks := FindPeakIntersections(statName, sm, th, valueFn, threshP95, tail)
	if len(p95Peaks) > 0 {
		fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, label, len(p95Peaks))
		fmt.Println(p95Peaks)
		return p95Peaks
	}
	return nil
}

func FindAllPeakIntersections(smoothed []SmoothedStats, thresholds []Thresholds, hasBothBulks bool, hasOneBulk bool) ([]PeakIntersection, []ConsolidatedQTL) {

	var peaks []PeakIntersection
	var qtls []ConsolidatedQTL
	byChr := make(map[string][]int)

	for i := range smoothed {
		byChr[smoothed[i].CHROM] = append(byChr[smoothed[i].CHROM], i)
	}

	for _, idx := range byChr {

		sm := make([]SmoothedStats, len(idx))
		th := make([]Thresholds, len(idx))

		for i, j := range idx {
			sm[i] = smoothed[j]
			th[i] = thresholds[j]
		}

		if hasBothBulks {
			// Two-bulk statistics
			statConfigs := []struct {
				name   string
				valueFn func(SmoothedStats) float64
				p99     func(Thresholds) float64
				p95     func(Thresholds) float64
				tail    Tail
				label   string
			}{
				{"Gstat", func(s SmoothedStats) float64 { return s.SmGstat }, func(t Thresholds) float64 { return t.TwoBulk.GstatP99 }, func(t Thresholds) float64 { return t.TwoBulk.GstatP95 }, UpperTail, "Gstat"},
				{"DeltaSI+", func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP99 }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP95 }, UpperTail, "DeltaSI+"},
				{"DeltaSI-", func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp99 }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp95 }, LowerTail, "DeltaSI-"},
				{"ED4", func(s SmoothedStats) float64 { return s.SmED }, func(t Thresholds) float64 { return t.TwoBulk.ED4P99 }, func(t Thresholds) float64 { return t.TwoBulk.ED4P95 }, UpperTail, "ED4"},
				{"LOD", func(s SmoothedStats) float64 { return s.SmLOD }, func(t Thresholds) float64 { return t.TwoBulk.LODP99 }, func(t Thresholds) float64 { return t.TwoBulk.LODP95 }, UpperTail, "LOD"},
				{"BBLogBF", func(s SmoothedStats) float64 { return s.SmBBLogBF }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP99 }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP95 }, UpperTail, "BBLogBF"},
			}

			for _, cfg := range statConfigs {
				foundPeaks := findPeaksWithFallback(cfg.name, sm, th, cfg.valueFn, cfg.p99, cfg.p95, cfg.tail, cfg.label)
				if foundPeaks != nil {
					peaks = append(peaks, foundPeaks...)
				}
			}
		}

		if hasOneBulk {
			// One-bulk statistics
			oneBulkConfigs := []struct {
				name   string
				valueFn func(SmoothedStats) float64
				p99     func(Thresholds) float64
				p95     func(Thresholds) float64
				tail    Tail
				label   string
			}{
				{"AFDev", func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevP99 }, func(t Thresholds) float64 { return t.OneBulk.AFDevP95 }, UpperTail, "AFDev"},
				{"AFDev", func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp99 }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp95 }, LowerTail, "AFDev-"},
				{"OneBulkGstat", func(s SmoothedStats) float64 { return s.SmOneBulkG }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP95 }, UpperTail, "OneBulkGstat"},
				{"OneBulkLOD", func(s SmoothedStats) float64 { return s.SmOneBulkLOD }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP95 }, UpperTail, "OneBulkLOD"},
				{"OneBulkBBLogBF", func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP95 }, UpperTail, "OneBulkBBLogBF"},
			}

			for _, cfg := range oneBulkConfigs {
				foundPeaks := findPeaksWithFallback(cfg.name, sm, th, cfg.valueFn, cfg.p99, cfg.p95, cfg.tail, cfg.label)
				if foundPeaks != nil {
					peaks = append(peaks, foundPeaks...)
				}
			}
		}
		
		chrQTLs := ConsolidateQTLs(peaks)
		qtls = append(qtls, chrQTLs...)
	}

	return peaks, qtls
}

func overlaps(a, b PeakIntersection) bool {
	return a.Start <= b.End && b.Start <= a.End
}

func addEvidence(q *ConsolidatedQTL, p PeakIntersection) {

	peak := &Peak{
		Pos:   p.PeakPos,
		Value: p.PeakValue,
	}

	switch p.Stat {

	case "HighSI":
		q.HighSI = peak

	case "LowSI":
		q.LowSI = peak

	case "DeltaSI+", "DeltaSI-":
		q.DeltaSI = peak

	case "Gstat":
		q.Gstat = peak

	case "ED", "ED4":
		q.ED = peak

	case "LOD":
		q.LOD = peak

	case "BBLogBF":
		q.BBLogBF = peak

	case "AFDev", "AFDev+", "AFDev-":
		q.AFDev = peak

	case "OneBulkGstat":
		q.OneBulkG = peak

	case "OneBulkLOD":
		q.OneBulkLOD = peak

	case "OneBulkBBLogBF":
		q.OneBulkBBLogBF = peak
	}
}

func ConsolidateQTLs(peaks []PeakIntersection) []ConsolidatedQTL {

	if len(peaks) == 0 {
		return nil
	}

	sort.Slice(peaks, func(i, j int) bool {

		if peaks[i].Chrom != peaks[j].Chrom {
			return peaks[i].Chrom < peaks[j].Chrom
		}

		return peaks[i].Start < peaks[j].Start
	})

	var qtls []ConsolidatedQTL

	var current ConsolidatedQTL

	for i, p := range peaks {

		if i == 0 {

			current = ConsolidatedQTL{
				Chrom: p.Chrom,
				Start: p.Start,
				End:   p.End,
			}

			addEvidence(&current, p)
			continue
		}

		// Same chromosome and overlapping interval?
		if p.Chrom == current.Chrom &&
			p.Start <= current.End &&
			p.End >= current.Start {

			if p.Start < current.Start {
				current.Start = p.Start
			}

			if p.End > current.End {
				current.End = p.End
			}

			addEvidence(&current, p)

			continue
		}

		qtls = append(qtls, current)

		current = ConsolidatedQTL{
			Chrom: p.Chrom,
			Start: p.Start,
			End:   p.End,
		}

		addEvidence(&current, p)
	}

	qtls = append(qtls, current)

	return qtls
}

// WriteConsolidatedQTLToExcel writes the consolidated QTLs to an Excel file
func WriteConsolidatedQTLToExcel(qtls []ConsolidatedQTL, filename string) error {
	f := excelize.NewFile()
	
	// Create header
	header := []interface{}{
		"Chrom", "Start", "End", "Threshold",
		"HighSI_Pos", "HighSI_Value",
		"LowSI_Pos", "LowSI_Value",
		"DeltaSI_Pos", "DeltaSI_Value",
		"Gstat_Pos", "Gstat_Value",
		"ED_Pos", "ED_Value",
		"LOD_Pos", "LOD_Value",
		"BBLogBF_Pos", "BBLogBF_Value",
		"AFDev_Pos", "AFDev_Value",
		"OneBulkG_Pos", "OneBulkG_Value",
		"OneBulkLOD_Pos", "OneBulkLOD_Value",
		"OneBulkBBLogBF_Pos", "OneBulkBBLogBF_Value",
	}
	
	// Write header to row 1
	f.SetSheetRow("Sheet1", "A1", &header)
	
	// Write data rows
	for rowIdx, qtl := range qtls {
		row := rowIdx + 2 // Start from row 2 (after header)
		
		// Build row data
		rowData := []interface{}{
			qtl.Chrom,
			qtl.Start,
			qtl.End,
			qtl.Threshold,
		}
		
		// Helper to add peak data
		addPeakData := func(peak *Peak) {
			if peak != nil {
				rowData = append(rowData, peak.Pos, peak.Value)
			} else {
				rowData = append(rowData, "", "")
			}
		}
		
		// Add all peak data in order
		addPeakData(qtl.HighSI)
		addPeakData(qtl.LowSI)
		addPeakData(qtl.DeltaSI)
		addPeakData(qtl.Gstat)
		addPeakData(qtl.ED)
		addPeakData(qtl.LOD)
		addPeakData(qtl.BBLogBF)
		addPeakData(qtl.AFDev)
		addPeakData(qtl.OneBulkG)
		addPeakData(qtl.OneBulkLOD)
		addPeakData(qtl.OneBulkBBLogBF)
		
		// Write the row
		cell, _ := excelize.CoordinatesToCellName(1, row)
		f.SetSheetRow("Sheet1", cell, &rowData)
	}
	
	// Save the file
	if err := f.SaveAs(filename); err != nil {
		return fmt.Errorf("failed to save Excel file: %w", err)
	}
	
	return nil
}

