package stats

import (
	"fmt"
	"math"
	"sort"

	"github.com/xuri/excelize/v2"
)

type MergedQTL struct {
	Chrom          string
	Start          int64
	Stop           int64
	Threshold      string
	CompositeZPos  int64
	CompositeZPeak float64
	BRMPeak        float64
	BRMThreshold   float64
	Source         string
}

func (m MergedQTL) geneSpaceRegion() (string, int, int) {
	return m.Chrom, int(m.Start), int(m.Stop)
}

type PeakIntersection struct {
	Stat           string
	Chrom          string
	Start          float64
	End            float64
	PeakPos        int64
	PeakValue      float64
	Threshold      float64
	ThresholdLevel string
	StartIndex     int
	EndIndex       int
	PeakIndex      int
}

type Peak struct {
	Pos   int64
	Value float64
}

type ConsolidatedQTL struct {
	Chrom     string
	Start     float64
	End       float64
	Threshold string

	DeltaSI *Peak
	Gstat   *Peak
	ED      *Peak
	LOD     *Peak
	BBLogBF *Peak

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

//func FindPeakIntersections(statName string, smoothed []SmoothedStats, thresholds []Thresholds, valueFn func(SmoothedStats) float64, threshFn func(Thresholds) float64, tail Tail) []PeakIntersection {
//
//	if len(smoothed) != len(thresholds) || len(smoothed) < 2 {
//		return nil
//	}
//
//	var threshSum float64
//	for _, t := range thresholds {
//		threshSum += threshFn(t)
//	}
//	avgThresh := threshSum / float64(len(thresholds))
//
//	var peaks []PeakIntersection
//
//	inPeak := false
//
//	var (
//		startPos   float64
//		startIndex int
//		peakPos    int64
//		peakValue  float64
//		peakThresh float64
//		peakIndex  int
//	)
//
//	for i := 0; i < len(smoothed)-1; i++ {
//
//		y1 := valueFn(smoothed[i])
//		y2 := valueFn(smoothed[i+1])
//
//		t1 := avgThresh
//		t2 := avgThresh
//
//		var d1, d2 float64
//
//		switch tail {
//		case UpperTail:
//			d1 = y1 - t1
//			d2 = y2 - t2
//
//		case LowerTail:
//			d1 = t1 - y1
//			d2 = t2 - y2
//		}
//
//		// Enter region
//		if !inPeak && d1 <= 0 && d2 > 0 {
//
//			f := d1 / (d1 - d2)
//
//			startPos = float64(smoothed[i].POS) +
//				f*float64(smoothed[i+1].POS-smoothed[i].POS)
//
//			startIndex = i + 1
//
//			inPeak = true
//			peakPos = smoothed[i+1].POS
//			peakValue = y2
//			peakThresh = t2
//			peakIndex = i + 1
//
//			continue
//		}
//
//		if !inPeak {
//			continue
//		}
//
//		// Update best point
//		better := false
//
//		switch tail {
//		case UpperTail:
//			better = y2 > peakValue
//
//		case LowerTail:
//			better = y2 < peakValue
//		}
//
//		if better {
//			peakValue = y2
//			peakPos = smoothed[i+1].POS
//			peakThresh = t2
//			peakIndex = i + 1
//		}
//
//		// Leave region
//		if d1 > 0 && d2 <= 0 {
//
//			f := d1 / (d1 - d2)
//
//			endPos := float64(smoothed[i].POS) +
//				f*float64(smoothed[i+1].POS-smoothed[i].POS)
//
//			peaks = append(peaks, PeakIntersection{
//				Stat:  statName,
//				Chrom: smoothed[i].CHROM,
//
//				Start: startPos,
//				End:   endPos,
//
//				PeakPos:   peakPos,
//				PeakValue: peakValue,
//				Threshold: peakThresh,
//
//				StartIndex: startIndex,
//				EndIndex:   i,
//				PeakIndex:  peakIndex,
//			})
//
//			inPeak = false
//		}
//	}
//
//	// Region continues to chromosome end
//	if inPeak {
//
//		last := len(smoothed) - 1
//
//		peaks = append(peaks, PeakIntersection{
//			Stat:  statName,
//			Chrom: smoothed[last].CHROM,
//
//			Start: startPos,
//			End:   float64(smoothed[last].POS),
//
//			PeakPos:   peakPos,
//			PeakValue: peakValue,
//			Threshold: peakThresh,
//
//			StartIndex: startIndex,
//			EndIndex:   last,
//			PeakIndex:  peakIndex,
//		})
//	}
//
//	return peaks
//}

func FindPeakIntersections(statName string, smoothed []SmoothedStats, thresholds []Thresholds, valueFn func(SmoothedStats) float64, threshFn func(Thresholds) float64, tail Tail) []PeakIntersection {

	if len(smoothed) != len(thresholds) || len(smoothed) < 2 {
		return nil
	}

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

	// Catch peaks that start at the very beginning of the chromosome
	y0 := valueFn(smoothed[0])
	var d0 float64
	switch tail {
	case UpperTail:
		d0 = y0 - avgThresh
	case LowerTail:
		d0 = avgThresh - y0
	}

	if d0 > 0 {
		inPeak = true
		startPos = float64(smoothed[0].POS)
		startIndex = 0
		peakPos = smoothed[0].POS
		peakValue = y0
		peakThresh = avgThresh
		peakIndex = 0
	}

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

		// Enter region
		if !inPeak && d1 <= 0 && d2 > 0 {

			f := d1 / (d1 - d2)

			startPos = float64(smoothed[i].POS) +
				f*float64(smoothed[i+1].POS-smoothed[i].POS)

			startIndex = i + 1

			inPeak = true
			peakPos = smoothed[i+1].POS
			peakValue = y2
			peakThresh = t2
			peakIndex = i + 1

			continue
		}

		if !inPeak {
			continue
		}

		// Update best point
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

		// Leave region
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

	// Region continues to chromosome end
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

func findPeaksWithFallback(statName string, sm []SmoothedStats, th []Thresholds, valueFn func(SmoothedStats) float64, threshP99 func(Thresholds) float64, threshP95 func(Thresholds) float64, tail Tail, label string) []PeakIntersection {
	chrom := sm[0].CHROM

	tryLevel := func(threshFn func(Thresholds) float64, level string) []PeakIntersection {
		peaks := FindPeakIntersections(statName, sm, th, valueFn, threshFn, tail)
		if len(peaks) == 0 {
			return nil
		}
		for i := range peaks {
			peaks[i].ThresholdLevel = level
			peaks[i].Chrom = chrom
		}
		return peaks
	}

	if peaks := tryLevel(threshP99, "p99"); peaks != nil {
		return peaks
	}
	return tryLevel(threshP95, "p95")
}

type peakStatConfig struct {
	name    string
	valueFn func(SmoothedStats) float64
	p99     func(Thresholds) float64
	p95     func(Thresholds) float64
	tail    Tail
	label   string
}

var twoBulkStatConfigs = []peakStatConfig{
	{"Gstat", func(s SmoothedStats) float64 { return s.SmGstat }, func(t Thresholds) float64 { return t.TwoBulk.GstatP99 }, func(t Thresholds) float64 { return t.TwoBulk.GstatP95 }, UpperTail, "Gstat"},
	{"DeltaSI+", func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP99 }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP95 }, UpperTail, "DeltaSI+"},
	{"DeltaSI-", func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp99 }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp95 }, LowerTail, "DeltaSI-"},
	{"ED4", func(s SmoothedStats) float64 { return s.SmED }, func(t Thresholds) float64 { return t.TwoBulk.ED4P99 }, func(t Thresholds) float64 { return t.TwoBulk.ED4P95 }, UpperTail, "ED4"},
	{"LOD", func(s SmoothedStats) float64 { return s.SmLOD }, func(t Thresholds) float64 { return t.TwoBulk.LODP99 }, func(t Thresholds) float64 { return t.TwoBulk.LODP95 }, UpperTail, "LOD"},
	{"BBLogBF", func(s SmoothedStats) float64 { return s.SmBBLogBF }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP99 }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP95 }, UpperTail, "BBLogBF"},
}

var oneBulkStatConfigs = []peakStatConfig{
	{"AFDev", func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevP99 }, func(t Thresholds) float64 { return t.OneBulk.AFDevP95 }, UpperTail, "AFDev"},
	{"AFDev", func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp99 }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp95 }, LowerTail, "AFDev-"},
	{"OneBulkGstat", func(s SmoothedStats) float64 { return s.SmOneBulkG }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP95 }, UpperTail, "OneBulkGstat"},
	{"OneBulkLOD", func(s SmoothedStats) float64 { return s.SmOneBulkLOD }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP95 }, UpperTail, "OneBulkLOD"},
	{"OneBulkBBLogBF", func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP95 }, UpperTail, "OneBulkBBLogBF"},
}

func collectPeaks(sm []SmoothedStats, th []Thresholds, configs []peakStatConfig) []PeakIntersection {
	var peaks []PeakIntersection
	for _, cfg := range configs {
		if found := findPeaksWithFallback(cfg.name, sm, th, cfg.valueFn, cfg.p99, cfg.p95, cfg.tail, cfg.label); found != nil {
			peaks = append(peaks, found...)
		}
	}
	return peaks
}

func addEvidence(q *ConsolidatedQTL, p PeakIntersection) {
	if p.ThresholdLevel == "p99" {
		q.Threshold = "p99"
	} else if q.Threshold != "p99" && p.ThresholdLevel == "p95" {
		q.Threshold = "p95"
	}

	peak := &Peak{
		Pos:   p.PeakPos,
		Value: p.PeakValue,
	}

	switch p.Stat {

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

func bpOverlap(startA, endA, startB, endB float64) bool {
	a0 := int64(math.Round(math.Min(startA, endA)))
	a1 := int64(math.Round(math.Max(startA, endA)))
	b0 := int64(math.Round(math.Min(startB, endB)))
	b1 := int64(math.Round(math.Max(startB, endB)))
	return b0 <= a1 && b1 >= a0
}

func qtlOverlap(a, b ConsolidatedQTL) bool {
	return a.Chrom == b.Chrom && bpOverlap(a.Start, a.End, b.Start, b.End)
}

func mergeConsolidatedQTL(a, b ConsolidatedQTL) ConsolidatedQTL {
	if b.Start < a.Start {
		a.Start = b.Start
	}
	if b.End > a.End {
		a.End = b.End
	}
	if b.Threshold == "p99" {
		a.Threshold = "p99"
	} else if a.Threshold != "p99" && b.Threshold == "p95" {
		a.Threshold = "p95"
	}

	mergePeak := func(dst **Peak, src *Peak) {
		if *dst == nil {
			*dst = src
		}
	}
	mergePeak(&a.DeltaSI, b.DeltaSI)
	mergePeak(&a.Gstat, b.Gstat)
	mergePeak(&a.ED, b.ED)
	mergePeak(&a.LOD, b.LOD)
	mergePeak(&a.BBLogBF, b.BBLogBF)
	mergePeak(&a.AFDev, b.AFDev)
	mergePeak(&a.OneBulkG, b.OneBulkG)
	mergePeak(&a.OneBulkLOD, b.OneBulkLOD)
	mergePeak(&a.OneBulkBBLogBF, b.OneBulkBBLogBF)

	return a
}

func mergeConsolidatedQTLPass(qtls []ConsolidatedQTL) ([]ConsolidatedQTL, bool) {
	if len(qtls) <= 1 {
		return qtls, false
	}

	sorted := append([]ConsolidatedQTL(nil), qtls...)
	sort.Slice(sorted, func(i, j int) bool {
		if sorted[i].Chrom != sorted[j].Chrom {
			return sorted[i].Chrom < sorted[j].Chrom
		}
		if sorted[i].Start != sorted[j].Start {
			return sorted[i].Start < sorted[j].Start
		}
		return sorted[i].End > sorted[j].End
	})

	out := make([]ConsolidatedQTL, 0, len(sorted))
	current := sorted[0]
	changed := false

	for i := 1; i < len(sorted); i++ {
		next := sorted[i]
		if qtlOverlap(current, next) {
			current = mergeConsolidatedQTL(current, next)
			changed = true
			continue
		}
		out = append(out, current)
		current = next
	}
	out = append(out, current)

	if len(out) < len(qtls) {
		changed = true
	}
	return out, changed
}

func ConsolidateQTLs(peaks []PeakIntersection) []ConsolidatedQTL {

	if len(peaks) == 0 {
		return nil
	}

	sorted := append([]PeakIntersection(nil), peaks...)
	sort.Slice(sorted, func(i, j int) bool {

		if sorted[i].Chrom != sorted[j].Chrom {
			return sorted[i].Chrom < sorted[j].Chrom
		}

		if sorted[i].Start != sorted[j].Start {
			return sorted[i].Start < sorted[j].Start
		}

		return sorted[i].End > sorted[j].End
	})

	var qtls []ConsolidatedQTL

	var current ConsolidatedQTL

	for i, p := range sorted {

		if i == 0 {

			current = ConsolidatedQTL{
				Chrom: p.Chrom,
				Start: p.Start,
				End:   p.End,
			}

			addEvidence(&current, p)
			continue
		}

		if p.Chrom == current.Chrom && bpOverlap(p.Start, p.End, current.Start, current.End) {

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

	for {
		next, changed := mergeConsolidatedQTLPass(qtls)
		qtls = next
		if !changed {
			break
		}
	}

	return qtls
}

func WriteIndividualQTLsToExcel(qtls []ConsolidatedQTL, bsaType string, filename string) error {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

	f := excelize.NewFile()

	header := []interface{}{"Chrom", "Start", "End", "Threshold"}
	if hasBothBulks {
		header = append(header,
			"DeltaSI_Pos", "DeltaSI_Value",
			"Gstat_Pos", "Gstat_Value",
			"ED_Pos", "ED_Value",
			"LOD_Pos", "LOD_Value",
			"BBLogBF_Pos", "BBLogBF_Value",
		)
	}
	if hasOneBulk {
		header = append(header,
			"AFDev_Pos", "AFDev_Value",
			"OneBulkG_Pos", "OneBulkG_Value",
			"OneBulkLOD_Pos", "OneBulkLOD_Value",
			"OneBulkBBLogBF_Pos", "OneBulkBBLogBF_Value",
		)
	}

	f.SetSheetRow("Sheet1", "A1", &header)

	addPeakData := func(rowData *[]interface{}, peak *Peak) {
		if peak != nil {
			*rowData = append(*rowData, peak.Pos, peak.Value)
		} else {
			*rowData = append(*rowData, "", "")
		}
	}

	for rowIdx, qtl := range qtls {
		row := rowIdx + 2

		rowData := []interface{}{
			qtl.Chrom,
			qtl.Start,
			qtl.End,
			qtl.Threshold,
		}

		if hasBothBulks {
			addPeakData(&rowData, qtl.DeltaSI)
			addPeakData(&rowData, qtl.Gstat)
			addPeakData(&rowData, qtl.ED)
			addPeakData(&rowData, qtl.LOD)
			addPeakData(&rowData, qtl.BBLogBF)
		}
		if hasOneBulk {
			addPeakData(&rowData, qtl.AFDev)
			addPeakData(&rowData, qtl.OneBulkG)
			addPeakData(&rowData, qtl.OneBulkLOD)
			addPeakData(&rowData, qtl.OneBulkBBLogBF)
		}

		cell, _ := excelize.CoordinatesToCellName(1, row)
		f.SetSheetRow("Sheet1", cell, &rowData)
	}

	if err := f.SaveAs(filename); err != nil {
		return fmt.Errorf("failed to save Excel file: %w", err)
	}

	return nil
}

func DetectIndividualQTLs(smoothed []SmoothedStats, thresholds []Thresholds, bsaType string) ([]PeakIntersection, []ConsolidatedQTL) {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)

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

		var chrPeaks []PeakIntersection
		if hasBothBulks {
			chrPeaks = append(chrPeaks, collectPeaks(sm, th, twoBulkStatConfigs)...)
		}
		if hasOneBulk {
			chrPeaks = append(chrPeaks, collectPeaks(sm, th, oneBulkStatConfigs)...)
		}

		peaks = append(peaks, chrPeaks...)
		qtls = append(qtls, ConsolidateQTLs(chrPeaks)...)
	}

	return peaks, qtls
}

func DetectCompositeZQTLs(smoothed []SmoothedStats, thresholds []Thresholds, bsaType string) []PeakIntersection {
	_ = bsaType

	var peaks []PeakIntersection
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

		if found := findPeaksWithFallback(
			"CompositeZ",
			sm,
			th,
			func(s SmoothedStats) float64 { return s.CompositeZ },
			func(t Thresholds) float64 { return t.Z.CompositeZP99 },
			func(t Thresholds) float64 { return t.Z.CompositeZP95 },
			UpperTail,
			"CompositeZ",
		); found != nil {
			peaks = append(peaks, found...)
		}
	}

	return peaks
}

func WriteCompositeZQTLsToExcel(peaks []PeakIntersection, filename string) error {
	f := excelize.NewFile()

	header := []interface{}{
		"Chrom", "Start", "End", "Threshold",
		"CompositeZ_Pos", "CompositeZ_Value",
	}
	f.SetSheetRow("Sheet1", "A1", &header)

	for rowIdx, peak := range peaks {
		row := rowIdx + 2
		rowData := []interface{}{
			peak.Chrom,
			peak.Start,
			peak.End,
			peak.ThresholdLevel,
			peak.PeakPos,
			peak.PeakValue,
		}
		cell, _ := excelize.CoordinatesToCellName(1, row)
		f.SetSheetRow("Sheet1", cell, &rowData)
	}

	if err := f.SaveAs(filename); err != nil {
		return fmt.Errorf("failed to save CompositeZ Excel file: %w", err)
	}
	return nil
}

func bestCompositePeakByChrom(peaks []PeakIntersection) map[string]PeakIntersection {
	best := make(map[string]PeakIntersection, len(peaks))
	for _, p := range peaks {
		prev, ok := best[p.Chrom]
		if !ok || p.PeakValue > prev.PeakValue {
			best[p.Chrom] = p
		}
	}
	return best
}

type brmChromSummary struct {
	start, stop int64
	peakPos     int64
	peak        float64
	threshold   float64
}

func summarizeBRMByChrom(blocks []BRMBlock) map[string]brmChromSummary {
	out := make(map[string]brmChromSummary, len(blocks))
	for _, b := range blocks {
		prev, exists := out[b.Chrom]
		if !exists {
			out[b.Chrom] = brmChromSummary{
				start:     b.Start,
				stop:      b.Stop,
				peakPos:   b.PeakPos,
				peak:      b.Peak,
				threshold: b.Threshold,
			}
			continue
		}
		s := prev
		if b.Start < s.start {
			s.start = b.Start
		}
		if b.Stop > s.stop {
			s.stop = b.Stop
		}
		if math.Abs(b.Peak) > math.Abs(s.peak) {
			s.peak = b.Peak
			s.peakPos = b.PeakPos
			s.threshold = b.Threshold
		}
		out[b.Chrom] = s
	}
	return out
}

func MergeCompositeBRM(compositePeaks []PeakIntersection, brmBlocks []BRMBlock) []MergedQTL {
	compositeByChrom := bestCompositePeakByChrom(compositePeaks)
	brmByChrom := summarizeBRMByChrom(brmBlocks)

	chromSet := make(map[string]struct{}, len(compositeByChrom)+len(brmByChrom))
	for chrom := range compositeByChrom {
		chromSet[chrom] = struct{}{}
	}
	for chrom := range brmByChrom {
		chromSet[chrom] = struct{}{}
	}

	chroms := make([]string, 0, len(chromSet))
	for chrom := range chromSet {
		chroms = append(chroms, chrom)
	}
	sort.Strings(chroms)

	nan := math.NaN()
	merged := make([]MergedQTL, 0, len(chroms))

	for _, chrom := range chroms {
		q, hasQ := compositeByChrom[chrom]
		b, hasB := brmByChrom[chrom]

		m := MergedQTL{
			Chrom:        chrom,
			BRMPeak:      nan,
			BRMThreshold: nan,
		}

		switch {
		case hasQ && hasB:
			m.Start = minI64(int64(math.Round(q.Start)), b.start)
			m.Stop = maxI64(int64(math.Round(q.End)), b.stop)
			m.Threshold = q.ThresholdLevel
			m.CompositeZPos = q.PeakPos
			m.CompositeZPeak = q.PeakValue
			m.BRMPeak = b.peak
			m.BRMThreshold = b.threshold
			m.Source = "CompositeZ+BRM"

		case hasQ:
			m.Start = int64(math.Round(q.Start))
			m.Stop = int64(math.Round(q.End))
			m.Threshold = q.ThresholdLevel
			m.CompositeZPos = q.PeakPos
			m.CompositeZPeak = q.PeakValue
			m.Source = "CompositeZ"

		case hasB:
			m.Start = b.start
			m.Stop = b.stop
			m.BRMPeak = b.peak
			m.BRMThreshold = b.threshold
			m.Source = "BRM"
		}

		merged = append(merged, m)
	}

	return merged
}

func WriteFinalQTLsToExcel(merged []MergedQTL, filename string) error {
	f := excelize.NewFile()

	header := []interface{}{
		"Chrom", "Start", "Stop", "Source", "Threshold",
		"CompositeZ_Pos", "CompositeZ_Value",
		"BRM_Peak", "BRM_Threshold",
	}
	f.SetSheetRow("Sheet1", "A1", &header)

	for rowIdx, qtl := range merged {
		row := rowIdx + 2
		rowData := []interface{}{
			qtl.Chrom,
			qtl.Start,
			qtl.Stop,
			qtl.Source,
			qtl.Threshold,
		}

		if qtl.Source == "BRM" {
			rowData = append(rowData, "", "")
		} else {
			rowData = append(rowData, qtl.CompositeZPos, qtl.CompositeZPeak)
		}

		if qtl.Source == "CompositeZ" {
			rowData = append(rowData, "", "")
		} else {
			rowData = append(rowData, nanOrFloat(qtl.BRMPeak), nanOrFloat(qtl.BRMThreshold))
		}

		cell, _ := excelize.CoordinatesToCellName(1, row)
		f.SetSheetRow("Sheet1", cell, &rowData)
	}

	if err := f.SaveAs(filename); err != nil {
		return fmt.Errorf("failed to save final QTL Excel file: %w", err)
	}
	return nil
}

func minI64(a, b int64) int64 {
	if a < b {
		return a
	}
	return b
}

func maxI64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}

func nanOrFloat(v float64) interface{} {
	if math.IsNaN(v) {
		return ""
	}
	return v
}
