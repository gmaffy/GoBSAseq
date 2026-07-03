package stats

import (
	"fmt"
	"sort"
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

			gstatP99Peaks := FindPeakIntersections("Gstat", sm, th, func(s SmoothedStats) float64 { return s.SmGstat }, func(t Thresholds) float64 { return t.TwoBulk.GstatP99 }, UpperTail)

			if len(gstatP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "Gstat", len(gstatP99Peaks))
				fmt.Println(gstatP99Peaks)
				peaks = append(peaks, gstatP99Peaks...)
			} else {
				gstatP95Peaks := FindPeakIntersections("Gstat", sm, th, func(s SmoothedStats) float64 { return s.SmGstat }, func(t Thresholds) float64 { return t.TwoBulk.GstatP95 }, UpperTail)
				if len(gstatP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "Gstat", len(gstatP95Peaks))
					fmt.Println(gstatP95Peaks)
					peaks = append(peaks, gstatP95Peaks...)
				}
			}

			dsiP99Peaks := FindPeakIntersections("DeltaSI+", sm, th, func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP99 }, UpperTail)

			if len(dsiP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "DeltaSI", len(dsiP99Peaks))
				fmt.Println(dsiP99Peaks)
				peaks = append(peaks, dsiP99Peaks...)
			} else {
				dsiP95Peaks := FindPeakIntersections("DeltaSI+", sm, th, func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP95 }, UpperTail)
				if len(dsiP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "DeltaSI+", len(dsiP95Peaks))
					fmt.Println(dsiP95Peaks)
					peaks = append(peaks, dsiP95Peaks...)
				}
			}

			dsiP99Troughs := FindPeakIntersections("DeltaSI-", sm, th, func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp99 }, LowerTail)

			if len(dsiP99Troughs) > 0 {
				fmt.Printf("%-10s %-10s : %2d troughs\n", sm[0].CHROM, "DeltaSI", len(dsiP99Troughs))
				fmt.Println(dsiP99Troughs)
				peaks = append(peaks, dsiP99Troughs...)
			} else {
				dsiP95Troughs := FindPeakIntersections("DeltaSI-", sm, th, func(s SmoothedStats) float64 { return s.SmDeltaSI }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp95 }, LowerTail)
				if len(dsiP95Troughs) > 0 {
					fmt.Printf("%-10s %-10s : %2d troughs\n", sm[0].CHROM, "DeltaSI-", len(dsiP95Troughs))
					fmt.Println(dsiP95Troughs)
					peaks = append(peaks, dsiP95Troughs...)
				}
			}

			ed4P99Peaks := FindPeakIntersections("ED4", sm, th, func(s SmoothedStats) float64 { return s.SmED }, func(t Thresholds) float64 { return t.TwoBulk.ED4P99 }, UpperTail)

			if len(ed4P99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "ED4", len(ed4P99Peaks))
				fmt.Println(ed4P99Peaks)
				peaks = append(peaks, ed4P99Peaks...)
			} else {
				ed4P95Peaks := FindPeakIntersections("ED4", sm, th, func(s SmoothedStats) float64 { return s.SmED }, func(t Thresholds) float64 { return t.TwoBulk.ED4P95 }, UpperTail)
				if len(ed4P95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "ED4", len(ed4P95Peaks))
					fmt.Println(ed4P95Peaks)
					peaks = append(peaks, ed4P95Peaks...)
				}
			}

			loP99Peaks := FindPeakIntersections("LOD", sm, th, func(s SmoothedStats) float64 { return s.SmLOD }, func(t Thresholds) float64 { return t.TwoBulk.LODP99 }, UpperTail)

			if len(loP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "LOD", len(loP99Peaks))
				fmt.Println(loP99Peaks)
				peaks = append(peaks, loP99Peaks...)
			} else {
				loP95Peaks := FindPeakIntersections("LOD", sm, th, func(s SmoothedStats) float64 { return s.SmLOD }, func(t Thresholds) float64 { return t.TwoBulk.LODP95 }, UpperTail)
				if len(loP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "LOD", len(loP95Peaks))
					fmt.Println(loP95Peaks)
					peaks = append(peaks, loP95Peaks...)
				}
			}

			bblogbfP99Peaks := FindPeakIntersections("BBLogBF", sm, th, func(s SmoothedStats) float64 { return s.SmBBLogBF }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP99 }, UpperTail)

			if len(bblogbfP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "BBLogBF", len(bblogbfP99Peaks))
				fmt.Println(bblogbfP99Peaks)
				peaks = append(peaks, bblogbfP99Peaks...)
			} else {
				bblogbfP95Peaks := FindPeakIntersections("BBLogBF", sm, th, func(s SmoothedStats) float64 { return s.SmBBLogBF }, func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP95 }, UpperTail)
				if len(bblogbfP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "BBLogBF", len(bblogbfP95Peaks))
					fmt.Println(bblogbfP95Peaks)
					peaks = append(peaks, bblogbfP95Peaks...)
				}
			}

		}

		if hasOneBulk {

			afDevP99Peaks := FindPeakIntersections("AFDev", sm, th, func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevP99 }, UpperTail)

			if len(afDevP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "AFDev", len(afDevP99Peaks))
				fmt.Println(afDevP99Peaks)
				peaks = append(peaks, afDevP99Peaks...)
			} else {
				afDevP95Peaks := FindPeakIntersections("AFDev", sm, th, func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevP95 }, UpperTail)
				if len(afDevP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "AFDev", len(afDevP95Peaks))
					fmt.Println(afDevP95Peaks)
					peaks = append(peaks, afDevP95Peaks...)
				}
			}

			afDevP99Troughs := FindPeakIntersections("AFDev", sm, th, func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp99 }, LowerTail)

			if len(afDevP99Troughs) > 0 {
				fmt.Printf("%-10s %-10s : %2d troughs\n", sm[0].CHROM, "AFDev", len(afDevP99Troughs))
				fmt.Println(afDevP99Troughs)
				peaks = append(peaks, afDevP99Troughs...)
			} else {
				afDevP95Troughs := FindPeakIntersections("AFDev", sm, th, func(s SmoothedStats) float64 { return s.SmAFDev }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp95 }, LowerTail)
				if len(afDevP95Troughs) > 0 {
					fmt.Printf("%-10s %-10s : %2d troughs\n", sm[0].CHROM, "AFDev", len(afDevP95Troughs))
					fmt.Println(afDevP95Troughs)
					peaks = append(peaks, afDevP95Troughs...)
				}
			}

			oneBulkGstatP99Peaks := FindPeakIntersections("OneBulkGstat", sm, th, func(s SmoothedStats) float64 { return s.SmOneBulkG }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP99 }, UpperTail)

			if len(oneBulkGstatP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "OneBulkGstat", len(oneBulkGstatP99Peaks))
				fmt.Println(oneBulkGstatP99Peaks)
				peaks = append(peaks, oneBulkGstatP99Peaks...)
			} else {
				oneBulkGstatP95Peaks := FindPeakIntersections("OneBulkGstat", sm, th, func(s SmoothedStats) float64 { return s.SmOneBulkG }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP95 }, UpperTail)
				if len(oneBulkGstatP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "OneBulkGstat", len(oneBulkGstatP95Peaks))
					fmt.Println(oneBulkGstatP95Peaks)
					peaks = append(peaks, oneBulkGstatP95Peaks...)
				}
			}

			oneBulkLODP99Peaks := FindPeakIntersections("OneBulkLOD", sm, th, func(s SmoothedStats) float64 { return s.SmOneBulkLOD }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP99 }, UpperTail)

			if len(oneBulkLODP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "OneBulkLOD", len(oneBulkLODP99Peaks))
				fmt.Println(oneBulkLODP99Peaks)
				peaks = append(peaks, oneBulkLODP99Peaks...)
			} else {
				oneBulkLODP95Peaks := FindPeakIntersections("OneBulkLOD", sm, th, func(s SmoothedStats) float64 { return s.SmOneBulkLOD }, func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP95 }, UpperTail)
				if len(oneBulkLODP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "OneBulkLOD", len(oneBulkLODP95Peaks))
					fmt.Println(oneBulkLODP95Peaks)
					peaks = append(peaks, oneBulkLODP95Peaks...)
				}
			}

			oneBulkBBLogBFP99Peaks := FindPeakIntersections("OneBulkBBLogBF", sm, th, func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP99 }, UpperTail)

			if len(oneBulkBBLogBFP99Peaks) > 0 {
				fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "OneBulkBBLogBF", len(oneBulkBBLogBFP99Peaks))
				fmt.Println(oneBulkBBLogBFP99Peaks)
				peaks = append(peaks, oneBulkBBLogBFP99Peaks...)
			} else {
				oneBulkBBLogBFP95Peaks := FindPeakIntersections("OneBulkBBLogBF", sm, th, func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF }, func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP95 }, UpperTail)
				if len(oneBulkBBLogBFP95Peaks) > 0 {
					fmt.Printf("%-10s %-10s : %2d peaks\n", sm[0].CHROM, "OneBulkBBLogBF", len(oneBulkBBLogBFP95Peaks))
					fmt.Println(oneBulkBBLogBFP95Peaks)
					peaks = append(peaks, oneBulkBBLogBFP95Peaks...)
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

	case "ED":
		q.ED = peak

	case "LOD":
		q.LOD = peak

	case "BBLogBF":
		q.BBLogBF = peak

	case "AFDev+", "AFDev-":
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
