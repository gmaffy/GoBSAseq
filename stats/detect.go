package stats

import "fmt"

// PeakIntersection describes one significant region where a smoothed statistic
// exceeds its empirical p99 threshold.
type PeakIntersection struct {
	Stat string

	Chrom string

	// Region boundaries (interpolated crossing points)
	Start float64
	End   float64

	// Peak within the region
	PeakPos   int64
	PeakValue float64
	Threshold float64

	// Useful indices into the original slices
	StartIndex int
	EndIndex   int
	PeakIndex  int
}

func FindPeakIntersections(
	statName string,
	smoothed []SmoothedStats,
	thresholds []Thresholds,
	valueFn func(SmoothedStats) float64,
	threshFn func(Thresholds) float64,
) []PeakIntersection {

	if len(smoothed) != len(thresholds) || len(smoothed) < 2 {
		return nil
	}

	// Calculate averaged threshold for the chromosome (to match plots behavior)
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

		d1 := y1 - t1
		d2 := y2 - t2

		//--------------------------------------------
		// Enter peak
		//--------------------------------------------
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

		//--------------------------------------------
		// Update peak maximum
		//--------------------------------------------
		if y2 > peakValue {
			peakValue = y2
			peakPos = smoothed[i+1].POS
			peakThresh = t2
			peakIndex = i + 1
		}

		//--------------------------------------------
		// Leave peak
		//--------------------------------------------
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

	// Peak continues to chromosome end.
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

func FindAllPeakIntersections(smoothed []SmoothedStats, thresholds []Thresholds, hasBothBulks bool, hasHighBulk bool,
	hasLowBulk bool,
	hasOneBulk bool,
) []PeakIntersection {

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

		if hasBothBulks {

			chrPeaks := FindPeakIntersections(
				"Gstat",
				sm,
				th,
				func(s SmoothedStats) float64 { return s.SmGstat },
				func(t Thresholds) float64 { return t.TwoBulk.GstatP99 },
			)

			fmt.Printf("%-10s %-10s : %2d peaks\n",
				sm[0].CHROM,
				"Gstat",
				len(chrPeaks),
			)

			peaks = append(peaks, chrPeaks...)

			//peaks = append(peaks,
			//	FindPeakIntersections("HighSI", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmHighSI },
			//		func(t Thresholds) float64 { return t.TwoBulk.HighSIP99 })...)
			//
			//peaks = append(peaks,
			//	FindPeakIntersections("LowSI", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmLowSI },
			//		func(t Thresholds) float64 { return t.TwoBulk.LowSIP99 })...)
			//
			//peaks = append(peaks,
			//	FindPeakIntersections("DeltaSI", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmDeltaSI },
			//		func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP99 })...)
			//
			//peaks = append(peaks,
			//	FindPeakIntersections("Gstat", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmGstat },
			//		func(t Thresholds) float64 { return t.TwoBulk.GstatP99 })...)
			//
			//peaks = append(peaks,
			//	FindPeakIntersections("ED", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmED },
			//		func(t Thresholds) float64 { return t.TwoBulk.ED4P99 })...)
			//
			//peaks = append(peaks,
			//	FindPeakIntersections("LOD", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmLOD },
			//		func(t Thresholds) float64 { return t.TwoBulk.LODP99 })...)
			//
			//peaks = append(peaks,
			//	FindPeakIntersections("BBLogBF", sm, th,
			//		func(s SmoothedStats) float64 { return s.SmBBLogBF },
			//		func(t Thresholds) float64 { return t.TwoBulk.BBLogBFP99 })...)
		}

		if hasHighBulk && !hasBothBulks {
			peaks = append(peaks,
				FindPeakIntersections("HighSI", sm, th,
					func(s SmoothedStats) float64 { return s.SmHighSI },
					func(t Thresholds) float64 { return t.TwoBulk.HighSIP99 })...)
		}

		if hasLowBulk && !hasBothBulks {
			peaks = append(peaks,
				FindPeakIntersections("LowSI", sm, th,
					func(s SmoothedStats) float64 { return s.SmLowSI },
					func(t Thresholds) float64 { return t.TwoBulk.LowSIP99 })...)
		}

		if hasOneBulk {

			peaks = append(peaks,
				FindPeakIntersections("AFDev", sm, th,
					func(s SmoothedStats) float64 { return s.SmAFDev },
					func(t Thresholds) float64 { return t.OneBulk.AFDevP99 })...)

			peaks = append(peaks,
				FindPeakIntersections("OneBulkG", sm, th,
					func(s SmoothedStats) float64 { return s.SmOneBulkLOD },
					func(t Thresholds) float64 { return t.OneBulk.OneBulkLODP99 })...)

			peaks = append(peaks,
				FindPeakIntersections("OneBulkBB", sm, th,
					func(s SmoothedStats) float64 { return s.SmOneBulkBBLogBF },
					func(t Thresholds) float64 { return t.OneBulk.OneBulkBBLogBFP99 })...)
		}
	}

	return peaks
}
