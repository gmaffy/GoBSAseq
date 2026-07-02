package stats

func detectQtlPeaks(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX int, rightIntersectX int, found bool) {

	abovePoints := make([]int, 0)
	for i := 0; i < len(yData); i++ {
		if yData[i] > y9Data[i] {
			abovePoints = append(abovePoints, i)
		}
	}

	if len(abovePoints) == 0 {
		return 0, 0, 0, 0, false
	}

	maxY := yData[abovePoints[0]]
	peakIndex := abovePoints[0]
	for _, idx := range abovePoints {
		if yData[idx] > maxY {
			maxY = yData[idx]
			peakIndex = idx
		}
	}
	peakX = x[peakIndex]
	peakY = maxY

	var intersections []int
	for i := 1; i < len(yData); i++ {

		if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) || (yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
			intersections = append(intersections, i)
		}
	}

	if len(intersections) < 2 {
		return peakX, peakY, 0, 0, false
	}
	leftIntersect := -1
	rightIntersect := -1

	for _, idx := range intersections {
		if idx < peakIndex && (leftIntersect == -1 || idx > leftIntersect) {
			leftIntersect = idx
		} else if idx > peakIndex && (rightIntersect == -1 || idx < rightIntersect) {
			rightIntersect = idx
		}
	}

	if leftIntersect == -1 || rightIntersect == -1 {
		return peakX, peakY, 0, 0, false
	}

	return peakX, peakY, x[leftIntersect], x[rightIntersect], true
}

func detectQtlValleys(x []int, yData, y9Data []float64) (peakX int, peakY float64, leftIntersectX int, rightIntersectX int, found bool) {

	belowPoints := make([]int, 0)
	for i := 0; i < len(yData); i++ {
		if yData[i] < y9Data[i] {
			belowPoints = append(belowPoints, i)
		}
	}

	if len(belowPoints) == 0 {
		return 0, 0, 0, 0, false
	}

	minY := yData[belowPoints[0]]
	peakIndex := belowPoints[0]
	for _, idx := range belowPoints {
		if yData[idx] < minY {
			minY = yData[idx]
			peakIndex = idx
		}
	}
	peakX = x[peakIndex]
	peakY = minY

	var intersections []int
	for i := 1; i < len(yData); i++ {
		if (yData[i-1] < y9Data[i-1] && yData[i] > y9Data[i]) ||
			(yData[i-1] > y9Data[i-1] && yData[i] < y9Data[i]) {
			intersections = append(intersections, i)
		}
	}

	if len(intersections) < 2 {
		return peakX, peakY, 0, 0, false
	}

	leftIntersect := -1
	rightIntersect := -1

	for _, idx := range intersections {
		if idx < peakIndex && (leftIntersect == -1 || idx > leftIntersect) {
			leftIntersect = idx
		} else if idx > peakIndex && (rightIntersect == -1 || idx < rightIntersect) {
			rightIntersect = idx
		}
	}

	if leftIntersect == -1 || rightIntersect == -1 {
		return peakX, peakY, 0, 0, false
	}

	return peakX, peakY, x[leftIntersect], x[rightIntersect], true
}
