package stats

import (
	"cmp"
	"math"
	"slices"
)

// QTL interval merging. Every merge in QTL detection (individual-statistic
// QTLs, CompositeZ regions, BRM blocks, consensus loci) goes through
// mergeRuns.

// mergeRuns sorts items by chromosome, then start (stable), and joins each
// item into the run before it when they are on the same chromosome, it
// starts within gap of the run's end, and canJoin allows it (nil: always).
// join folds b into a and must extend a's end to cover b. One pass over the
// sorted items is complete. The result is sorted by chromosome, then start.
func mergeRuns[T any](items []T, gap int64, span func(*T) (chrom string, start, end float64),
	canJoin func(a, b *T) bool, join func(a *T, b T)) []T {
	out := append([]T(nil), items...)
	slices.SortStableFunc(out, func(a, b T) int {
		ca, sa, _ := span(&a)
		cb, sb, _ := span(&b)
		return cmp.Or(cmp.Compare(ca, cb), cmp.Compare(sa, sb))
	})
	merged := out[:0]
	for _, it := range out {
		if n := len(merged); n > 0 {
			last := &merged[n-1]
			lc, _, le := span(last)
			c, s, _ := span(&it)
			if lc == c && s <= le+float64(gap) && (canJoin == nil || canJoin(last, &it)) {
				join(last, it)
				continue
			}
		}
		merged = append(merged, it)
	}
	return merged
}

// addEvidence adds a region to a QTL: the threshold level (p99 if any
// evidence is p99), and its peak in its statistic's slot if stronger. Regions
// of statistics without a track (not in the track tables) add only the level.
func addEvidence(q *ConsolidatedQTL, p PeakIntersection) {
	if p.ThresholdLevel == "p99" {
		q.Threshold = "p99"
	} else if q.Threshold != "p99" && p.ThresholdLevel == "p95" {
		q.Threshold = "p95"
	}
	if tr, ok := trackByName[p.Stat]; ok {
		keepStrongest(&q.Peaks[tr.slot], &Peak{Pos: p.PeakPos, Value: p.PeakValue, Threshold: p.Threshold, Stat: p.Stat})
	}
}

// keepStrongest sets *dst to p when there is no peak yet or p is stronger
// (larger |value|: lower-tail peaks are negative deviations).
func keepStrongest(dst **Peak, p *Peak) {
	if p != nil && (*dst == nil || math.Abs(p.Value) > math.Abs((*dst).Value)) {
		*dst = p
	}
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

	for s := range a.Peaks {
		keepStrongest(&a.Peaks[s], b.Peaks[s])
	}

	return a
}

// mergePeakIntersections joins one track's regions that overlap or lie
// within gap, keeping the higher peak; the merged region is p99 if any part
// is, and counts all its markers.
func mergePeakIntersections(peaks []PeakIntersection, gap int64) []PeakIntersection {
	return mergeRuns(peaks, gap, func(p *PeakIntersection) (string, float64, float64) { return p.Chrom, p.Start, p.End }, nil,
		func(a *PeakIntersection, b PeakIntersection) {
			a.End = math.Max(a.End, b.End)
			if b.PeakValue > a.PeakValue {
				a.PeakPos, a.PeakValue, a.Threshold = b.PeakPos, b.PeakValue, b.Threshold
			}
			if b.ThresholdLevel == "p99" {
				a.ThresholdLevel = "p99"
			}
			a.NMarkers += b.NMarkers
		})
}

// ConsolidateQTLs joins the individual statistics' regions into QTLs, by
// direction:
//  1. Signed regions (ΔSI or AFDev tails) are joined with regions of the same
//     direction on the same chromosome that overlap or lie within gap.
//  2. Each unsigned region (G′, ED⁴, |ΔSI|, G1, HMM LOD) joins the signed QTL
//     it overlaps most (or, failing overlap, the nearest within gap): it takes
//     that QTL's direction.
//  3. Unsigned regions near no signed QTL are joined among themselves
//     (direction "none").
//
// QTLs of opposite direction are kept apart even where they overlap: the
// high-parent allele is enriched in different bulks. A QTL keeps each
// statistic's strongest peak and is p99 if any evidence is p99. The result
// is sorted by chromosome, then start.
func ConsolidateQTLs(peaks []PeakIntersection, gap int64) []ConsolidatedQTL {
	var signed, unsigned []PeakIntersection
	for _, p := range peaks {
		if p.Direction != "" {
			signed = append(signed, p)
		} else {
			unsigned = append(unsigned, p)
		}
	}
	qtls := mergeByDirection(toQTLs(signed), gap)

	g := float64(gap)
	var rest []PeakIntersection
	for _, p := range unsigned {
		best, bestScore := -1, math.Inf(-1)
		for i := range qtls {
			q := &qtls[i]
			if q.Chrom != p.Chrom || p.Start > q.End+g || q.Start > p.End+g {
				continue
			}
			// overlap length; negative distance when they only meet within gap
			if score := math.Min(p.End, q.End) - math.Max(p.Start, q.Start); score > bestScore {
				best, bestScore = i, score
			}
		}
		if best < 0 {
			rest = append(rest, p)
			continue
		}
		q := &qtls[best]
		q.Start, q.End = math.Min(q.Start, p.Start), math.Max(q.End, p.End)
		addEvidence(q, p)
	}
	// Extending signed QTLs by unsigned evidence can make same-direction
	// QTLs meet: join them again, with the unsigned-only QTLs.
	return mergeByDirection(append(qtls, toQTLs(rest)...), gap)
}

// toQTLs makes one QTL of each region (direction "none" for unsigned ones).
func toQTLs(peaks []PeakIntersection) []ConsolidatedQTL {
	qtls := make([]ConsolidatedQTL, len(peaks))
	for i, p := range peaks {
		dir := p.Direction
		if dir == "" {
			dir = "none"
		}
		qtls[i] = ConsolidatedQTL{Chrom: p.Chrom, Start: p.Start, End: p.End, Direction: dir}
		addEvidence(&qtls[i], p)
	}
	return qtls
}

// mergeByDirection joins QTLs of the same chromosome and direction that
// overlap or lie within gap (opposite directions are never joined, and do
// not block joins between QTLs of one direction), sorted by chromosome, then
// start.
func mergeByDirection(qtls []ConsolidatedQTL, gap int64) []ConsolidatedQTL {
	byDir := map[string][]ConsolidatedQTL{}
	var dirs []string
	for _, q := range qtls {
		if _, ok := byDir[q.Direction]; !ok {
			dirs = append(dirs, q.Direction)
		}
		byDir[q.Direction] = append(byDir[q.Direction], q)
	}
	slices.Sort(dirs)
	var out []ConsolidatedQTL
	for _, d := range dirs {
		out = append(out, mergeRuns(byDir[d], gap, func(q *ConsolidatedQTL) (string, float64, float64) { return q.Chrom, q.Start, q.End }, nil,
			func(a *ConsolidatedQTL, b ConsolidatedQTL) { *a = mergeConsolidatedQTL(*a, b) })...)
	}
	slices.SortStableFunc(out, func(a, b ConsolidatedQTL) int {
		return cmp.Or(cmp.Compare(a.Chrom, b.Chrom), cmp.Compare(a.Start, b.Start), cmp.Compare(a.Direction, b.Direction))
	})
	return out
}
