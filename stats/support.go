package stats

import (
	"math"
	"sort"
)

// SupportZDrop is how far smoothed LOD may fall below its peak before the
// support interval ends, in units of the null LOD scale (NullScales, i.e.
// Z_LOD units). The drop in LOD units is SupportZDrop × that scale, stored in
// Thresholds.SupportDrop.
//
// Calibration (CHANGES.md §14): 216 simulated F2 datasets (depth 30/60/150×,
// bulks of 30/50/100, strong/moderate/weak QTL at random positions, default
// 2 Mb window); 6.5 gives 95% coverage of the true QTL position overall
// (92–98% per depth, bulk size and strength). A fixed drop in raw LOD units
// does not transfer between depths (smoothed LOD grows with depth): the
// classical 1.5 gave 99% coverage but intervals up to whole chromosomes at
// low depth. Null-scaled units varied least across conditions.
const SupportZDrop = 6.5

// SupportInterval localises a QTL. The threshold decides WHETHER a region
// holds a QTL; the region above the threshold can span a whole chromosome
// when the QTL is strong (few meioses ⇒ long linkage blocks), so it says
// little about WHERE the QTL is. The support interval is the stretch around
// the highest smoothed LOD in the region where LOD stays within
// Thresholds.SupportDrop of that peak.
type SupportInterval struct {
	Start   int64   // first marker with LOD ≥ peak − drop
	End     int64   // last marker with LOD ≥ peak − drop
	PeakPos int64   // position of the highest smoothed LOD
	PeakLOD float64 // smoothed LOD at PeakPos
}

// Valid reports whether the interval was computed (positions are ≥ 1).
func (s SupportInterval) Valid() bool { return s.End > 0 }

// union combines the support intervals of two merged regions: the span of
// both, with the peak of the higher one.
func (s SupportInterval) union(o SupportInterval) SupportInterval {
	if !o.Valid() {
		return s
	}
	if !s.Valid() {
		return o
	}
	out := s
	if o.PeakLOD > s.PeakLOD {
		out.PeakPos, out.PeakLOD = o.PeakPos, o.PeakLOD
	}
	out.Start = minI64(s.Start, o.Start)
	out.End = maxI64(s.End, o.End)
	return out
}

// supportLOD returns the smoothed LOD used for support intervals: the
// two-bulk LOD when both bulks exist, otherwise the one-bulk LOD.
func supportLOD(bsaType string) func(SmoothedStats) float64 {
	_, _, hasBoth, _ := BulkFlags(bsaType)
	if hasBoth {
		return func(s SmoothedStats) float64 { return s.SmLOD }
	}
	return func(s SmoothedStats) float64 { return s.SmOneBulkLOD }
}

// supportInChrom computes the support interval for the region [start, end]
// on one chromosome. sm must hold one chromosome's markers sorted by POS.
// The peak is searched inside the region; the interval may extend beyond it
// (it never does in practice, since LOD outside a called region is low).
func supportInChrom(sm []SmoothedStats, lodFn func(SmoothedStats) float64, start, end float64, drop float64) SupportInterval {
	lo := sort.Search(len(sm), func(i int) bool { return float64(sm[i].POS) >= math.Floor(start) })
	hi := sort.Search(len(sm), func(i int) bool { return float64(sm[i].POS) > math.Ceil(end) }) - 1
	if lo >= len(sm) || hi < lo {
		return SupportInterval{}
	}

	peak := lo
	for i := lo + 1; i <= hi; i++ {
		if lodFn(sm[i]) > lodFn(sm[peak]) {
			peak = i
		}
	}
	peakLOD := lodFn(sm[peak])
	if math.IsNaN(peakLOD) {
		return SupportInterval{}
	}
	floor := peakLOD - drop

	l := peak
	for l > 0 && lodFn(sm[l-1]) >= floor {
		l--
	}
	r := peak
	for r < len(sm)-1 && lodFn(sm[r+1]) >= floor {
		r++
	}
	return SupportInterval{Start: sm[l].POS, End: sm[r].POS, PeakPos: sm[peak].POS, PeakLOD: peakLOD}
}

// sortedByPos returns sm sorted by POS (a copy only if it is not already).
func sortedByPos(sm []SmoothedStats) []SmoothedStats {
	if sort.SliceIsSorted(sm, func(i, j int) bool { return sm[i].POS < sm[j].POS }) {
		return sm
	}
	out := append([]SmoothedStats(nil), sm...)
	sort.Slice(out, func(i, j int) bool { return out[i].POS < out[j].POS })
	return out
}
