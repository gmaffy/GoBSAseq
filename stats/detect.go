package stats

import (
	"math"
)

type PeakIntersection struct {
	Stat           string
	Chrom          string
	Start          float64
	End            float64
	PeakPos        int64
	PeakValue      float64
	Threshold      float64
	ThresholdLevel string
	NMarkers       int    // markers beyond the threshold line in the region
	Direction      string // "+" or "−" for a signed statistic's upper or lower tail; "" for unsigned
}

// Peak is a statistic's strongest point in a QTL: position, value, the
// threshold it was called at, and which statistic it is (a slot can hold
// either tail of a signed statistic).
type Peak struct {
	Pos       int64
	Value     float64
	Threshold float64
	Stat      string
}

// qtlSlot is an individual-QTL evidence column: one peak per slot.
type qtlSlot int

const (
	slotDeltaSI qtlSlot = iota // ΔSI upper or lower tail, or |ΔSI|
	slotGstat
	slotED
	slotAFDev // AFDev upper or lower tail
	slotOneBulkG
	slotHMMLOD
	nSlots
)

// slotColumns are the slots' output column prefixes.
var slotColumns = [nSlots]string{"DeltaSI", "Gstat", "ED", "AFDev", "OneBulkG", "HMMLOD"}

// ConsolidatedQTL is an individual-statistic QTL: the joined regions of the
// statistics supporting it.
type ConsolidatedQTL struct {
	Chrom     string
	Start     float64
	End       float64
	Threshold string
	// Direction: "+" or "−" from its signed evidence (ΔSI or AFDev upper or
	// lower tail: the high-parent allele enriched in the high bulk, or
	// depleted); "none" when only unsigned statistics support it.
	Direction string
	// Peaks: each slot's strongest peak; nil where no statistic of the slot
	// supports the QTL.
	Peaks [nSlots]*Peak
	// NStats: the number of statistics supporting the QTL.
	NStats int
	// LeadStat, LeadPos: the lead position, the midpoint of the HMM LOD's
	// 1.5-LOD support interval when the HMM LOD supports the QTL, else the
	// peak of the statistic furthest beyond its threshold (largest
	// |value / threshold|).
	LeadStat string
	LeadPos  int64
}

type Tail int

const (
	UpperTail Tail = iota
	LowerTail
)

// FindPeakIntersections returns the regions of one chromosome's smoothed
// markers (sorted by position) at or beyond the threshold (above it for
// UpperTail, below it for LowerTail), with interpolated ends, each region's
// peak (most extreme value) and its number of markers (scanCrossings).
func FindPeakIntersections(statName string, smoothed []SmoothedStats, threshold float64, valueFn func(SmoothedStats) float64, tail Tail) []PeakIntersection {
	if len(smoothed) == 0 || math.IsNaN(threshold) || math.IsInf(threshold, 0) {
		return nil
	}
	sign := 1.0
	if tail == LowerTail {
		sign = -1
	}
	margin := func(i int) float64 { return sign * (valueFn(smoothed[i]) - threshold) }
	var peaks []PeakIntersection
	for _, c := range scanCrossings(len(smoothed), func(i int) float64 { return float64(smoothed[i].POS) }, margin, margin) {
		pk := &smoothed[c.peak]
		peaks = append(peaks, PeakIntersection{
			Stat: statName, Chrom: pk.CHROM, Start: c.start, End: c.end,
			PeakPos: pk.POS, PeakValue: valueFn(*pk), Threshold: threshold, NMarkers: c.last - c.first + 1,
		})
	}
	return peaks
}

// scanGenome finds one statistic's regions on every chromosome: at its p99
// threshold, or, only if the statistic has no p99 region anywhere in the
// genome, at its p95 threshold (labelled p95, "suggestive"). Deciding the
// fallback genome-wide, not per chromosome, keeps the error rate genome-wide.
// Regions with fewer than minMarkers markers beyond the line are dropped
// (--min-qtl-markers) before the fallback is decided. dir labels the
// regions' direction ("+", "−" or "" for unsigned).
func scanGenome(chroms []chromSlice, statName string, p99, p95 float64, valueFn func(SmoothedStats) float64, tail Tail, minMarkers int, dir string) []PeakIntersection {
	scan := func(threshold float64, level string) []PeakIntersection {
		var out []PeakIntersection
		for _, c := range chroms {
			for _, p := range FindPeakIntersections(statName, c.sm, threshold, valueFn, tail) {
				if p.NMarkers >= minMarkers {
					p.ThresholdLevel, p.Chrom, p.Direction = level, c.chrom, dir
					out = append(out, p)
				}
			}
		}
		return out
	}
	if peaks := scan(p99, "p99"); len(peaks) > 0 {
		return peaks
	}
	return scan(p95, "p95")
}

// statTrack is a scanned statistic: its smoothed value, its genome-wide
// thresholds, the tail scanned, the direction its regions carry ("+" or "−"
// for a signed statistic's tail, "" for unsigned) and its QTL slot. The track
// tables below are the single place a statistic is defined.
type statTrack struct {
	name     string
	value    func(SmoothedStats) float64
	p99, p95 func(Thresholds) float64
	tail     Tail
	dir      string
	slot     qtlSlot
}

// hmmLODTrack scans the pool-HMM LOD (not smoothed) against its genome-wide
// null thresholds.
var hmmLODTrack = statTrack{"HMMLOD", func(s SmoothedStats) float64 { return s.HMMLOD },
	func(t Thresholds) float64 { return t.HMM.LODP99 }, func(t Thresholds) float64 { return t.HMM.LODP95 }, UpperTail, "", slotHMMLOD}

var twoBulkTracks = []statTrack{
	{"Gstat", func(s SmoothedStats) float64 { return s.SmGstat },
		func(t Thresholds) float64 { return t.TwoBulk.GstatP99 }, func(t Thresholds) float64 { return t.TwoBulk.GstatP95 }, UpperTail, "", slotGstat},
	{"DeltaSI+", func(s SmoothedStats) float64 { return s.SmDeltaSI },
		func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP99 }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIP95 }, UpperTail, "+", slotDeltaSI},
	{"DeltaSI-", func(s SmoothedStats) float64 { return s.SmDeltaSI },
		func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp99 }, func(t Thresholds) float64 { return t.TwoBulk.DeltaSIMp95 }, LowerTail, "−", slotDeltaSI},
	{"ED4", func(s SmoothedStats) float64 { return s.SmED },
		func(t Thresholds) float64 { return t.TwoBulk.ED4P99 }, func(t Thresholds) float64 { return t.TwoBulk.ED4P95 }, UpperTail, "", slotED},
	hmmLODTrack,
}

// absDeltaSITracks replaces twoBulkTracks when polarity is unknown (2b,
// --ref-parent none): SmDeltaSI holds smoothed |ΔSI|, so there is one
// unsigned upper-tail ΔSI track, no lower tail and no HMM LOD.
var absDeltaSITracks = func() []statTrack {
	var out []statTrack
	for _, tr := range twoBulkTracks {
		switch tr.name {
		case "DeltaSI-", "HMMLOD":
			continue
		case "DeltaSI+":
			tr.name, tr.dir = "AbsDeltaSI", ""
		}
		out = append(out, tr)
	}
	return out
}()

var oneBulkTracks = []statTrack{
	{"AFDev", func(s SmoothedStats) float64 { return s.SmAFDev },
		func(t Thresholds) float64 { return t.OneBulk.AFDevP99 }, func(t Thresholds) float64 { return t.OneBulk.AFDevP95 }, UpperTail, "+", slotAFDev},
	{"AFDev-", func(s SmoothedStats) float64 { return s.SmAFDev },
		func(t Thresholds) float64 { return t.OneBulk.AFDevMp99 }, func(t Thresholds) float64 { return t.OneBulk.AFDevMp95 }, LowerTail, "−", slotAFDev},
	{"OneBulkGstat", func(s SmoothedStats) float64 { return s.SmOneBulkG },
		func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP99 }, func(t Thresholds) float64 { return t.OneBulk.OneBulkGstatP95 }, UpperTail, "", slotOneBulkG},
	hmmLODTrack,
}

// trackByName finds a statistic's track (and so its slot) by name.
var trackByName = func() map[string]statTrack {
	m := map[string]statTrack{}
	for _, list := range [][]statTrack{twoBulkTracks, absDeltaSITracks, oneBulkTracks} {
		for _, tr := range list {
			m[tr.name] = tr
		}
	}
	return m
}()

// modeSlots are the QTL slots of a mode's statistics, in slot order (the
// individual-QTL output columns).
func modeSlots(bsaType string) []qtlSlot {
	_, _, both, one := BulkFlags(bsaType)
	var used [nSlots]bool
	if both {
		for _, tr := range twoBulkTracks {
			used[tr.slot] = true
		}
	}
	if one {
		for _, tr := range oneBulkTracks {
			used[tr.slot] = true
		}
	}
	var out []qtlSlot
	for s := range nSlots {
		if used[s] {
			out = append(out, s)
		}
	}
	return out
}

// collectPeaks scans every track across the genome.
func collectPeaks(chroms []chromSlice, th Thresholds, tracks []statTrack, minMarkers int) []PeakIntersection {
	var peaks []PeakIntersection
	for _, tr := range tracks {
		peaks = append(peaks, scanGenome(chroms, tr.name, tr.p99(th), tr.p95(th), tr.value, tr.tail, minMarkers, tr.dir)...)
	}
	return peaks
}

// DetectIndividualQTLs finds each statistic's regions (genome-wide p99, or
// p95 when the statistic has no p99 region anywhere; scanGenome) and joins
// them into QTLs by direction (ConsolidateQTLs). signedDSI is
// cfg.SignedDeltaSI(bsaType): false selects the |ΔSI| track set. minMarkers
// (--min-qtl-markers) is the fewest markers beyond the line a region needs.
func DetectIndividualQTLs(smoothed []SmoothedStats, th Thresholds, bsaType string, signedDSI bool, mergeGap int64, minMarkers int) ([]PeakIntersection, []ConsolidatedQTL) {
	_, _, hasBothBulks, hasOneBulk := BulkFlags(bsaType)
	twoBulk := twoBulkTracks
	if !signedDSI {
		twoBulk = absDeltaSITracks
	}
	chroms := byChromosome(smoothed)
	var peaks []PeakIntersection
	if hasBothBulks {
		peaks = append(peaks, collectPeaks(chroms, th, twoBulk, minMarkers)...)
	}
	if hasOneBulk {
		peaks = append(peaks, collectPeaks(chroms, th, oneBulkTracks, minMarkers)...)
	}
	qtls := ConsolidateQTLs(peaks, mergeGap)
	byName := map[string][]SmoothedStats{}
	for _, c := range chroms {
		byName[c.chrom] = c.sm
	}
	for i := range qtls {
		setLead(&qtls[i], byName[qtls[i].Chrom])
	}
	return peaks, qtls
}

// setLead sets a QTL's support count and lead (see ConsolidatedQTL): the
// HMM LOD support-interval midpoint, else the peak furthest beyond its
// threshold. sm is the QTL's chromosome, sorted by position.
func setLead(q *ConsolidatedQTL, sm []SmoothedStats) {
	q.NStats, q.LeadStat, q.LeadPos = 0, "", 0
	best := -1.0
	for _, p := range q.Peaks {
		if p == nil {
			continue
		}
		q.NStats++
		if r := math.Abs(p.Value / p.Threshold); p.Threshold != 0 && r > best {
			best, q.LeadStat, q.LeadPos = r, p.Stat, p.Pos
		}
	}
	if h := q.Peaks[slotHMMLOD]; h != nil && sm != nil {
		s, e := lodSupport(sm, h.Pos, hmmSupportDrop)
		q.LeadStat, q.LeadPos = "HMMLOD", (s+e)/2
	}
}

// DetectCompositeZQTLs finds CompositeZ regions (genome-wide p99, or p95 when
// there is no p99 region anywhere) and joins those within mergeGap.
func DetectCompositeZQTLs(smoothed []SmoothedStats, th Thresholds, mergeGap int64, minMarkers int) []PeakIntersection {
	found := scanGenome(byChromosome(smoothed), "CompositeZ", th.Z.CompositeZP99, th.Z.CompositeZP95,
		func(s SmoothedStats) float64 { return s.CompositeZ }, UpperTail, minMarkers, "")
	return mergePeakIntersections(found, mergeGap)
}
