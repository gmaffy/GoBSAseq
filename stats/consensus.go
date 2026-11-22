package stats

import (
	"cmp"
	"math"
	"slices"
	"strings"
)

// Final QTLs: consensus of the HMM LOD, BRM and DsiGconsensusZ
// (qtl_detection_todo.md, 1.6).
//
// Each method contributes its significant regions on a chromosome: HMM LOD
// and DsiGconsensusZ regions beyond their genome-wide p99 thresholds, and BRM
// blocks at the genome-wide --brm-alpha (default 0.01, the same level);
// regions need --min-qtl-markers markers. Regions of any method that overlap,
// or lie within the merge distance of each other, form a locus. A locus is a
// consensus (final) QTL when every available method has a region in it: all
// three, or BRM and DsiGconsensusZ in 2b without --ref-parent (no HMM LOD).
// Loci with partial support are reported too, labelled by their support, but
// are not final.
//
// BRM and DsiGconsensusZ both test the bulks' allele-frequency contrast (BRM's
// z is DsiGconsensusZ's ΔSI component), so they are one line of evidence, the
// "ΔSI family"; the HMM LOD, which models linkage along the chromosome
// instead of averaging over a window, is the other. Evidence names the lines
// present (FinalQTL.Evidence).
//
// A locus reports its span (the union of its regions) and its core: where
// all its methods are above their thresholds at once (the intersection; the
// piece holding the lead peak, else the longest). If the methods' regions
// only meet within the merge distance, without overlapping, the core is the
// lead method's region. The lead position is the midpoint of the HMM LOD's
// 1.5-LOD support interval (the markers around its maximum within 1.5 LOD of
// it; the standard QTL support interval, and robust where the LOD is a flat
// plateau and its maximum is arbitrary), else the BRM peak, else the
// DsiGconsensusZ peak.

// Method names, in report order.
const (
	methodHMM  = "HMMLOD"
	methodBRM  = "BRM"
	methodDsiG = "DsiGconsensusZ"
)

// Evidence lines (FinalQTL.Evidence).
const (
	evidenceHMM = "HMM"
	evidenceDSI = "ΔSI-family"
)

// MethodPeak is one method's strongest peak in a locus. Z is the
// standardised value where the method has one (BRM's signal/σ), else NaN.
type MethodPeak struct {
	Pos              int64
	Value, Threshold float64
	Z                float64
}

// FinalQTL is a locus of the final QTL list.
type FinalQTL struct {
	Chrom string
	// Start, Stop: the core, where all supporting methods are above threshold.
	Start, Stop int64
	// SpanStart, SpanStop: the union of the supporting regions.
	SpanStart, SpanStop int64
	// Consensus: every available method supports the locus.
	Consensus bool
	// Support: the supporting methods, e.g. "HMMLOD+BRM+DsiGconsensusZ".
	Support string
	// Evidence: the independent lines of evidence among them: "HMM" (the HMM
	// LOD) and "ΔSI-family" (BRM and/or DsiGconsensusZ), e.g.
	// "HMM+ΔSI-family". In 2b without --ref-parent a consensus rests on the
	// ΔSI family alone.
	Evidence string
	// Direction: "+" or "−" from the BRM peak (ΔSI or AFDev sign: the
	// high-parent allele enriched in the high bulk, or depleted); "none" when
	// BRM is polarity-free or absent; "mixed" when BRM blocks of both signs
	// fall in the locus.
	Direction  string
	LeadMethod string
	LeadPos    int64
	// HMMSupportStart/Stop: the HMM LOD's 1.5-LOD support interval (0 when
	// there is no HMM region).
	HMMSupportStart, HMMSupportStop int64
	// Per-method strongest peaks; nil when the method has no region here.
	HMM, BRM, DsiG *MethodPeak
	// Statistics: the individual statistics whose regions overlap the span.
	Statistics string
}

func (q FinalQTL) geneSpaceRegion() (string, int, int) { return q.Chrom, int(q.Start), int(q.Stop) }

// consensusRegion is one method's region in the consensus.
type consensusRegion struct {
	method     string
	start, end float64
	peak       MethodPeak
}

// ConsensusQTLs builds the final QTL list from the smoothed markers (HMM LOD
// and DsiGconsensusZ regions are found here at p99), the BRM blocks and the
// individual-statistic regions (for the Statistics column). signedDSI is
// cfg.SignedDeltaSI(bsaType); hmmOn whether the HMM LOD was computed.
func ConsensusQTLs(smoothed []SmoothedStats, th Thresholds, brm []BRMBlock, individual []PeakIntersection,
	signedDSI, hmmOn bool, gap int64, minMarkers int) []FinalQTL {
	methods := []string{methodBRM, methodDsiG}
	if hmmOn {
		methods = []string{methodHMM, methodBRM, methodDsiG}
	}
	brmBy := map[string][]BRMBlock{}
	for _, b := range brm {
		brmBy[b.Chrom] = append(brmBy[b.Chrom], b)
	}
	indBy := map[string][]PeakIntersection{}
	for _, p := range individual {
		indBy[p.Chrom] = append(indBy[p.Chrom], p)
	}
	nan := math.NaN()
	// regions scans one track at its p99 threshold.
	regions := func(sm []SmoothedStats, method string, threshold float64, get func(SmoothedStats) float64) []consensusRegion {
		var out []consensusRegion
		var found []PeakIntersection
		for _, p := range FindPeakIntersections(method, sm, threshold, get, UpperTail) {
			if p.NMarkers >= minMarkers {
				found = append(found, p)
			}
		}
		for _, p := range mergePeakIntersections(found, gap) {
			out = append(out, consensusRegion{method, p.Start, p.End, MethodPeak{p.PeakPos, p.PeakValue, threshold, nan}})
		}
		return out
	}

	var final []FinalQTL
	for _, c := range byChromosome(smoothed) {
		var rs []consensusRegion
		if hmmOn && th.HMM.LODP99 > 0 {
			rs = append(rs, regions(c.sm, methodHMM, th.HMM.LODP99, func(s SmoothedStats) float64 { return s.HMMLOD })...)
		}
		rs = append(rs, regions(c.sm, methodDsiG, th.Consensus.P99, func(s SmoothedStats) float64 { return s.DsiGconsensusZ })...)
		for _, b := range brmBy[c.chrom] {
			rs = append(rs, consensusRegion{methodBRM, float64(b.Start), float64(b.Stop), MethodPeak{b.PeakPos, b.Peak, b.Threshold, b.PeakZ}})
		}
		// Loci: regions connected by overlap or within gap.
		type locus struct {
			start, end float64
			rs         []consensusRegion
		}
		loci := make([]locus, len(rs))
		for i, r := range rs {
			loci[i] = locus{r.start, r.end, []consensusRegion{r}}
		}
		loci = mergeRuns(loci, gap, func(l *locus) (string, float64, float64) { return "", l.start, l.end }, nil,
			func(a *locus, b locus) { a.end, a.rs = math.Max(a.end, b.end), append(a.rs, b.rs...) })
		for _, l := range loci {
			final = append(final, buildLocus(c.chrom, l.rs, methods, signedDSI, indBy[c.chrom], c.sm))
		}
	}
	return final
}

// buildLocus summarises one locus's regions.
func buildLocus(chrom string, rs []consensusRegion, methods []string, signedDSI bool, individual []PeakIntersection, sm []SmoothedStats) FinalQTL {
	q := FinalQTL{Chrom: chrom, SpanStart: math.MaxInt64, SpanStop: math.MinInt64}
	byMethod := map[string][]consensusRegion{}
	for _, r := range rs {
		byMethod[r.method] = append(byMethod[r.method], r)
		q.SpanStart = min(q.SpanStart, int64(math.Round(r.start)))
		q.SpanStop = max(q.SpanStop, int64(math.Round(r.end)))
	}
	strongest := func(m string) *MethodPeak {
		var best *MethodPeak
		for i := range byMethod[m] {
			p := byMethod[m][i].peak
			if best == nil || math.Abs(p.Value) > math.Abs(best.Value) {
				best = &p
			}
		}
		return best
	}
	q.HMM, q.BRM, q.DsiG = strongest(methodHMM), strongest(methodBRM), strongest(methodDsiG)

	var support []string
	q.Consensus = true
	for _, m := range methods {
		if len(byMethod[m]) > 0 {
			support = append(support, m)
		} else {
			q.Consensus = false
		}
	}
	q.Support = strings.Join(support, "+")
	var lines []string
	if q.HMM != nil {
		lines = append(lines, evidenceHMM)
	}
	if q.BRM != nil || q.DsiG != nil {
		lines = append(lines, evidenceDSI)
	}
	q.Evidence = strings.Join(lines, "+")

	switch {
	case q.HMM != nil:
		q.HMMSupportStart, q.HMMSupportStop = lodSupport(sm, q.HMM.Pos, hmmSupportDrop)
		q.LeadMethod, q.LeadPos = methodHMM, (q.HMMSupportStart+q.HMMSupportStop)/2
	case q.BRM != nil:
		q.LeadMethod, q.LeadPos = methodBRM, q.BRM.Pos
	default:
		q.LeadMethod, q.LeadPos = methodDsiG, q.DsiG.Pos
	}

	q.Direction = "none"
	if signedDSI && q.BRM != nil {
		pos, neg := false, false
		for _, r := range byMethod[methodBRM] {
			pos, neg = pos || r.peak.Value > 0, neg || r.peak.Value < 0
		}
		switch {
		case pos && neg:
			q.Direction = "mixed"
		case pos:
			q.Direction = "+"
		case neg:
			q.Direction = "−"
		}
	}

	q.Start, q.Stop = locusCore(byMethod, support, q.LeadMethod, float64(q.LeadPos))

	seen := map[string]bool{}
	var stats []string
	for _, p := range individual {
		if p.End >= float64(q.SpanStart) && p.Start <= float64(q.SpanStop) && !seen[p.Stat] {
			seen[p.Stat] = true
			stats = append(stats, p.Stat)
		}
	}
	slices.Sort(stats)
	q.Statistics = strings.Join(stats, ",")
	return q
}

// hmmSupportDrop is the LOD drop defining the HMM LOD support interval.
const hmmSupportDrop = 1.5

// lodSupport returns the HMM LOD support interval around the marker at pos
// (one chromosome's markers, sorted by position): the contiguous markers
// whose LOD is within drop of the LOD at pos.
func lodSupport(sm []SmoothedStats, pos int64, drop float64) (int64, int64) {
	i, ok := slices.BinarySearchFunc(sm, pos, func(s SmoothedStats, p int64) int { return cmp.Compare(s.POS, p) })
	if !ok {
		return pos, pos
	}
	floor := sm[i].HMMLOD - drop
	lo, hi := i, i
	for lo > 0 && sm[lo-1].HMMLOD >= floor {
		lo--
	}
	for hi < len(sm)-1 && sm[hi+1].HMMLOD >= floor {
		hi++
	}
	return sm[lo].POS, sm[hi].POS
}

// locusCore is the intersection of the supporting methods' coverage (each
// method's regions in the locus): the piece containing the lead peak, else
// the longest. Without overlap, the lead method's region nearest the lead
// peak.
func locusCore(byMethod map[string][]consensusRegion, support []string, lead string, leadPos float64) (int64, int64) {
	type iv struct{ s, e float64 }
	cover := func(m string) []iv {
		var out []iv
		for _, r := range byMethod[m] {
			out = append(out, iv{r.start, r.end})
		}
		slices.SortFunc(out, func(a, b iv) int { return cmp.Compare(a.s, b.s) })
		return out
	}
	core := cover(support[0])
	for _, m := range support[1:] {
		var next []iv
		for _, a := range core {
			for _, b := range cover(m) {
				if s, e := math.Max(a.s, b.s), math.Min(a.e, b.e); s <= e {
					next = append(next, iv{s, e})
				}
			}
		}
		core = next
	}
	pick := func(ivs []iv) (int64, int64) {
		best := ivs[0]
		for _, v := range ivs {
			if v.s <= leadPos && leadPos <= v.e {
				best = v
				break
			}
			if v.e-v.s > best.e-best.s {
				best = v
			}
		}
		return int64(math.Round(best.s)), int64(math.Round(best.e))
	}
	if len(core) > 0 {
		return pick(core)
	}
	return pick(cover(lead))
}
