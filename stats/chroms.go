package stats

import (
	"cmp"
	"slices"
)

// chromSlice is one chromosome's smoothed markers, sorted by position.
type chromSlice struct {
	chrom string
	sm    []SmoothedStats
}

// byChromosome splits sm by chromosome, chromosomes sorted by name and each
// chromosome's markers by position (stable, so duplicate positions keep their
// order). A chromosome whose markers are one contiguous, position-sorted run
// of sm (the usual case: smoothed output keeps VCF order) is a subslice of sm,
// not a copy; only other chromosomes are copied and sorted. Callers must not
// modify the returned markers.
func byChromosome(sm []SmoothedStats) []chromSlice {
	type span struct {
		start, end int // first run [start, end)
		runs       int
		sorted     bool
	}
	spans := map[string]*span{}
	var names []string
	for i := 0; i < len(sm); {
		j, sorted := i+1, true
		for ; j < len(sm) && sm[j].CHROM == sm[i].CHROM; j++ {
			sorted = sorted && sm[j].POS >= sm[j-1].POS
		}
		s := spans[sm[i].CHROM]
		if s == nil {
			s = &span{start: i, end: j, sorted: sorted}
			spans[sm[i].CHROM] = s
			names = append(names, sm[i].CHROM)
		}
		s.runs++
		i = j
	}
	slices.Sort(names)
	out := make([]chromSlice, 0, len(names))
	for _, name := range names {
		s := spans[name]
		if s.runs == 1 && s.sorted {
			out = append(out, chromSlice{name, sm[s.start:s.end:s.end]})
			continue
		}
		var c []SmoothedStats
		for i := s.start; i < len(sm); i++ {
			if sm[i].CHROM == name {
				c = append(c, sm[i])
			}
		}
		slices.SortStableFunc(c, func(a, b SmoothedStats) int { return cmp.Compare(a.POS, b.POS) })
		out = append(out, chromSlice{name, c})
	}
	return out
}
