package stats

import (
	"fmt"
	"math"
	"sort"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// RegionFilter drops called regions that are too narrow or contain too few
// markers to be a credible QTL (--min-qtl-width, --min-qtl-markers).
//
// Smoothing spreads every real signal over roughly ±σ, so a region above the
// threshold that is much narrower than σ, or supported by only a handful of
// markers, is more likely a local artefact than a QTL. The genome-wide
// thresholds are calibrated without this filter, so it only ever makes calls
// more conservative.
type RegionFilter struct {
	MinWidth   int64 // bp; region End − Start must be at least this
	MinMarkers int   // markers inside [Start, End]

	pos map[string][]int64 // sorted marker positions per chromosome
}

// NewRegionFilter builds the filter from the config. MinQTLWidth = 0 means
// the smoothing σ (WindowSize / 2); negative values disable a criterion.
func NewRegionFilter(cfg utils.AnalysisConfig, smoothed []SmoothedStats) *RegionFilter {
	f := &RegionFilter{MinWidth: cfg.MinQTLWidth, MinMarkers: cfg.MinQTLMarkers, pos: map[string][]int64{}}
	if f.MinWidth == 0 {
		f.MinWidth = int64(cfg.WindowSize / 2)
	}
	for _, s := range smoothed {
		f.pos[s.CHROM] = append(f.pos[s.CHROM], s.POS)
	}
	for c := range f.pos {
		p := f.pos[c]
		if !sort.SliceIsSorted(p, func(i, j int) bool { return p[i] < p[j] }) {
			sort.Slice(p, func(i, j int) bool { return p[i] < p[j] })
		}
	}
	return f
}

// markers returns the number of markers on chrom within [start, end].
func (f *RegionFilter) markers(chrom string, start, end float64) int {
	p := f.pos[chrom]
	lo := sort.Search(len(p), func(i int) bool { return float64(p[i]) >= math.Floor(start) })
	hi := sort.Search(len(p), func(i int) bool { return float64(p[i]) > math.Ceil(end) })
	return hi - lo
}

// Keep reports whether the region passes both criteria. A nil filter keeps
// everything.
func (f *RegionFilter) Keep(chrom string, start, end float64) bool {
	if f == nil {
		return true
	}
	if f.MinWidth > 0 && end-start < float64(f.MinWidth) {
		return false
	}
	if f.MinMarkers > 0 && f.markers(chrom, start, end) < f.MinMarkers {
		return false
	}
	return true
}

// report prints how many regions a method lost to the filter.
func (f *RegionFilter) report(method string, dropped int) {
	if f == nil || dropped == 0 {
		return
	}
	color.Yellow("%s: dropped %d region(s) narrower than %s or with fewer than %d markers (--min-qtl-width, --min-qtl-markers)",
		method, dropped, fmtBp(f.MinWidth), f.MinMarkers)
}

func fmtBp(bp int64) string {
	if bp >= 1_000_000 {
		return fmt.Sprintf("%.2f Mb", float64(bp)/1e6)
	}
	return fmt.Sprintf("%d bp", bp)
}
