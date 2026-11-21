package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
)

// BSAstats holds one site's raw statistics. Genotypes and AD strings are not
// stored: only the raw TSV needs them, and it reads them from the filter.Site.
type BSAstats struct {
	CHROM string
	POS   int64
	REF   string
	ALT   string
	// HighIsAlt: the high-parent allele is ALT (from the filter's polarity).
	HighIsAlt bool

	// Reads on other ALTs of a split multi-allelic record (not in H/L).
	HighBulkOther int
	LowBulkOther  int

	HighBulkL int
	HighBulkH int
	LowBulkL  int
	LowBulkH  int
	HighSI    float64
	LowSI     float64

	DeltaSI      float64
	Gstat        float64
	ED           float64
	OneBulkP0    float64
	OneBulkAFDev float64
	OneBulkGstat float64

	// HMMLOD is the linkage-aware LOD from the pool HMM (see hmm.go),
	// computed along each chromosome after the per-site statistics (RawStats).
	// It is a separate track: not smoothed, normalised or consolidated.
	HMMLOD float64

	// HighDepth and LowDepth are each bulk's informative depth (REF + target
	// ALT reads; 0 when the mode lacks that bulk).
	HighDepth int
	LowDepth  int
	// Depth is the depth used for smoothing weights and the CompositeZ null:
	// the bulk's depth in one-bulk modes, and in two-bulk modes the effective
	// depth 2·dH·dL/(dH+dL), i.e. the equal depth that gives ΔSI the same
	// read-sampling variance, p(1−p)(1/dH + 1/dL).
	Depth int
}

func RawStats(cfg utils.AnalysisConfig, bsaType string, idxs []int, passedVariants []filter.Site) ([]BSAstats, error) {
	if len(idxs) == 0 {
		return nil, fmt.Errorf("no sample indices supplied for %s raw stats", bsaType)
	}

	highParentIdx, lowParentIdx, highBulkIdx, lowBulkIdx := -1, -1, -1, -1
	switch bsaType {
	case "2p2b":
		highParentIdx, lowParentIdx, highBulkIdx, lowBulkIdx = 0, 1, 2, 3
	case "2phb":
		highParentIdx, lowParentIdx, highBulkIdx = 0, 1, 2
	case "2plb":
		highParentIdx, lowParentIdx, lowBulkIdx = 0, 1, 2
	case "hp2b":
		highParentIdx, highBulkIdx, lowBulkIdx = 0, 1, 2
	case "lp2b":
		lowParentIdx, highBulkIdx, lowBulkIdx = 0, 1, 2
	case "hphb":
		highParentIdx, highBulkIdx = 0, 1
	case "hplb":
		highParentIdx, lowBulkIdx = 0, 1
	case "lphb":
		lowParentIdx, highBulkIdx = 0, 1
	case "lplb":
		lowParentIdx, lowBulkIdx = 0, 1
	case "2b":
		highBulkIdx, lowBulkIdx = 0, 1
	default:
		return nil, fmt.Errorf("unsupported bsaseq type %q", bsaType)
	}

	// The filter emits sites with exactly the kept samples, so a mismatch is a
	// bug upstream; fail loudly instead of silently dropping sites.
	for i := range passedVariants {
		if n := int(passedVariants[i].NSamples); n != len(idxs) {
			s := &passedVariants[i]
			return nil, fmt.Errorf("%s raw stats: site %s:%d has %d samples, expected %d", bsaType, s.Chrom, s.Pos, n, len(idxs))
		}
	}

	layout := rawLayout{
		highParent: highParentIdx >= 0, lowParent: lowParentIdx >= 0,
		highBulk: highBulkIdx >= 0, lowBulk: lowBulkIdx >= 0,
		// Without known polarity (2b, --ref-parent none) the sign of ΔSI is
		// arbitrary at each SNP and would cancel when smoothed, so DeltaSI holds |ΔSI|.
		signedDSI: cfg.SignedDeltaSI(bsaType),
		hmm:       HMMAvailable(cfg, bsaType),
	}
	hasOneBulk := layout.highBulk != layout.lowBulk
	signedDSI := layout.signedDSI
	// One-bulk null frequency of the high-parent allele, by polarity:
	// p0[0] where the high allele is REF, p0[1] where it is ALT. See
	// oneBulkP0ByPolarity.
	var oneBulkP0 [2]float64
	if hasOneBulk {
		p0, err := ExpectedAF(cfg.Population)
		if err != nil {
			return nil, err
		}
		bulk := highBulkIdx
		if bulk < 0 {
			bulk = lowBulkIdx
		}
		oneBulkP0 = oneBulkP0ByPolarity(cfg.OneBulkP0Mode, passedVariants, bulk, p0)
	}

	outDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(outDir, 0775); err != nil {
		return nil, err
	}

	color.Cyan("============================ Calculating Raw Statistics (%s) ============================\n\n", bsaType)

	calc := rawCalc{
		roles:     siteRoles{hp: highParentIdx, lp: lowParentIdx, hb: highBulkIdx, lb: lowBulkIdx},
		signedDSI: signedDSI,
		oneBulk:   hasOneBulk,
		p0:        oneBulkP0,
	}
	results := computeRaw(passedVariants, calc)

	// HMM LOD: a linkage- and pool-aware statistic computed along each
	// chromosome from the per-site read counts (hmm.go). It stays a separate
	// track in LOD units: carried through smoothing unsmoothed, not
	// normalised, not in CompositeZ, with its own genome-wide thresholds.
	if layout.hmm {
		setup, err := newHMMSetup(cfg, bsaType)
		if err != nil {
			return nil, err
		}
		if cfg.CMPerMb > 0 {
			color.Cyan("HMM LOD: recombination %.3g cM/Mb\n", cfg.CMPerMb)
		} else {
			color.Cyan("HMM LOD: recombination assumed 1 Morgan per chromosome (set --cm-per-mb to change)\n")
		}
		computeHMMLOD(setup, results, newMarkerOrder(results))
	} else {
		color.Yellow("HMM LOD not computed: 2b without --ref-parent has no allele polarity\n")
	}

	rawPath := filepath.Join(outDir, fmt.Sprintf("GoBSAseq.%s.raw.tsv", bsaType))
	if err := writeRaw(rawPath, results, passedVariants, calc.roles, layout); err != nil {
		return nil, fmt.Errorf("writing %s: %w", rawPath, err)
	}

	color.Green("Raw stats written to %s (%d variants)", rawPath, len(results))
	return results, nil
}

// siteRoles are the positions of each role in filter.Site.Samples (-1 = absent).
type siteRoles struct{ hp, lp, hb, lb int }

// rawCalc computes one site's raw statistics.
type rawCalc struct {
	roles     siteRoles
	signedDSI bool       // false: DeltaSI holds |ΔSI| (2b without --ref-parent)
	oneBulk   bool       // exactly one bulk: compute the one-bulk statistics
	p0        [2]float64 // one-bulk null frequency: [high allele REF, high allele ALT]
}

// site returns the raw statistics of one site. Values are stored at full
// precision; rounding happens only in the TSV. (Rounding stored values to
// 6 dp zeroed ED⁴ whenever |ΔSI| < 0.027, about 16% of null sites.)
func (c *rawCalc) site(site *filter.Site) BSAstats {
	s := BSAstats{CHROM: site.Chrom, POS: site.Pos, REF: site.Ref, ALT: site.Alt, HighIsAlt: site.HighIsAlt}
	// The filter has already decided polarity (which allele the high parent
	// carries; in 2b mode from --ref-parent).
	highIsAlt := site.HighIsAlt

	// bulk returns a bulk's high-allele and low-allele reads, other-ALT reads,
	// high-allele fraction and informative depth.
	bulk := func(b *filter.SampleAD) (h, l, other int, si float64, depth int) {
		h, l = int(b.Ref), int(b.Alt)
		if highIsAlt {
			h, l = l, h
		}
		if depth = h + l; depth > 0 {
			si = float64(h) / float64(depth)
		}
		return h, l, int(b.Other), si, depth
	}
	if c.roles.hb >= 0 {
		s.HighBulkH, s.HighBulkL, s.HighBulkOther, s.HighSI, s.HighDepth = bulk(&site.Samples[c.roles.hb])
	}
	if c.roles.lb >= 0 {
		s.LowBulkH, s.LowBulkL, s.LowBulkOther, s.LowSI, s.LowDepth = bulk(&site.Samples[c.roles.lb])
	}
	s.Depth = effectiveDepth(s.HighDepth, s.LowDepth)

	if c.roles.hb >= 0 && c.roles.lb >= 0 {
		s.DeltaSI = s.HighSI - s.LowSI
		if !c.signedDSI {
			s.DeltaSI = math.Abs(s.DeltaSI)
		}
		s.Gstat = GStatistic(s.HighBulkH, s.HighBulkL, s.LowBulkH, s.LowBulkL)
		s.ED = euclideanDistance4(s.HighSI, s.LowSI)
	}

	if c.oneBulk {
		p0 := c.p0[0]
		if highIsAlt {
			p0 = c.p0[1]
		}
		h, l, si := s.HighBulkH, s.HighBulkL, s.HighSI
		if c.roles.lb >= 0 {
			h, l, si = s.LowBulkH, s.LowBulkL, s.LowSI
		}
		s.OneBulkP0 = p0
		s.OneBulkAFDev = si - p0
		s.OneBulkGstat = oneBulkGStatistic(h, l, p0)
	}
	return s
}

// rawChunk is the number of sites a worker computes or formats at a time.
const rawChunk = 8192

// computeRaw computes every site's statistics, in parallel chunks.
func computeRaw(sites []filter.Site, calc rawCalc) []BSAstats {
	results := make([]BSAstats, len(sites))
	nChunks := (len(sites) + rawChunk - 1) / rawChunk
	bar := progressbar.Default(int64(len(sites)), "Processing variants")
	var next atomic.Int64
	var wg sync.WaitGroup
	for range max(runtime.GOMAXPROCS(0), 1) {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for {
				c := int(next.Add(1)) - 1
				if c >= nChunks {
					return
				}
				lo, hi := c*rawChunk, min((c+1)*rawChunk, len(sites))
				for i := lo; i < hi; i++ {
					results[i] = calc.site(&sites[i])
				}
				_ = bar.Add(hi - lo)
			}
		}()
	}
	wg.Wait()
	_ = bar.Finish()
	return results
}

// writeRaw writes the raw TSV. Workers take chunks of sites in order and
// format their rows into a buffer; this goroutine writes the buffers in chunk
// order. At most 2 chunks per worker are in flight, which bounds the buffered
// output. Every write, flush and close error is returned.
func writeRaw(path string, results []BSAstats, sites []filter.Site, roles siteRoles, lay rawLayout) (err error) {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	w := bufio.NewWriterSize(f, 1<<20)
	if _, err := w.Write(rawHeader(lay)); err != nil {
		return err
	}

	nChunks := (len(sites) + rawChunk - 1) / rawChunk
	workers := max(runtime.GOMAXPROCS(0), 1)
	done := make([]chan []byte, nChunks)
	for c := range done {
		done[c] = make(chan []byte, 1)
	}
	tokens := make(chan struct{}, 2*workers)
	bufs := sync.Pool{New: func() any { return make([]byte, 0, 1<<20) }}

	// A chunk index is claimed only after taking a token, so the chunk this
	// goroutine waits for next is always claimed: no deadlock.
	var next atomic.Int64
	var wg sync.WaitGroup
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for {
				tokens <- struct{}{}
				c := int(next.Add(1)) - 1
				if c >= nChunks {
					<-tokens
					return
				}
				lo, hi := c*rawChunk, min((c+1)*rawChunk, len(sites))
				buf := bufs.Get().([]byte)[:0]
				for i := lo; i < hi; i++ {
					buf = appendRawRow(buf, &results[i], &sites[i], roles, lay)
				}
				done[c] <- buf
			}
		}()
	}

	var werr error
	for c := range done {
		buf := <-done[c]
		if werr == nil {
			_, werr = w.Write(buf)
		}
		bufs.Put(buf) //nolint:staticcheck // []byte header copy is fine here
		<-tokens
	}
	wg.Wait()
	if werr != nil {
		return werr
	}
	return w.Flush()
}

// minOneBulkP0Sites is the fewest informative sites each polarity class needs
// before the empirical mapping-bias estimate is trusted.
const minOneBulkP0Sites = 200

// oneBulkP0ByPolarity returns the one-bulk null frequency of the high-parent
// allele for sites where that allele is REF ([0]) and where it is ALT ([1]).
//
// Reference mapping bias adds a shift b to the REF read fraction. The
// high-allele fraction is therefore theory + b where the high allele is REF
// and theory − b where it is ALT, and a single theoretical p0 turns b into a
// genome-wide AFDev/G1/LOD1/BF1 offset. With mode "empirical" (the default)
// b is estimated as half the difference between the two classes' typical
// high-allele fractions, and p0 = theory ± b. A QTL, or a wrong population
// setting, moves both classes the same way, so it cancels in b and does not
// leak into p0. Mode "theoretical" uses theory for both.
//
// Each class's typical fraction is the interquartile mean (mean of the middle
// 50%), which is robust to QTL regions and, unlike a median, does not snap to
// the 1/depth grid of read fractions. Only sites where the bulk shows at least
// 2 reads of each allele are used, so non-segregating sites (bulk fixed for
// one allele, common in one-parent one-bulk modes) are left out. If either
// class has fewer than minOneBulkP0Sites such sites, b cannot be separated
// from other shifts and the theoretical value is used.
func oneBulkP0ByPolarity(mode string, sites []filter.Site, bulk int, theory float64) [2]float64 {
	p0 := [2]float64{theory, theory}
	if mode == "theoretical" {
		color.Blue("One-bulk p₀: theoretical %.4f (--one-bulk-p0 theoretical)", theory)
		return p0
	}
	var si [2][]float64
	for i := range sites {
		b := sites[i].Samples[bulk]
		if b.Ref < 2 || b.Alt < 2 {
			continue
		}
		k, h := 0, b.Ref
		if sites[i].HighIsAlt {
			k, h = 1, b.Alt
		}
		si[k] = append(si[k], float64(h)/float64(b.Ref+b.Alt))
	}
	if len(si[0]) < minOneBulkP0Sites || len(si[1]) < minOneBulkP0Sites {
		color.Yellow("One-bulk p₀: too few informative sites to estimate mapping bias (%d with the high allele on REF, %d on ALT; need %d each); using theoretical %.4f",
			len(si[0]), len(si[1]), minOneBulkP0Sites, theory)
		return p0
	}
	sort.Float64s(si[0])
	sort.Float64s(si[1])
	refClass, altClass := interquartileMean(si[0]), interquartileMean(si[1])
	bias := (refClass - altClass) / 2
	p0 = [2]float64{clampP0(theory + bias), clampP0(theory - bias)}
	color.Blue("One-bulk p₀: REF mapping bias %+.4f (from %d + %d sites), so p₀ = %.4f where the high allele is REF and %.4f where it is ALT (theoretical %.4f)",
		bias, len(si[0]), len(si[1]), p0[0], p0[1], theory)
	if math.Abs(bias) > 0.1 {
		color.Yellow("One-bulk p₀: estimated mapping bias %+.3f is unusually large; check the alignments, or use --one-bulk-p0 theoretical", bias)
	}
	return p0
}

func clampP0(p float64) float64 { return math.Min(math.Max(p, 0.01), 0.99) }

// interquartileMean is the mean of the middle half of sorted v.
func interquartileMean(v []float64) float64 {
	lo, hi := len(v)/4, len(v)-len(v)/4
	sum := 0.0
	for _, x := range v[lo:hi] {
		sum += x
	}
	return sum / float64(hi-lo)
}

// rawLayout records which samples a mode has (from the role indices RawStats
// already resolved) and whether ΔSI is signed; it decides the TSV columns.
type rawLayout struct {
	highParent, lowParent, highBulk, lowBulk bool
	signedDSI                                bool
	hmm                                      bool // HMMLOD column (HMMAvailable)
}

// rawHeader returns the raw TSV header line.
//
// AD columns are the biallelic view (REF, target ALT). Reads on any other ALT
// of a split multi-allelic record are in the *Other columns, so AD + Other
// accounts for every allele read.
func rawHeader(lay rawLayout) []byte {
	header := []string{"CHROM", "POS", "REF", "ALT"}
	if lay.highParent {
		header = append(header, "HighParGT")
	}
	if lay.lowParent {
		header = append(header, "LowParGT")
	}
	if lay.highBulk {
		header = append(header, "HighBulkGT", "HighBulkAD", "HighBulkOther", "HighBulkL", "HighBulkH", "HighSI")
	}
	if lay.lowBulk {
		header = append(header, "LowBulkGT", "LowBulkAD", "LowBulkOther", "LowBulkL", "LowBulkH", "LowSI")
	}
	if lay.highBulk && lay.lowBulk {
		header = append(header, DeltaSIName(lay.signedDSI), "Gstat", "ED4")
		if lay.hmm {
			header = append(header, "HMMLOD")
		}
		header = append(header, "HighDepth", "LowDepth")
	}
	if lay.highBulk != lay.lowBulk {
		header = append(header, "P0", "AFDev", "Gstat1")
		if lay.hmm {
			header = append(header, "HMMLOD")
		}
	}
	header = append(header, "Depth")
	return []byte(strings.Join(header, "\t") + "\n")
}

// appendRawRow appends one site's TSV row to buf, without allocating.
// Floats use 6 decimals, except ED⁴ (6 significant digits: it is routinely
// below 1e-6, which fixed decimals would print as 0).
func appendRawRow(buf []byte, s *BSAstats, site *filter.Site, r siteRoles, lay rawLayout) []byte {
	tab := func(b []byte) []byte { return append(b, '\t') }
	f6 := func(b []byte, v float64) []byte { return strconv.AppendFloat(tab(b), v, 'f', 6, 64) }
	itoa := func(b []byte, v int) []byte { return strconv.AppendInt(tab(b), int64(v), 10) }

	buf = append(buf, s.CHROM...)
	buf = strconv.AppendInt(tab(buf), s.POS, 10)
	buf = append(tab(buf), s.REF...)
	buf = append(tab(buf), s.ALT...)
	if lay.highParent {
		buf = appendGT(tab(buf), &site.Samples[r.hp])
	}
	if lay.lowParent {
		buf = appendGT(tab(buf), &site.Samples[r.lp])
	}
	bulkCols := func(b []byte, sa *filter.SampleAD, other, l, h int, si float64) []byte {
		b = appendGT(tab(b), sa)
		b = strconv.AppendInt(tab(b), int64(sa.Ref), 10)
		b = strconv.AppendInt(append(b, ','), int64(sa.Alt), 10)
		return f6(itoa(itoa(itoa(b, other), l), h), si)
	}
	if lay.highBulk {
		buf = bulkCols(buf, &site.Samples[r.hb], s.HighBulkOther, s.HighBulkL, s.HighBulkH, s.HighSI)
	}
	if lay.lowBulk {
		buf = bulkCols(buf, &site.Samples[r.lb], s.LowBulkOther, s.LowBulkL, s.LowBulkH, s.LowSI)
	}
	if lay.highBulk && lay.lowBulk {
		buf = f6(f6(buf, s.DeltaSI), s.Gstat)
		buf = strconv.AppendFloat(tab(buf), s.ED, 'g', 6, 64)
		if lay.hmm {
			buf = f6(buf, s.HMMLOD)
		}
		buf = itoa(itoa(buf, s.HighDepth), s.LowDepth)
	}
	if lay.highBulk != lay.lowBulk {
		buf = f6(f6(f6(buf, s.OneBulkP0), s.OneBulkAFDev), s.OneBulkGstat)
		if lay.hmm {
			buf = f6(buf, s.HMMLOD)
		}
	}
	buf = itoa(buf, s.Depth)
	return append(buf, '\n')
}

// appendGT writes a sample's biallelic genotype as filter.SampleAD.GT would
// print with %v: allele indices (-1 missing, 0 REF, 1 ALT, 2 other ALT),
// sorted, space-separated, in brackets.
func appendGT(buf []byte, sa *filter.SampleAD) []byte {
	buf = append(buf, '[')
	first := true
	for _, a := range [4]struct {
		tok string
		n   uint8
	}{{"-1", sa.GTMissing}, {"0", sa.GTRef}, {"1", sa.GTAlt}, {"2", sa.GTOther}} {
		for range a.n {
			if !first {
				buf = append(buf, ' ')
			}
			buf = append(buf, a.tok...)
			first = false
		}
	}
	return append(buf, ']')
}

// effectiveDepth combines the two bulk depths into the single equal depth with
// the same ΔSI read-sampling variance (harmonic mean). With one bulk (the
// other depth 0) it is that bulk's depth.
func effectiveDepth(high, low int) int {
	switch {
	case high <= 0:
		return max(low, 0)
	case low <= 0:
		return high
	}
	return int(math.Round(2 * float64(high) * float64(low) / float64(high+low)))
}

// DeltaSIName is the column name for ΔSI: "DeltaSI", or "AbsDeltaSI" when
// polarity is unknown and the column holds |ΔSI|.
func DeltaSIName(signed bool) string {
	if signed {
		return "DeltaSI"
	}
	return "AbsDeltaSI"
}

// GStatistic is the standard 2×2 G-test (Magwene et al. 2011) of allele counts
// between bulks: a1/b1 are the two allele counts in bulk 1, a2/b2 in bulk 2.
// G = 2 Σ o·ln(o/e) over the four cells, with 0·ln 0 = 0 and no pseudo-counts.
// It is symmetric in the allele order, so callers may pass (H, L) or (ALT, REF).
// A bulk or allele with no reads gives G = 0 (no information).
func GStatistic(a1, b1, a2, b2 int) float64 {
	o := [4]float64{float64(a1), float64(b1), float64(a2), float64(b2)}
	n1, n2 := o[0]+o[1], o[2]+o[3]
	ca, cb := o[0]+o[2], o[1]+o[3]
	total := n1 + n2
	if n1 == 0 || n2 == 0 || ca == 0 || cb == 0 {
		return 0
	}
	e := [4]float64{n1 * ca / total, n1 * cb / total, n2 * ca / total, n2 * cb / total}
	g := 0.0
	for k := range o {
		if o[k] > 0 {
			g += o[k] * math.Log(o[k]/e[k])
		}
	}
	return 2 * g
}

// euclideanDistance4 is ED⁴ (Hill et al. 2013) for a biallelic site: the
// Euclidean distance between the bulks' allele-frequency vectors is
// √((ΔSI)² + (−ΔSI)²) = √2·|ΔSI|, so ED⁴ = 4·ΔSI⁴.
func euclideanDistance4(hSI, lSI float64) float64 {
	d := hSI - lSI
	return 4 * d * d * d * d
}

// ExpectedAF is the population's expected high-parent allele frequency
// under H₀ (utils.ParsePopulation).
func ExpectedAF(population string) (float64, error) {
	return utils.ExpectedAF(population)
}

func oneBulkGStatistic(success, fail int, p0 float64) float64 {
	n := float64(success + fail)
	if n == 0 || p0 <= 0 || p0 >= 1 {
		return 0
	}

	obsSuccess := float64(success)
	obsFail := float64(fail)
	expSuccess := n * p0
	expFail := n * (1 - p0)

	g := 0.0
	if obsSuccess > 0 {
		g += obsSuccess * math.Log(obsSuccess/expSuccess)
	}
	if obsFail > 0 {
		g += obsFail * math.Log(obsFail/expFail)
	}
	return 2 * g
}
