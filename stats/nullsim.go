package stats

// Genome-wide null simulation for significance thresholds.
//
// Smoothed BSA statistics are dominated, under the null, by genetic drift: the
// handful of individuals in each bulk carry long parental segments, so linked
// markers deviate together. A null built from independent single SNPs cannot
// reproduce that, so thresholds are instead derived by simulating whole null
// genomes and running them through the same smoothing, normalisation and
// CompositeZ formula as the real data:
//
//  1. For every bulk individual, simulate the cross (F2/F3/F4/RIL selfing or
//     BCn backcrossing) chromosome by chromosome with Haldane crossovers.
//  2. At every real marker, the pool allele frequency is the simulated dosage;
//     reads are Binomial(observed depth at that marker, pool frequency).
//  3. Compute the raw statistics with the pipeline's own functions, smooth with
//     the same depth-weighted Gaussian kernel, robust-Z normalise genome-wide,
//     and combine with compositeZ().
//  4. Record the genome-wide maximum (and minimum, for signed tracks).
//
// Thresholds are quantiles of those per-genome extremes, so "P99" controls the
// genome-wide false-positive rate at 1%: across null genomes, only 1% produce
// any call anywhere.

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/stat/distuv"
)

// morgansPerChromosome is the assumed genetic length of every chromosome:
// 1 Morgan (100 cM), i.e. on average one crossover per chromosome per meiosis,
// a typical value for plant chromosomes. Thresholds depend on it: shorter
// genetic lengths mean longer co-inherited segments, more drift-driven
// signal, and higher thresholds.
const morgansPerChromosome = 1.0

// rilSelfingGenerations approximates a RIL as F8 (>99% homozygous).
const rilSelfingGenerations = 7

// ── Cross design ───────────────────────────────────────────────────────────

// crossDesign describes how each bulk individual is produced from the F1.
type crossDesign struct {
	backcrosses int // number of backcrosses to the recurrent parent
	recurrentH  bool
	selfings    int // number of selfing generations after backcrossing
}

func parseCrossDesign(pop string) (crossDesign, error) {
	p := strings.ToUpper(strings.TrimSpace(pop))
	switch p {
	case "F2":
		return crossDesign{selfings: 1}, nil
	case "F3":
		return crossDesign{selfings: 2}, nil
	case "F4":
		return crossDesign{selfings: 3}, nil
	case "RIL":
		return crossDesign{selfings: rilSelfingGenerations}, nil
	}
	if strings.HasPrefix(p, "BC") {
		suffix := p[2:]
		d := crossDesign{}
		switch {
		case strings.HasSuffix(suffix, "F2H"), strings.HasSuffix(suffix, "F2L"):
			d.recurrentH = strings.HasSuffix(suffix, "H")
			d.selfings = 1
			suffix = suffix[:len(suffix)-3]
		case strings.HasSuffix(suffix, "H"), strings.HasSuffix(suffix, "L"):
			d.recurrentH = strings.HasSuffix(suffix, "H")
			suffix = suffix[:len(suffix)-1]
		default:
			return d, fmt.Errorf("unknown population structure: %s", pop)
		}
		n, err := strconv.Atoi(suffix)
		if err != nil || n < 1 {
			return d, fmt.Errorf("unknown population structure: %s", pop)
		}
		d.backcrosses = n
		return d, nil
	}
	return crossDesign{}, fmt.Errorf("unknown population structure: %s", pop)
}

// ── Haplotypes ─────────────────────────────────────────────────────────────

// haplotype is a parental mosaic along one chromosome: the allele (1 = high
// parent, 0 = low parent) at position 0, flipping at each switch point.
type haplotype struct {
	start    uint8
	switches []float64 // sorted, in bp
}

func (h haplotype) alleleAt(pos float64) uint8 {
	// Number of switches at or before pos.
	n := sort.Search(len(h.switches), func(i int) bool { return h.switches[i] > pos })
	return h.start ^ uint8(n&1)
}

// gamete returns a recombinant of a and b with Haldane (Poisson) crossovers.
func gamete(rng *rand.Rand, a, b haplotype, lengthBp, morgans float64) haplotype {
	src := [2]haplotype{a, b}
	cur := rng.Intn(2)

	nCO := int(distuv.Poisson{Lambda: morgans, Src: rng}.Rand())
	cos := make([]float64, nCO)
	for i := range cos {
		cos[i] = rng.Float64() * lengthBp
	}
	sort.Float64s(cos)

	// Candidate change points: every crossover and every switch in either parent.
	cands := make([]float64, 0, nCO+len(a.switches)+len(b.switches))
	cands = append(cands, cos...)
	cands = append(cands, a.switches...)
	cands = append(cands, b.switches...)
	sort.Float64s(cands)

	out := haplotype{start: src[cur].alleleAt(0)}
	allele := out.start
	ci := 0
	for _, x := range cands {
		for ci < len(cos) && cos[ci] <= x {
			cur ^= 1
			ci++
		}
		if a2 := src[cur].alleleAt(x); a2 != allele {
			out.switches = append(out.switches, x)
			allele = a2
		}
	}
	return out
}

// simulateIndividual returns the two haplotypes of one bulk individual.
func simulateIndividual(rng *rand.Rand, d crossDesign, lengthBp, morgans float64) [2]haplotype {
	H := haplotype{start: 1}
	L := haplotype{start: 0}
	ind := [2]haplotype{H, L} // F1
	for range d.backcrosses {
		rec := L
		if d.recurrentH {
			rec = H
		}
		ind = [2]haplotype{gamete(rng, ind[0], ind[1], lengthBp, morgans), rec}
	}
	for range d.selfings {
		ind = [2]haplotype{
			gamete(rng, ind[0], ind[1], lengthBp, morgans),
			gamete(rng, ind[0], ind[1], lengthBp, morgans),
		}
	}
	return ind
}

// poolHighAlleleCount simulates a bulk of size n and returns, at each marker
// position (sorted), the number of high-parent alleles among the pool's 2n
// (pool frequency = count / 2n).
func poolHighAlleleCount(rng *rand.Rand, d crossDesign, n int, pos []float64, lengthBp, morgans float64, diff []int, out []int) {
	for i := range diff {
		diff[i] = 0
	}
	addRange := func(a, b float64) {
		lo := sort.SearchFloat64s(pos, a)
		hi := sort.SearchFloat64s(pos, b)
		if lo < hi {
			diff[lo]++
			diff[hi]--
		}
	}
	for range n {
		ind := simulateIndividual(rng, d, lengthBp, morgans)
		for _, h := range ind {
			allele := h.start
			prev := math.Inf(-1)
			for _, sw := range h.switches {
				if allele == 1 {
					addRange(prev, sw)
				}
				allele ^= 1
				prev = sw
			}
			if allele == 1 {
				addRange(prev, math.Inf(1))
			}
		}
	}
	c := 0
	for i := range pos {
		c += diff[i]
		out[i] = c
	}
}

// ── Marker layout shared by every replicate ─────────────────────────────────

type nullChrom struct {
	name     string
	lengthBp float64
	pos      []float64 // raw marker positions (sorted)
	hDepth   []int
	lDepth   []int
	weight   []float64 // depth weight used by the smoother: max(Depth, 1)
	bin      []int     // smoothing bin of each raw marker
	nBins    int
	// Bins at which smoothed values are needed: each bin holding at least one
	// evaluated (real smoothed) marker, once, with the number of markers in
	// it. Markers in the same bin get bit-identical smoothed values, so each
	// bin is computed once and weighted by its marker count when normalising.
	evalBins  []int
	evalCount []int
}

type nullLayout struct {
	chroms     []nullChrom
	binWidth   float64
	sigma      float64
	kernel     []float64 // kernel weight by bin offset, index = offset + radius
	radius     int
	nEval      int   // number of evaluated markers (Σ evalCount)
	counts     []int // evalCount of every evaluated bin, genome-wide, in track order
	hasHigh    bool
	hasLow     bool
	hasBoth    bool
	hasOne     bool
	randomPol  bool // no parent: allele polarity is arbitrary per marker
	p0         float64
	design     crossDesign
	highSize   int
	lowSize    int
	oneSize    int

	hReads, lReads *readSampler // read samplers for the high / low bulk
	lg             *lgammaTables
	brmScale       float64 // PopulationVarianceScale, for BRM's σ
}

func buildNullLayout(raw []BSAstats, smoothed []SmoothedStats, bsaType string, sigma float64, contigLen map[string]float64,
	design crossDesign, p0 float64, highSize, lowSize, oneSize int) *nullLayout {

	hasHigh, hasLow, hasBoth, hasOne := BulkFlags(bsaType)
	L := &nullLayout{
		sigma: sigma, hasHigh: hasHigh, hasLow: hasLow, hasBoth: hasBoth, hasOne: hasOne,
		randomPol: bsaType == "2b", p0: p0, design: design,
		highSize: highSize, lowSize: lowSize, oneSize: oneSize,
	}
	// Bins are ~1% of σ: positional error is negligible next to the kernel width.
	L.binWidth = math.Max(1000, sigma/100)
	L.radius = int(math.Ceil(kernelCutoffSigmas * sigma / L.binWidth))
	L.kernel = make([]float64, 2*L.radius+1)
	for d := -L.radius; d <= L.radius; d++ {
		x := float64(d) * L.binWidth / sigma
		L.kernel[d+L.radius] = math.Exp(-0.5 * x * x)
	}

	hSize, lSize := highSize, lowSize
	if hasOne {
		hSize, lSize = oneSize, oneSize
	}
	maxH, maxL := 0, 0
	var hDepths, lDepths []int
	for _, s := range raw {
		hd, ld := s.HighBulkH+s.HighBulkL, s.LowBulkH+s.LowBulkL
		hDepths, lDepths = append(hDepths, hd), append(lDepths, ld)
		maxH, maxL = max(maxH, hd), max(maxL, ld)
	}
	if hasHigh {
		L.hReads = newReadSampler(2*hSize, hDepths)
	}
	if hasLow {
		L.lReads = newReadSampler(2*lSize, lDepths)
	}
	L.lg = newLgammaTables(maxH + maxL)

	byChrom := map[string][]BSAstats{}
	for _, s := range raw {
		byChrom[s.CHROM] = append(byChrom[s.CHROM], s)
	}
	evalByChrom := map[string][]int64{}
	for _, s := range smoothed {
		evalByChrom[s.CHROM] = append(evalByChrom[s.CHROM], s.POS)
	}
	names := make([]string, 0, len(byChrom))
	for n := range byChrom {
		names = append(names, n)
	}
	sort.Strings(names)

	for _, name := range names {
		rs := byChrom[name]
		sort.Slice(rs, func(i, j int) bool { return rs[i].POS < rs[j].POS })
		c := nullChrom{name: name}
		maxPos := 0.0
		for _, s := range rs {
			p := float64(s.POS)
			c.pos = append(c.pos, p)
			c.hDepth = append(c.hDepth, s.HighBulkH+s.HighBulkL)
			c.lDepth = append(c.lDepth, s.LowBulkH+s.LowBulkL)
			c.weight = append(c.weight, math.Max(float64(s.Depth), 1))
			c.bin = append(c.bin, int(p/L.binWidth))
			if p > maxPos {
				maxPos = p
			}
		}
		c.lengthBp = math.Max(maxPos+1, contigLen[name])
		c.nBins = int(maxPos/L.binWidth) + 1
		eb := make([]int, 0, len(evalByChrom[name]))
		for _, p := range evalByChrom[name] {
			eb = append(eb, int(float64(p)/L.binWidth))
		}
		sort.Ints(eb)
		for i, b := range eb {
			if i > 0 && b == eb[i-1] {
				c.evalCount[len(c.evalCount)-1]++
				continue
			}
			c.evalBins = append(c.evalBins, b)
			c.evalCount = append(c.evalCount, 1)
		}
		L.nEval += len(eb)
		L.counts = append(L.counts, c.evalCount...)
		L.chroms = append(L.chroms, c)
	}
	return L
}

// ── One null genome ────────────────────────────────────────────────────────

// nullExtremes holds the genome-wide extremes of one simulated null genome.
type nullExtremes struct {
	HighSIMax, HighSIMin   float64
	LowSIMax, LowSIMin     float64
	DeltaSIMax, DeltaSIMin float64
	GstatMax, ED4Max       float64
	LODMax, BBLogBFMax     float64
	AFDevMax, AFDevMin     float64
	OneGMax, OneLODMax     float64
	OneBBMax               float64
	CompositeZMax          float64
	BRMZMax                float64 // max over the genome of |ΔSI|/σ (or |AFDev|/σ), see brm.go
}

// Indices of the smoothed tracks within a replicate's work arrays.
const (
	tHighSI = iota
	tLowSI
	tDeltaSI
	tG
	tED4
	tLOD
	tBB
	tAFDev
	tG1
	tLOD1
	tBB1
	nTracks
)

type nullWorker struct {
	diff  []int
	afH   []int // high-parent allele count per marker, high bulk
	afL   []int
	wsum  []float64
	wv    [nTracks][]float64
	track [nTracks][]float64 // per evaluated bin, genome-wide (see nullChrom.evalBins)
	rz    robustZScratch
}

func (L *nullLayout) activeTracks() []int {
	var t []int
	if L.hasHigh {
		t = append(t, tHighSI)
	}
	if L.hasLow {
		t = append(t, tLowSI)
	}
	if L.hasBoth {
		t = append(t, tDeltaSI, tG, tED4, tLOD, tBB)
	}
	if L.hasOne {
		t = append(t, tAFDev, tG1, tLOD1, tBB1)
	}
	return t
}

func grow[T any](s []T, n int) []T {
	if cap(s) < n {
		return make([]T, n)
	}
	return s[:n]
}

func (L *nullLayout) simulateGenome(rng *rand.Rand, w *nullWorker) nullExtremes {
	tracks := L.activeTracks()
	for _, t := range tracks {
		w.track[t] = grow(w.track[t], 0)
	}

	for ci := range L.chroms {
		c := &L.chroms[ci]
		n := len(c.pos)
		w.diff = grow(w.diff, n+1)
		w.afH = grow(w.afH, n)
		w.afL = grow(w.afL, n)

		if L.hasHigh {
			size := L.highSize
			if L.hasOne {
				size = L.oneSize
			}
			poolHighAlleleCount(rng, L.design, size, c.pos, c.lengthBp, morgansPerChromosome, w.diff, w.afH)
		}
		if L.hasLow {
			size := L.lowSize
			if L.hasOne {
				size = L.oneSize
			}
			poolHighAlleleCount(rng, L.design, size, c.pos, c.lengthBp, morgansPerChromosome, w.diff, w.afL)
		}

		w.wsum = grow(w.wsum, c.nBins)
		for b := range w.wsum {
			w.wsum[b] = 0
		}
		for _, t := range tracks {
			w.wv[t] = grow(w.wv[t], c.nBins)
			for b := range w.wv[t] {
				w.wv[t][b] = 0
			}
		}

		// Raw statistics per marker, binned with depth weights.
		for k := 0; k < n; k++ {
			flip := L.randomPol && rng.Intn(2) == 1
			var hH, hL, lH, lL int
			var hSI, lSI float64
			if L.hasHigh {
				hH, hL = L.hReads.draw(rng, c.hDepth[k], w.afH[k], flip)
				if t := hH + hL; t > 0 {
					hSI = float64(hH) / float64(t)
				}
			}
			if L.hasLow {
				lH, lL = L.lReads.draw(rng, c.lDepth[k], w.afL[k], flip)
				if t := lH + lL; t > 0 {
					lSI = float64(lH) / float64(t)
				}
			}

			var v [nTracks]float64
			v[tHighSI], v[tLowSI] = hSI, lSI
			if L.hasBoth {
				v[tDeltaSI] = hSI - lSI
				v[tG] = GStatistic(hH, hL, lH, lL)
				v[tED4] = euclideanDistance4(hSI, lSI)
				v[tLOD] = lod(hL, hH, lL, lH)
				v[tBB] = L.lg.betaBinomialLogBF(hH, hL, lH, lL)
			}
			if L.hasOne {
				h, l, si := hH, hL, hSI
				if L.hasLow {
					h, l, si = lH, lL, lSI
				}
				v[tAFDev] = si - L.p0
				v[tG1] = oneBulkGStatistic(h, l, L.p0)
				v[tLOD1] = oneBulkLOD(h, l, L.p0)
				v[tBB1] = oneBulkBetaBinomialLogBF(h, l, L.p0)
			}

			b := c.bin[k]
			wk := c.weight[k]
			w.wsum[b] += wk
			for _, t := range tracks {
				w.wv[t][b] += wk * v[t]
			}
		}

		// Gaussian smoothing evaluated at the bins of the real smoothed markers.
		for _, e := range c.evalBins {
			lo, hi := e-L.radius, e+L.radius
			if lo < 0 {
				lo = 0
			}
			if hi >= c.nBins {
				hi = c.nBins - 1
			}
			den := 0.0
			for b := lo; b <= hi; b++ {
				den += L.kernel[b-e+L.radius] * w.wsum[b]
			}
			for _, t := range tracks {
				num := 0.0
				wv := w.wv[t]
				for b := lo; b <= hi; b++ {
					num += L.kernel[b-e+L.radius] * wv[b]
				}
				val := 0.0
				if den > 0 {
					val = num / den
				}
				w.track[t] = append(w.track[t], val)
			}
		}
	}

	ex := nullExtremes{}
	maxMin := func(xs []float64) (float64, float64) {
		mx, mn := math.Inf(-1), math.Inf(1)
		for _, x := range xs {
			mx = math.Max(mx, x)
			mn = math.Min(mn, x)
		}
		return mx, mn
	}
	if L.hasHigh {
		ex.HighSIMax, ex.HighSIMin = maxMin(w.track[tHighSI])
	}
	if L.hasLow {
		ex.LowSIMax, ex.LowSIMin = maxMin(w.track[tLowSI])
	}

	// CompositeZ: identical normalisation (robustZScores) and formula (compositeZ).
	var zS, zG, zL []float64
	if L.hasBoth {
		ex.DeltaSIMax, ex.DeltaSIMin = maxMin(w.track[tDeltaSI])
		ex.GstatMax, _ = maxMin(w.track[tG])
		ex.ED4Max, _ = maxMin(w.track[tED4])
		ex.LODMax, _ = maxMin(w.track[tLOD])
		ex.BBLogBFMax, _ = maxMin(w.track[tBB])
		zS = w.rz.robustZ(w.track[tDeltaSI], L.counts, 0)
		zG = w.rz.robustZ(w.track[tG], L.counts, 1)
		zL = w.rz.robustZ(w.track[tLOD], L.counts, 2)
	} else if L.hasOne {
		ex.AFDevMax, ex.AFDevMin = maxMin(w.track[tAFDev])
		ex.OneGMax, _ = maxMin(w.track[tG1])
		ex.OneLODMax, _ = maxMin(w.track[tLOD1])
		ex.OneBBMax, _ = maxMin(w.track[tBB1])
		zS = w.rz.robustZ(w.track[tAFDev], L.counts, 0)
		zG = w.rz.robustZ(w.track[tG1], L.counts, 1)
		zL = w.rz.robustZ(w.track[tLOD1], L.counts, 2)
	}
	// BRM: largest |ΔSI|/σ (two bulks) or |SI − p0|/σ (one bulk), with σ
	// computed exactly as RunBRM does on the real data.
	if L.brmScale > 0 {
		if L.hasBoth {
			hi, lo, d := w.track[tHighSI], w.track[tLowSI], w.track[tDeltaSI]
			for i := range d {
				z := math.Abs(d[i]) / brmSigmaTwoBulk((hi[i]+lo[i])/2, L.highSize, L.lowSize, L.brmScale)
				ex.BRMZMax = math.Max(ex.BRMZMax, z)
			}
		} else if L.hasOne {
			s := brmSigmaOneBulk(L.p0, L.oneSize, L.brmScale)
			for _, v := range w.track[tAFDev] {
				ex.BRMZMax = math.Max(ex.BRMZMax, math.Abs(v)/s)
			}
		}
	}

	ex.CompositeZMax = math.Inf(-1)
	for i := range zS {
		ex.CompositeZMax = math.Max(ex.CompositeZMax, compositeZ(zS[i], zG[i], zL[i]))
	}
	return ex
}

// robustZScratch holds reusable buffers for robustZ.
type robustZScratch struct {
	idx []int
	dev []float64
	out [3][]float64
}

// robustZ returns robustZScores of the expanded vector in which vals[i]
// appears counts[i] times, but only once per distinct entry: the result is
// bit-identical to robustZScores on the expanded vector (same median, MAD,
// fallback and rounding). slot selects which output buffer to reuse.
func (s *robustZScratch) robustZ(vals []float64, counts []int, slot int) []float64 {
	const consistencyFactor = 1.4826
	n := len(vals)
	out := grow(s.out[slot], n)
	s.out[slot] = out
	if n == 0 {
		return out
	}
	total := 0
	for _, c := range counts {
		total += c
	}

	med := s.weightedMedian(vals, counts, total)
	s.dev = grow(s.dev, n)
	for i, v := range vals {
		s.dev[i] = math.Abs(v - med)
	}
	mad := s.weightedMedian(s.dev, counts, total) * consistencyFactor

	if mad > 0 {
		for i, v := range vals {
			out[i] = math.Round(((v-med)/mad)*1e6) / 1e6
		}
		return out
	}
	// MAD = 0 fallback: classical Z-score, summed in expanded order.
	var mean float64
	for i, v := range vals {
		for range counts[i] {
			mean += v
		}
	}
	mean /= float64(total)
	var variance float64
	for i, v := range vals {
		d := v - mean
		for range counts[i] {
			variance += d * d
		}
	}
	sd := math.Sqrt(variance / float64(total))
	for i, v := range vals {
		if sd == 0 {
			out[i] = 0
		} else {
			out[i] = math.Round(((v-mean)/sd)*1e6) / 1e6
		}
	}
	return out
}

// weightedMedian is medianOf on the expanded vector (vals[i] repeated
// counts[i] times; total = Σ counts).
func (s *robustZScratch) weightedMedian(vals []float64, counts []int, total int) float64 {
	s.idx = grow(s.idx, len(vals))
	for i := range s.idx {
		s.idx[i] = i
	}
	sort.Slice(s.idx, func(a, b int) bool { return vals[s.idx[a]] < vals[s.idx[b]] })
	at := func(rank int) float64 {
		cum := 0
		for _, i := range s.idx {
			cum += counts[i]
			if rank < cum {
				return vals[i]
			}
		}
		return vals[s.idx[len(s.idx)-1]]
	}
	if total%2 == 1 {
		return at(total / 2)
	}
	return (at(total/2-1) + at(total/2)) / 2
}

// ── Read sampling ──────────────────────────────────────────────────────────

// maxTableDepth caps the depths that get a precomputed CDF; deeper markers
// (rare) fall back to gonum's sampler.
const maxTableDepth = 500

// binomCDF is the CDF of Binomial(depth, c/2n) over k = lo … lo+len(cdf)-1;
// tails below 1e-16 are dropped.
type binomCDF struct {
	lo  int
	cdf []float64
}

// readSampler draws Binomial(depth, c/2n) reads. Pool frequencies can only be
// c/2n (c = 0 … 2n alleles) and depths are the observed ones, so every
// distribution needed is tabulated once and sampled with one uniform draw and
// a binary search, instead of gonum's rejection sampler (lgamma + tan per
// draw, about half the simulation time before).
type readSampler struct {
	twoN int
	tab  [][]binomCDF // [depth][c]; nil for depths not observed
}

func newReadSampler(twoN int, depths []int) *readSampler {
	s := &readSampler{twoN: twoN}
	maxD := 0
	for _, d := range depths {
		if d <= maxTableDepth {
			maxD = max(maxD, d)
		}
	}
	s.tab = make([][]binomCDF, maxD+1)
	pmf := make([]float64, maxD+1)
	for _, d := range depths {
		if d <= 0 || d > maxTableDepth || s.tab[d] != nil {
			continue
		}
		row := make([]binomCDF, twoN+1)
		lgD, _ := math.Lgamma(float64(d) + 1)
		for c := 1; c < twoN; c++ { // c = 0 and 2n are deterministic
			p := float64(c) / float64(twoN)
			lp, lq := math.Log(p), math.Log1p(-p)
			lo, hi := -1, 0
			for k := 0; k <= d; k++ {
				a, _ := math.Lgamma(float64(k) + 1)
				b, _ := math.Lgamma(float64(d-k) + 1)
				pmf[k] = math.Exp(lgD - a - b + float64(k)*lp + float64(d-k)*lq)
				if pmf[k] >= 1e-16 {
					if lo < 0 {
						lo = k
					}
					hi = k
				}
			}
			cdf := make([]float64, hi-lo+1)
			sum := 0.0
			for k := lo; k <= hi; k++ {
				sum += pmf[k]
				cdf[k-lo] = sum
			}
			for i := range cdf {
				cdf[i] /= sum
			}
			cdf[len(cdf)-1] = 1
			row[c] = binomCDF{lo: lo, cdf: cdf}
		}
		s.tab[d] = row
	}
	return s
}

// draw returns (high-allele reads, low-allele reads) at a marker of the given
// depth whose pool carries c high-parent alleles out of 2n. With flip, the
// observed "high" allele is actually the low-parent allele, mimicking
// parent-less runs where polarity is arbitrary per marker.
func (s *readSampler) draw(rng *rand.Rand, depth, c int, flip bool) (int, int) {
	if depth <= 0 {
		return 0, 0
	}
	if flip {
		c = s.twoN - c
	}
	var h int
	switch {
	case c <= 0:
		h = 0
	case c >= s.twoN:
		h = depth
	case depth < len(s.tab) && s.tab[depth] != nil:
		t := s.tab[depth][c]
		h = t.lo + sort.SearchFloat64s(t.cdf, rng.Float64())
	default:
		p := float64(c) / float64(s.twoN)
		h = int(distuv.Binomial{N: float64(depth), P: p, Src: rng}.Rand())
	}
	return h, depth - h
}

// lgammaTables evaluate betaBinomialLogBF from tabulated log-gamma values.
// Every argument it needs is k + 0.5 or a whole number, so the tables give
// bit-identical results without six math.Lgamma calls per marker.
type lgammaTables struct {
	half  []float64 // lgamma(0.5 + k)
	whole []float64 // lgamma(1 + n)
}

func newLgammaTables(maxCount int) *lgammaTables {
	t := &lgammaTables{half: make([]float64, maxCount+1), whole: make([]float64, maxCount+1)}
	for k := 0; k <= maxCount; k++ {
		t.half[k], _ = math.Lgamma(0.5 + float64(k))
		t.whole[k], _ = math.Lgamma(1 + float64(k))
	}
	return t
}

// logBeta(0.5+k, 0.5+m), exactly as logBeta computes it.
func (t *lgammaTables) logBetaHalf(k, m int) float64 {
	return t.half[k] + t.half[m] - t.whole[k+m]
}

// betaBinomialLogBF is betaBinomialLogBF (stats.go) with tabulated lgamma.
func (t *lgammaTables) betaBinomialLogBF(highSucc, highFail, lowSucc, lowFail int) float64 {
	if highSucc+lowSucc+highFail+lowFail >= len(t.half) || highSucc+highFail+lowSucc+lowFail < 0 {
		return betaBinomialLogBF(highSucc, highFail, lowSucc, lowFail)
	}
	b0 := t.logBetaHalf(0, 0)
	logAlt := t.logBetaHalf(highSucc, highFail) - b0
	logAlt += t.logBetaHalf(lowSucc, lowFail) - b0
	logNull := t.logBetaHalf(highSucc+lowSucc, highFail+lowFail) - b0
	return (logAlt - logNull) / math.Log(10)
}

// ── Driver ─────────────────────────────────────────────────────────────────

// simulateNullExtremes runs reps null genomes in parallel. Each replicate's RNG
// is seeded from (seed, replicate index), so results do not depend on goroutine
// scheduling.
func simulateNullExtremes(L *nullLayout, reps int, seed int64) []nullExtremes {
	out := make([]nullExtremes, reps)
	bar := progressbar.NewOptions(reps,
		progressbar.OptionSetDescription("Null genome simulation"),
		progressbar.OptionSetWidth(40),
		progressbar.OptionShowCount(),
	)
	var next atomic.Int64
	var wg sync.WaitGroup
	workers := runtime.GOMAXPROCS(0)
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			w := &nullWorker{}
			for {
				r := int(next.Add(1)) - 1
				if r >= reps {
					return
				}
				rng := rand.New(rand.NewSource(deriveSeed(seed, fmt.Sprintf("null-genome:%d", r))))
				out[r] = L.simulateGenome(rng, w)
				_ = bar.Add(1)
			}
		}()
	}
	wg.Wait()
	fmt.Println()
	return out
}
