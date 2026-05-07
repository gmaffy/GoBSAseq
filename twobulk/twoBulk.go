package twobulk

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"sync"

	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/utils"
)

// -----------------------------------------------------------------------------
// Structures
// -----------------------------------------------------------------------------

type BSAstats struct {
	CHROM string
	POS   int64

	HighSI float64
	LowSI  float64

	DeltaSI float64
	Gstat   float64
	BBLogBF float64

	Depth  int
	Weight float64
}

type Smoothed struct {
	CHROM string
	POS   int64

	ZCombined float64
}

// -----------------------------------------------------------------------------
// GATK Hard Filtering
// -----------------------------------------------------------------------------

func getInfoFloat(v *vcfgo.Variant, key string) (float64, bool) {
	val, err := v.Info().Get(key)
	if err != nil || val == nil {
		return 0, false
	}

	switch t := val.(type) {
	case float64:
		return t, true
	case []float64:
		if len(t) > 0 {
			return t[0], true
		}
	}
	return 0, false
}

func passGATK(v *vcfgo.Variant) bool {
	if qd, ok := getInfoFloat(v, "QD"); ok && qd < 2.0 {
		return false
	}
	if fs, ok := getInfoFloat(v, "FS"); ok && fs > 60.0 {
		return false
	}
	if mq, ok := getInfoFloat(v, "MQ"); ok && mq < 40.0 {
		return false
	}
	if mqrs, ok := getInfoFloat(v, "MQRankSum"); ok && mqrs < -12.5 {
		return false
	}
	if rprs, ok := getInfoFloat(v, "ReadPosRankSum"); ok && rprs < -8.0 {
		return false
	}
	return true
}

// -----------------------------------------------------------------------------
// Variant Weight (variance + quality aware)
// -----------------------------------------------------------------------------

func variantWeight(v *vcfgo.Variant, depth int, hSI, lSI float64) float64 {
	if depth <= 0 {
		return 0
	}

	p := (hSI + lSI) / 2.0
	base := 1.0
	if p > 0 && p < 1 {
		base = float64(depth) / (p * (1 - p))
	}

	w := math.Sqrt(base)

	if qd, ok := getInfoFloat(v, "QD"); ok && qd < 5 {
		w *= qd / 5.0
	}
	if fs, ok := getInfoFloat(v, "FS"); ok && fs > 30 {
		w *= 30.0 / fs
	}
	if mq, ok := getInfoFloat(v, "MQ"); ok && mq < 50 {
		w *= mq / 50.0
	}

	if w < 0.05 {
		w = 0.05
	}
	return w
}

// -----------------------------------------------------------------------------
// Statistics
// -----------------------------------------------------------------------------

func GStatistic(a1, r1, a2, r2 int) float64 {
	h1 := float64(a1 + r1)
	h2 := float64(a2 + r2)
	if h1 == 0 || h2 == 0 {
		return 0
	}

	p1 := float64(a1) / h1
	p2 := float64(a2) / h2
	return math.Abs(p1-p2) * (h1 + h2)
}

func betaBinomialLogBF(a1, r1, a2, r2 int) float64 {
	return math.Log(float64(a1+r1+1)) - math.Log(float64(a2+r2+1))
}

// -----------------------------------------------------------------------------
// Core Pipeline
// -----------------------------------------------------------------------------

func RunTwoBulkV3(rdr *vcfgo.Reader, highIdx, lowIdx int) ([]Smoothed, []QTLPeak) {

	variantChan := make(chan *vcfgo.Variant, 1000)
	statChan := make(chan BSAstats, 1000)

	var wg sync.WaitGroup

	// Workers
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()

			for v := range variantChan {

				if !passGATK(v) || len(v.Alt()) != 1 {
					continue
				}

				h := v.Samples[highIdx]
				l := v.Samples[lowIdx]

				if h.DP == 0 || l.DP == 0 {
					continue
				}

				hRef, _ := h.RefDepth()
				hAlt, _ := h.AltDepths()
				lRef, _ := l.RefDepth()
				lAlt, _ := l.AltDepths()

				if len(hAlt) == 0 || len(lAlt) == 0 {
					continue
				}

				hSI := float64(hAlt[0]) / float64(h.DP)
				lSI := float64(lAlt[0]) / float64(l.DP)

				delta := hSI - lSI
				depth := h.DP
				if l.DP < depth {
					depth = l.DP
				}

				wgt := variantWeight(v, depth, hSI, lSI)

				statChan <- BSAstats{
					CHROM: v.Chromosome,
					POS:   int64(v.Pos),

					HighSI: hSI,
					LowSI:  lSI,

					DeltaSI: delta,
					Gstat:   GStatistic(hAlt[0], hRef, lAlt[0], lRef),
					BBLogBF: betaBinomialLogBF(hAlt[0], hRef, lAlt[0], lRef),

					Depth:  depth,
					Weight: wgt,
				}
			}
		}()
	}

	// Feed
	go func() {
		for v := rdr.Read(); v != nil; v = rdr.Read() {
			variantChan <- v
		}
		close(variantChan)
	}()

	go func() {
		wg.Wait()
		close(statChan)
	}()

	// Collect
	byChr := map[string][]BSAstats{}
	for s := range statChan {
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}

	// -------------------------------------------------------------------------
	// Adaptive smoothing (k-NN)
	// -------------------------------------------------------------------------

	const k = 100

	var smoothed []Smoothed

	for chrom, stats := range byChr {

		sort.Slice(stats, func(i, j int) bool {
			return stats[i].POS < stats[j].POS
		})

		// Precompute arrays
		dsi := make([]float64, len(stats))
		g := make([]float64, len(stats))
		bb := make([]float64, len(stats))

		for i := range stats {
			dsi[i] = stats[i].DeltaSI
			g[i] = stats[i].Gstat
			bb[i] = stats[i].BBLogBF
		}

		// Robust normalization (median + MAD)
		median := func(x []float64) float64 {
			tmp := append([]float64{}, x...)
			sort.Float64s(tmp)
			return tmp[len(tmp)/2]
		}

		mad := func(x []float64, med float64) float64 {
			dev := make([]float64, len(x))
			for i := range x {
				dev[i] = math.Abs(x[i] - med)
			}
			sort.Float64s(dev)
			return dev[len(dev)/2]
		}

		dsiMed := median(dsi)
		dsiMAD := mad(dsi, dsiMed)

		gMed := median(g)
		gMAD := mad(g, gMed)

		bbMed := median(bb)
		bbMAD := mad(bb, bbMed)

		for i := range stats {

			start := i - k/2
			end := i + k/2
			if start < 0 {
				start = 0
			}
			if end >= len(stats) {
				end = len(stats) - 1
			}

			var sum, wsum float64

			for j := start; j <= end; j++ {

				d := math.Abs(float64(stats[j].POS - stats[i].POS))
				kernel := 1.0 / (1.0 + d) // cheap fast kernel

				// Z normalize inline
				z1 := (stats[j].DeltaSI - dsiMed) / (dsiMAD*1.4826 + 1e-6)
				z2 := (stats[j].Gstat - gMed) / (gMAD*1.4826 + 1e-6)
				z3 := (stats[j].BBLogBF - bbMed) / (bbMAD*1.4826 + 1e-6)

				z := (z1 + z2 + z3) / 3.0 // combined Z

				w := kernel * stats[j].Weight

				sum += z * w
				wsum += w
			}

			if wsum == 0 {
				continue
			}

			smoothed = append(smoothed, Smoothed{
				CHROM:     chrom,
				POS:       stats[i].POS,
				ZCombined: sum / wsum,
			})
		}
	}

	// -------------------------------------------------------------------------
	// Peak-based QTL detection
	// -------------------------------------------------------------------------

	var peaks []QTLPeak

	for chrom := range byChr {

		var chr []Smoothed
		for _, s := range smoothed {
			if s.CHROM == chrom {
				chr = append(chr, s)
			}
		}

		for i := 1; i < len(chr)-1; i++ {

			// Local maximum
			if chr[i].ZCombined > chr[i-1].ZCombined &&
				chr[i].ZCombined > chr[i+1].ZCombined &&
				chr[i].ZCombined > 3.0 {

				// Expand interval
				left := i
				for left > 0 && chr[left].ZCombined > 2.0 {
					left--
				}

				right := i
				for right < len(chr)-1 && chr[right].ZCombined > 2.0 {
					right++
				}

				peaks = append(peaks, QTLPeak{
					Chrom: chrom,
					Start: chr[left].POS,
					Stop:  chr[right].POS,
					Peak:  chr[i].ZCombined,
				})
			}
		}
	}

	return smoothed, peaks
}

func RunTwoBulkTwoParents(cfg utils.AnalysisConfig) {
	smoothed, peaks := RunTwoBulkV3(cfg.Rdr, cfg.HighBulkIdx, cfg.LowBulkIdx)
	fmt.Printf("TwoBulk: smoothed=%d peaks=%d\n", len(smoothed), len(peaks))
}

// -----------------------------------------------------------------------------
// QTL structure
// -----------------------------------------------------------------------------

type QTLPeak struct {
	Chrom string
	Start int64
	Stop  int64
	Peak  float64
}
