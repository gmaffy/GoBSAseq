package twobulk

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"sync"
	"time"
)

type ThresholdsResult struct {
	Alphas []float64
	Opt1   map[float64]Thresholds            // [alpha]
	Opt2   map[string]map[float64]Thresholds // [chrom][alpha]
	Opt4   map[float64]Thresholds            // [alpha]
	Opt5   map[float64]Thresholds            // [alpha]

	// Option 3 memoization
	opt3Cache sync.Map // map[string]map[float64]Thresholds
	nReps     int
}

func calcThresholdsMulti(chromStats map[string][]BSAstats, allSmoothed []SmoothedStats, windowSize, stepSize int64, nReps int, alphas []float64) *ThresholdsResult {
	res := &ThresholdsResult{
		Alphas: alphas,
		Opt1:   make(map[float64]Thresholds),
		Opt2:   make(map[string]map[float64]Thresholds),
		Opt4:   make(map[float64]Thresholds),
		Opt5:   make(map[float64]Thresholds),
		nReps:  nReps,
	}

	// -------------------------------------------------------------------------
	// Option 5: Z-Score (Mean + Z*SD) excluding Chr "0"
	// -------------------------------------------------------------------------
	var sumD, sumG, sumE, sumL, sumB float64
	var count float64
	var valsD, valsG, valsE, valsL, valsB []float64

	for _, sm := range allSmoothed {
		if sm.CHROM == "0" || sm.CHROM == "chr0" || sm.CHROM == "Un" {
			continue
		}
		sumD += math.Abs(sm.DeltaSI)
		sumG += sm.Gstat
		sumE += sm.ED
		sumL += sm.LOD
		sumB += sm.BBLogBF

		valsD = append(valsD, math.Abs(sm.DeltaSI))
		valsG = append(valsG, sm.Gstat)
		valsE = append(valsE, sm.ED)
		valsL = append(valsL, sm.LOD)
		valsB = append(valsB, sm.BBLogBF)
		count++
	}

	if count > 0 {
		meanD, meanG, meanE, meanL, meanB := sumD/count, sumG/count, sumE/count, sumL/count, sumB/count
		var sqD, sqG, sqE, sqL, sqB float64
		for i := 0; i < int(count); i++ {
			sqD += (valsD[i] - meanD) * (valsD[i] - meanD)
			sqG += (valsG[i] - meanG) * (valsG[i] - meanG)
			sqE += (valsE[i] - meanE) * (valsE[i] - meanE)
			sqL += (valsL[i] - meanL) * (valsL[i] - meanL)
			sqB += (valsB[i] - meanB) * (valsB[i] - meanB)
		}
		sdD, sdG, sdE, sdL, sdB := math.Sqrt(sqD/count), math.Sqrt(sqG/count), math.Sqrt(sqE/count), math.Sqrt(sqL/count), math.Sqrt(sqB/count)

		for _, alpha := range alphas {
			// Z-score approximation for Normal distribution one-tailed
			z := 0.0
			if alpha == 0.05 {
				z = 1.645 // 95th percentile
			} else if alpha == 0.01 {
				z = 2.326 // 99th percentile
			} else {
				// crude approximation for other alphas
				z = 2.0
			}

			res.Opt5[alpha] = Thresholds{
				DeltaSI: meanD + z*sdD,
				Gstat:   meanG + z*sdG,
				ED:      meanE + z*sdE,
				LOD:     meanL + z*sdL,
				BBLogBF: meanB + z*sdB,
			}
		}
	}

	// -------------------------------------------------------------------------
	// Option 1, 2, 4: Permutation Testing
	// -------------------------------------------------------------------------
	type maxima struct {
		deltaSI, gstat, ed, lod, bbLogBF float64
	}

	var allStats []BSAstats
	for _, s := range chromStats {
		allStats = append(allStats, s...)
	}

	// For Option 1: genome-wide max excluding Chr 0
	opt1Results := make([]maxima, nReps)
	// For Option 2: max per chromosome
	opt2Results := make(map[string][]maxima)
	for chrom := range chromStats {
		opt2Results[chrom] = make([]maxima, nReps)
	}
	// For Option 4: all permuted window values to compute FDR empirical p-values
	var opt4AllD, opt4AllG, opt4AllE, opt4AllL, opt4AllB []float64
	var opt4Mu sync.Mutex

	repChan := make(chan int, nReps)
	for r := 0; r < nReps; r++ {
		repChan <- r
	}
	close(repChan)

	var wg sync.WaitGroup
	// 8 workers
	for w := 0; w < 8; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			rng := rand.New(rand.NewSource(time.Now().UnixNano()))
			for r := range repChan {
				permAll := permuteSNPs(allStats, rng)

				permChrom := make(map[string][]BSAstats)
				for _, s := range permAll {
					permChrom[s.CHROM] = append(permChrom[s.CHROM], s)
				}

				var mxOpt1 maxima
				chrMax := make(map[string]maxima)
				var localD, localG, localE, localL, localB []float64

				for chrom, stats := range permChrom {
					var mxChr maxima
					for _, sm := range smoothChromosome(stats, windowSize, stepSize) {
						ad := math.Abs(sm.DeltaSI)
						if ad > mxChr.deltaSI {
							mxChr.deltaSI = ad
						}
						if sm.Gstat > mxChr.gstat {
							mxChr.gstat = sm.Gstat
						}
						if sm.ED > mxChr.ed {
							mxChr.ed = sm.ED
						}
						if sm.LOD > mxChr.lod {
							mxChr.lod = sm.LOD
						}
						if sm.BBLogBF > mxChr.bbLogBF {
							mxChr.bbLogBF = sm.BBLogBF
						}

						localD = append(localD, ad)
						localG = append(localG, sm.Gstat)
						localE = append(localE, sm.ED)
						localL = append(localL, sm.LOD)
						localB = append(localB, sm.BBLogBF)
					}
					chrMax[chrom] = mxChr

					// Option 1 excludes Chr 0
					if chrom != "0" && chrom != "chr0" && chrom != "Un" {
						if mxChr.deltaSI > mxOpt1.deltaSI {
							mxOpt1.deltaSI = mxChr.deltaSI
						}
						if mxChr.gstat > mxOpt1.gstat {
							mxOpt1.gstat = mxChr.gstat
						}
						if mxChr.ed > mxOpt1.ed {
							mxOpt1.ed = mxChr.ed
						}
						if mxChr.lod > mxOpt1.lod {
							mxOpt1.lod = mxChr.lod
						}
						if mxChr.bbLogBF > mxOpt1.bbLogBF {
							mxOpt1.bbLogBF = mxChr.bbLogBF
						}
					}
				}

				opt1Results[r] = mxOpt1
				opt4Mu.Lock()
				for chrom, mx := range chrMax {
					if resArr, ok := opt2Results[chrom]; ok {
						resArr[r] = mx
					}
				}
				opt4AllD = append(opt4AllD, localD...)
				opt4AllG = append(opt4AllG, localG...)
				opt4AllE = append(opt4AllE, localE...)
				opt4AllL = append(opt4AllL, localL...)
				opt4AllB = append(opt4AllB, localB...)
				opt4Mu.Unlock()
			}
		}()
	}
	wg.Wait()

	// Compute quantiles for Option 1
	var o1d, o1g, o1e, o1l, o1b []float64
	for _, m := range opt1Results {
		o1d = append(o1d, m.deltaSI)
		o1g = append(o1g, m.gstat)
		o1e = append(o1e, m.ed)
		o1l = append(o1l, m.lod)
		o1b = append(o1b, m.bbLogBF)
	}
	sort.Float64s(o1d)
	sort.Float64s(o1g)
	sort.Float64s(o1e)
	sort.Float64s(o1l)
	sort.Float64s(o1b)
	for _, alpha := range alphas {
		q := 1.0 - alpha
		res.Opt1[alpha] = Thresholds{
			DeltaSI: quantile(o1d, q),
			Gstat:   quantile(o1g, q),
			ED:      quantile(o1e, q),
			LOD:     quantile(o1l, q),
			BBLogBF: quantile(o1b, q),
		}
	}

	// Compute quantiles for Option 2
	for chrom, maxArr := range opt2Results {
		res.Opt2[chrom] = make(map[float64]Thresholds)
		var o2d, o2g, o2e, o2l, o2b []float64
		for _, m := range maxArr {
			o2d = append(o2d, m.deltaSI)
			o2g = append(o2g, m.gstat)
			o2e = append(o2e, m.ed)
			o2l = append(o2l, m.lod)
			o2b = append(o2b, m.bbLogBF)
		}
		sort.Float64s(o2d)
		sort.Float64s(o2g)
		sort.Float64s(o2e)
		sort.Float64s(o2l)
		sort.Float64s(o2b)
		for _, alpha := range alphas {
			q := 1.0 - alpha
			res.Opt2[chrom][alpha] = Thresholds{
				DeltaSI: quantile(o2d, q),
				Gstat:   quantile(o2g, q),
				ED:      quantile(o2e, q),
				LOD:     quantile(o2l, q),
				BBLogBF: quantile(o2b, q),
			}
		}
	}

	// Compute FDR thresholds for Option 4
	sort.Float64s(opt4AllD)
	sort.Float64s(opt4AllG)
	sort.Float64s(opt4AllE)
	sort.Float64s(opt4AllL)
	sort.Float64s(opt4AllB)

	// Since we want the threshold where FDR < alpha:
	// FDR for a stat value v is approx: (number of permuted stats >= v) / (total permuted) * (total real) / (number of real stats >= v)
	// For simplicity and speed, we will approximate the FDR threshold by finding the value at the (1 - (alpha/TotalRealWindows)) quantile of the null, or applying BH.
	// Actually, applying BH to find a single threshold cutoff across the genome:
	// 1. Sort real values descending.
	// 2. For each real value, find its p-value = count(perm >= real) / len(perm).
	// 3. BH adjusted p-value = p-value * len(real) / rank.
	// 4. Find the largest rank where BH adjusted p-value < alpha.
	// 5. The real value at this rank is our threshold.

	findFDRThreshold := func(realVals []float64, permVals []float64, alpha float64) float64 {
		sort.Sort(sort.Reverse(sort.Float64Slice(realVals)))
		nPerm := float64(len(permVals))
		nReal := float64(len(realVals))

		for rank, rv := range realVals {
			// binary search for count(perm >= rv)
			idx := sort.Search(len(permVals), func(i int) bool { return permVals[i] >= rv })
			countGreater := float64(len(permVals) - idx)
			pVal := countGreater / nPerm
			fdr := pVal * nReal / float64(rank+1)

			if fdr > alpha {
				// Since we are iterating from largest stat (smallest pval) downwards,
				// the first time FDR exceeds alpha, the PREVIOUS value is our threshold.
				if rank == 0 {
					return rv // none passed, return max
				}
				return realVals[rank-1]
			}
		}
		if len(realVals) > 0 {
			return realVals[len(realVals)-1]
		}
		return 0
	}

	for _, alpha := range alphas {
		res.Opt4[alpha] = Thresholds{
			DeltaSI: findFDRThreshold(valsD, opt4AllD, alpha),
			Gstat:   findFDRThreshold(valsG, opt4AllG, alpha),
			ED:      findFDRThreshold(valsE, opt4AllE, alpha),
			LOD:     findFDRThreshold(valsL, opt4AllL, alpha),
			BBLogBF: findFDRThreshold(valsB, opt4AllB, alpha),
		}
	}

	return res
}

// getOpt3Threshold simulates depth CIs on the fly (with memoization)
func (tr *ThresholdsResult) getOpt3Threshold(meanHighDP, meanLowDP int) map[float64]Thresholds {
	key := fmt.Sprintf("%d_%d", meanHighDP, meanLowDP)
	if val, ok := tr.opt3Cache.Load(key); ok {
		return val.(map[float64]Thresholds)
	}

	res := make(map[float64]Thresholds)
	if meanHighDP == 0 || meanLowDP == 0 {
		for _, alpha := range tr.Alphas {
			res[alpha] = Thresholds{}
		}
		tr.opt3Cache.Store(key, res)
		return res
	}

	nReps := 10000 // typical for QTL-seq simulation
	var ds, gs, es, ls, bs []float64

	for i := 0; i < nReps; i++ {
		// simulate Binomial under null p=0.5
		var hbH, lbH int
		for d := 0; d < meanHighDP; d++ {
			if rand.Float32() < 0.5 {
				hbH++
			}
		}
		for d := 0; d < meanLowDP; d++ {
			if rand.Float32() < 0.5 {
				lbH++
			}
		}
		hbL := meanHighDP - hbH
		lbL := meanLowDP - lbH

		hSI := float64(hbH) / float64(meanHighDP)
		lSI := float64(lbH) / float64(meanLowDP)

		ds = append(ds, math.Abs(hSI-lSI))
		gs = append(gs, GStatistic(hbH, hbL, lbH, lbL))
		es = append(es, math.Abs(hSI-lSI))
		ls = append(ls, lod(hbL, hbH, lbL, lbH))
		bs = append(bs, betaBinomialLogBF(hbH, hbL, lbH, lbL))
	}

	sort.Float64s(ds)
	sort.Float64s(gs)
	sort.Float64s(es)
	sort.Float64s(ls)
	sort.Float64s(bs)

	for _, alpha := range tr.Alphas {
		q := 1.0 - alpha
		res[alpha] = Thresholds{
			DeltaSI: quantile(ds, q),
			Gstat:   quantile(gs, q),
			ED:      quantile(es, q),
			LOD:     quantile(ls, q),
			BBLogBF: quantile(bs, q),
		}
	}

	tr.opt3Cache.Store(key, res)
	return res
}
