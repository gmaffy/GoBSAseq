package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// BRMBlock is a BRM block: [Start, Stop] where |signal| ≥ u·σ, and its peak
// (the largest |signal|): position, signal, the signal threshold u·σ there
// and the standardised z = signal/σ, which is what is compared with u.
type BRMBlock struct {
	Chrom     string
	Start     int64
	Stop      int64
	PeakPos   int64
	Peak      float64
	Threshold float64
	PeakZ     float64
	NMarkers  int // markers beyond the line (--min-qtl-markers)
}

// brmModel is a mode's BRM test at a marker: the signal and its null
// variance. BRM calls markers where |signal| ≥ u·σ, i.e. |z| ≥ u for the
// standardised z = signal/σ; u comes from the genome-wide null (BRMU in
// Thresholds), so the error rate --brm-alpha is genome-wide.
//
// Signals: the smoothed signed ΔSI (two-bulk, polarity known); √max(smoothed
// bias-corrected ΔSI², 0), a polarity-free |ΔAF| estimate (2b without
// --ref-parent; one-sided); the smoothed AF deviation SmAFDev (one-bulk, from
// the bias-corrected P0).
//
// Null variance σ² = pool variance + smoothed read variance:
//   - pool: p(1−p)/(V_scale·n) per bulk, p the null frequency (p0; one-bulk:
//     the site's P0), V_scale the population variance scale (2 for F2), n the
//     bulk size. It is a single position's variance, a slight upper bound for
//     its kernel average (linked pool frequencies vary little in the kernel).
//   - read: ReadVar0 + ρ·ReadVar1 from the kernel pass, the exact variance of
//     the smoothed value under beta-binomial reads with overdispersion ρ
//     (RunBRM requires it).
type brmModel struct {
	signal    func(*SmoothedStats) float64
	sigma2    func(s *SmoothedStats, rho float64) float64
	polarFree bool // 2b without --ref-parent: the signal is √(ΔSI²c) ≥ 0
}

func newBRMModel(cfg utils.AnalysisConfig, bsaType string) (*brmModel, error) {
	_, _, both, one := BulkFlags(bsaType)
	if !both && !one {
		return nil, fmt.Errorf("BRM: mode %q has no bulk", bsaType)
	}
	high, low, oneN, err := cfg.BulkSizes(bsaType)
	if err != nil {
		return nil, err
	}
	pop, err := utils.ParsePopulation(cfg.Population)
	if err != nil {
		return nil, err
	}
	p0, c := pop.P0, pop.VarScale
	read := func(s *SmoothedStats, rho float64) float64 { return s.ReadVar0 + rho*s.ReadVar1 }
	m := &brmModel{}
	if both {
		v := p0 * (1 - p0)
		pool := v/(c*float64(high)) + v/(c*float64(low))
		m.signal = func(s *SmoothedStats) float64 { return s.SmDeltaSI }
		if !cfg.SignedDeltaSI(bsaType) {
			m.signal = func(s *SmoothedStats) float64 { return math.Sqrt(math.Max(s.SmDeltaSI2, 0)) }
			m.polarFree = true
		}
		m.sigma2 = func(s *SmoothedStats, rho float64) float64 { return pool + read(s, rho) }
		return m, nil
	}
	pAt := func(s *SmoothedStats) float64 {
		if s.OneBulkP0 > 0 && s.OneBulkP0 < 1 {
			return s.OneBulkP0
		}
		return p0
	}
	m.signal = func(s *SmoothedStats) float64 { return s.SmAFDev }
	m.sigma2 = func(s *SmoothedStats, rho float64) float64 {
		p := pAt(s)
		return p*(1-p)/(c*float64(oneN)) + read(s, rho)
	}
	return m, nil
}

// z is the standardised BRM statistic at a marker, |signal|/σ (the
// polarity-free signal is ≥ 0, so this is one-sided there).
func (m *brmModel) z(s *SmoothedStats, rho float64) float64 {
	sd := math.Sqrt(m.sigma2(s, rho))
	if !(sd > 0) {
		return 0
	}
	return math.Abs(m.signal(s)) / sd
}

// blocks scans one chromosome's markers (sorted by position) for runs with
// |signal| ≥ u·σ.
func (m *brmModel) blocks(chrom string, stats []SmoothedStats, rho, u float64) []BRMBlock {
	return scanBRMBlocks(chrom, stats, u, func(s *SmoothedStats) (float64, float64) {
		return m.signal(s), math.Sqrt(m.sigma2(s, rho))
	})
}

// scanBRMBlocks returns the runs of one chromosome's markers (sorted by
// position) where |value| ≥ u·σ, with value and σ from eval (scanCrossings:
// interpolated ends; the peak is the largest |value|).
func scanBRMBlocks(chrom string, stats []SmoothedStats, u float64, eval func(*SmoothedStats) (value, sigma float64)) []BRMBlock {
	n := len(stats)
	value, sigma := make([]float64, n), make([]float64, n)
	for i := range stats {
		value[i], sigma[i] = eval(&stats[i])
	}
	margin := func(i int) float64 {
		if !(u*sigma[i] > 0) {
			return math.Inf(-1)
		}
		return math.Abs(value[i]) - u*sigma[i]
	}
	var blocks []BRMBlock
	for _, c := range scanCrossings(n, func(i int) float64 { return float64(stats[i].POS) }, margin,
		func(i int) float64 { return math.Abs(value[i]) }) {
		p := c.peak
		blocks = append(blocks, BRMBlock{Chrom: chrom, Start: int64(math.Round(c.start)), Stop: int64(math.Round(c.end)),
			PeakPos: stats[p].POS, Peak: value[p], Threshold: u * sigma[p], PeakZ: value[p] / sigma[p], NMarkers: c.last - c.first + 1})
	}
	return blocks
}

// RunBRM runs BRM block detection across chromosomes (in sorted order) and
// writes a TSV of blocks. th supplies the genome-wide threshold BRMU (|z| ≥ u)
// and the read overdispersion Rho from the null simulation.
func RunBRM(cfg utils.AnalysisConfig, bsaType string, sm []SmoothedStats, th Thresholds) ([]BRMBlock, error) {

	color.Cyan("===================================== Running BRM block detection =========================================")
	if len(sm) == 0 {
		return nil, nil
	}
	if !(th.BRMU > 0) || math.IsInf(th.BRMU, 0) {
		return nil, fmt.Errorf("BRM: no genome-wide threshold (the null simulation must run first)")
	}
	m, err := newBRMModel(cfg, bsaType)
	if err != nil {
		return nil, err
	}
	for i := range sm {
		if !(sm[i].ReadVar0 > 0) {
			return nil, fmt.Errorf("BRM: smoothed marker %s:%d has no read variance (smooth with the scorer)", sm[i].CHROM, sm[i].POS)
		}
	}
	if m.polarFree {
		color.Cyan("BRM: polarity-free form (2b without --ref-parent): √(smoothed bias-corrected ΔSI²); block peaks have no direction")
	}
	color.Cyan("BRM: |z| ≥ %.3f (genome-wide α = %g from %d null genomes; ρ = %.4f)", th.BRMU, th.BRMAlpha, th.NullReps, th.Rho)

	// Chromosomes in sorted order (stable output), scanned in place; blocks
	// with fewer than --min-qtl-markers markers are dropped before merging.
	minMarkers := cfg.EffectiveMinQTLMarkers()
	var allBlocks []BRMBlock
	for _, c := range byChromosome(sm) {
		blocks := m.blocks(c.chrom, c.sm, th.Rho, th.BRMU)
		kept := blocks[:0]
		for _, b := range blocks {
			if b.NMarkers >= minMarkers {
				kept = append(kept, b)
			}
		}
		if len(kept) > 0 {
			allBlocks = append(allBlocks, mergeBRMBlocks(kept, cfg.EffectiveMergeDistance())...)
		}
	}

	outPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.brm_blocks.tsv", bsaType))
	if err := writeBRMBlocks(outPath, allBlocks); err != nil {
		return nil, fmt.Errorf("writing %s: %w", outPath, err)
	}
	color.Green("BRM blocks written to %s", outPath)
	return allBlocks, nil
}

// writeBRMBlocks writes the blocks TSV, returning any write, flush or close
// error.
func writeBRMBlocks(path string, blocks []BRMBlock) (err error) {
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	w := bufio.NewWriter(f)
	fmt.Fprintln(w, "CHROM\tSTART\tSTOP\tPEAK_POS\tPEAK\tTHRESHOLD\tPEAK_Z")
	for _, b := range blocks {
		fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.3f\n", b.Chrom, b.Start, b.Stop, b.PeakPos, b.Peak, b.Threshold, b.PeakZ)
	}
	return w.Flush()
}

// mergeBRMBlocks joins one chromosome's blocks (in position order, as
// scanBRMBlocks returns them) whose peaks have the same sign and which are at
// most gap bp apart (consecutive blocks only, as before); the merged block
// keeps the larger peak and counts all its markers.
func mergeBRMBlocks(blocks []BRMBlock, gap int64) []BRMBlock {
	return mergeRuns(blocks, gap, func(b *BRMBlock) (string, float64, float64) { return b.Chrom, float64(b.Start), float64(b.Stop) },
		func(a, b *BRMBlock) bool { return a.Peak == 0 || b.Peak == 0 || (a.Peak > 0) == (b.Peak > 0) },
		func(a *BRMBlock, b BRMBlock) {
			a.Stop = max(a.Stop, b.Stop)
			a.NMarkers += b.NMarkers
			if math.Abs(b.Peak) > math.Abs(a.Peak) {
				a.PeakPos, a.Peak, a.Threshold, a.PeakZ = b.PeakPos, b.Peak, b.Threshold, b.PeakZ
			}
		})
}
