package utils

import (
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/tabix"
	"github.com/brentp/vcfgo"
	"github.com/schollz/progressbar/v3"
)

func singleRealAlt(v *vcfgo.Variant) int {
	realIdx := -1
	count := 0
	for i, alt := range v.Alt() {
		if !(alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<')) {
			realIdx = i
			count++
		}
	}
	if count != 1 {
		return -1
	}
	return realIdx + 1 // vcfgo uses 1-based allele numbering
}

func classifyVariant(v *vcfgo.Variant) (isSNP, isIndel bool) {
	refLen := len(v.Ref())
	isSNP = refLen == 1
	for _, alt := range v.Alt() {
		if alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<') {
			continue
		}
		if len(alt) != 1 {
			isSNP = false
		}
		if len(alt) != refLen {
			isIndel = true
		}
	}
	return
}

func sampleHasOnlyRefOrAlt(v *vcfgo.Variant, idx, targetAlt int) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]
	if s == nil || len(s.GT) == 0 {
		return false
	}
	for _, allele := range s.GT {
		if allele < 0 { // missing ('.')
			return false
		}
		if allele != 0 && allele != targetAlt {
			return false
		}
	}
	return true
}

func PassesHardFilter(v *vcfgo.Variant, hfcfg HardFilterConfig) bool {
	isSNP, isIndel := classifyVariant(v)

	switch {
	case isSNP:
		if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
			return false
		}
		if qd, ok := getFloat(v, "QD"); ok && qd < hfcfg.SNP_QD_Min {
			return false
		}
		if fs, ok := getFloat(v, "FS"); ok && fs > hfcfg.SNP_FS_Max {
			return false
		}
		if sor, ok := getFloat(v, "SOR"); ok && sor > hfcfg.SNP_SOR_Max {
			return false
		}
		if mq, ok := getFloat(v, "MQ"); ok && mq < hfcfg.SNP_MQ_Min {
			return false
		}
		if mqrs, ok := getFloat(v, "MQRankSum"); ok && mqrs < hfcfg.SNP_MQRankSum_Min {
			return false
		}
		if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.SNP_ReadPosRankSum_Min {
			return false
		}
		return true

	case isIndel:
		if float64(v.Quality) < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if qd, ok := getFloat(v, "QD"); ok && qd < hfcfg.INDEL_QD_Min {
			return false
		}
		if fs, ok := getFloat(v, "FS"); ok && fs > hfcfg.INDEL_FS_Max {
			return false
		}
		if sor, ok := getFloat(v, "SOR"); ok && sor > hfcfg.INDEL_SOR_Max {
			return false
		}
		if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.INDEL_ReadPosRankSum_Min {
			return false
		}
		return true

	default:
		return false // mixed / MNP / symbolic-only
	}
}

func BsaSeqFilter(v *vcfgo.Variant, cfg AnalysisConfig, bsaType int) bool {
	targetAlt := singleRealAlt(v)
	if targetAlt < 0 {
		return false
	}

	switch bsaType {
	case 0:
		// bulks only
		for _, idx := range []int{cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		return v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case 1:
		// 2 parents 2 bulks filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]

		if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) {
			return false
		}
		if hp.GT[0] == lp.GT[0] {
			return false // same allele in both parents — not informative
		}

		return hp.DP >= cfg.HighParentDepth &&
			lp.DP >= cfg.LowParentDepth &&
			v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth &&
			v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	case 2:
		// 2 parents low bulk filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]

		if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) {
			return false
		}
		if hp.GT[0] == lp.GT[0] {
			return false // same allele in both parents — not informative
		}

		return hp.DP >= cfg.HighParentDepth && lp.DP >= cfg.LowParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	case 3:
		// 2 parents high bulk filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]

		if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) {
			return false
		}
		if hp.GT[0] == lp.GT[0] {
			return false // same allele in both parents — not informative
		}

		return hp.DP >= cfg.HighParentDepth && lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth
	case 4:
		// high parent 2 bulks filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	case 5:
		// low parent  2 bulks filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	default:
		return false
	}

}

// ---------------------------------------------------------------------------
// countingWriter wraps an io.Writer and tracks total bytes written.
// Required to compute bgzf.Offset.File (the compressed file byte position)
// when building the tabix index, since bgzf.Writer does not expose this.
// ---------------------------------------------------------------------------

type countingWriter struct {
	io.Writer
	n int64
}

func (w *countingWriter) Write(p []byte) (int, error) {
	n, err := w.Writer.Write(p)
	w.n += int64(n)
	return n, err
}

func newWriterWithSampleNames(w io.Writer, h *vcfgo.Header, sampleNames []string) (*vcfgo.Writer, error) {
	originalSampleNames := h.SampleNames
	h.SampleNames = sampleNames
	writer, err := vcfgo.NewWriter(w, h)
	h.SampleNames = originalSampleNames
	return writer, err
}

// ---------------------------------------------------------------------------
// HardFilterVcf
// ---------------------------------------------------------------------------

func getKeepIndices(cfg AnalysisConfig, bsaType int) []int {
	var idxs []int
	switch bsaType {
	case 0:
		idxs = []int{cfg.HighBulkIdx, cfg.LowBulkIdx}
	case 1:
		idxs = []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}
	case 2:
		idxs = []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.LowBulkIdx}
	case 3:
		idxs = []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx}
	case 4:
		idxs = []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}
	case 5:
		idxs = []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}
	}

	var kept []int
	for _, idx := range idxs {
		if idx >= 0 {
			kept = append(kept, idx)
		}
	}
	return kept
}

func sanitizeVariant(v *vcfgo.Variant, keepIndices []int) {
	// Subset samples
	newSamples := make([]*vcfgo.SampleGenotype, len(keepIndices))
	for i, idx := range keepIndices {
		if idx < len(v.Samples) {
			newSamples[i] = v.Samples[idx]
		}
	}
	v.Samples = newSamples

	// Clean Format string - remove PGT, PID which are being deleted from header
	newFormat := make([]string, 0, len(v.Format))
	for _, f := range v.Format {
		if f != "PGT" && f != "PID" {
			newFormat = append(newFormat, f)
		}
	}
	v.Format = newFormat
}

func HardFilterVcf(rdr *vcfgo.Reader, hardFilteredVcfPath string, badVcfPath string, cfg AnalysisConfig, hfcfg HardFilterConfig, bsaseqType int) ([]*vcfgo.Variant, int, int, error) {

	// Get original sample names
	origSampleNames := rdr.Header.SampleNames

	// Subset samples in header for output
	keepIndices := getKeepIndices(cfg, bsaseqType)
	var newSampleNames []string
	for _, idx := range keepIndices {
		if idx >= 0 && idx < len(origSampleNames) {
			newSampleNames = append(newSampleNames, origSampleNames[idx])
		}
	}

	for _, id := range []string{"PGT", "PID"} {
		delete(rdr.Header.SampleFormats, id)
	}

	// ── Open output files ─────────────────────────────────────────────────────
	hfFile, err := os.Create(hardFilteredVcfPath)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("create hard-filtered VCF: %w", err)
	}
	defer hfFile.Close()

	hfCounting := &countingWriter{Writer: hfFile}
	hfBgzf := bgzf.NewWriter(hfCounting, 1)
	hfWriter, err := newWriterWithSampleNames(hfBgzf, rdr.Header, newSampleNames)
	if err != nil {
		hfBgzf.Close()
		return nil, 0, 0, fmt.Errorf("create hard-filtered VCF writer: %w", err)
	}

	badFile, err := os.Create(badVcfPath)
	if err != nil {
		hfBgzf.Close()
		return nil, 0, 0, fmt.Errorf("create rejected VCF: %w", err)
	}
	defer badFile.Close()

	badBgzf := bgzf.NewWriter(badFile, 1)
	badWriter, err := newWriterWithSampleNames(badBgzf, rdr.Header, newSampleNames)
	if err != nil {
		hfBgzf.Close()
		badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("create rejected VCF writer: %w", err)
	}

	hfIdx := newTabixIndex()
	bar := progressbar.Default(-1, "Hard filtering variants")

	// ── Pipeline channels ─────────────────────────────────────────────────────
	type variantResult struct {
		v      *vcfgo.Variant
		passed bool
	}

	const chanBuf = 512
	filterCh := make(chan *vcfgo.Variant, chanBuf) // reader  → filter workers
	resultCh := make(chan variantResult, chanBuf)  // workers → writer
	var readerErr atomic.Pointer[error]

	// ── Stage 1: Reader goroutine ─────────────────────────────────────────────
	var (
		originalCount atomic.Int64
		skippedCount  atomic.Int64
	)

	go func() {
		defer close(filterCh)
		for {
			v := rdr.Read()
			if v == nil {
				break
			}
			if err := rdr.Error(); err != nil {
				if strings.Contains(err.Error(), "bad sample string") {
					rdr.Clear()
				} else {
					e := fmt.Errorf("VCF parse error at line %d: %w", v.LineNumber, err)
					readerErr.Store(&e)
					return
				}
			}
			alts := v.Alt()
			if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
				skippedCount.Add(1)
				continue
			}
			originalCount.Add(1)
			_ = bar.Add(1)
			filterCh <- v
		}

		if err := rdr.Error(); err != nil && !strings.Contains(err.Error(), "bad sample string") {
			e := fmt.Errorf("VCF read error: %w", err)
			readerErr.Store(&e)
		}
	}()

	// ── Stage 2: Filter worker pool ───────────────────────────────────────────
	// Filtering is CPU-bound (no shared state); run one worker per core.
	// vcfgo.Variant fields accessed here are read-only after Read() returns,
	// so no mutex is needed inside the workers.
	numWorkers := runtime.GOMAXPROCS(0)
	var wg sync.WaitGroup
	wg.Add(numWorkers)
	for range numWorkers {
		go func() {
			defer wg.Done()
			for v := range filterCh {
				passed := PassesHardFilter(v, hfcfg) && BsaSeqFilter(v, cfg, bsaseqType)
				resultCh <- variantResult{v: v, passed: passed}
			}
		}()
	}

	// Close resultCh once all workers are done.
	go func() {
		wg.Wait()
		close(resultCh)
	}()

	// ── Stage 3: Writer goroutine (single, preserves order not required) ──────
	// NOTE: bgzf / vcfgo writers are NOT goroutine-safe; keep all writes here.
	var (
		passedVariants []*vcfgo.Variant
		passed         int
		writerErr      error
	)

	for res := range resultCh {
		// Sanitize variant (subset samples and clean format)
		sanitizeVariant(res.v, keepIndices)

		if !res.passed {
			badWriter.WriteVariant(res.v)
			continue
		}

		blockOffset, _ := hfBgzf.Next()
		startOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffset)}

		hfWriter.WriteVariant(res.v)

		blockOffsetEnd, _ := hfBgzf.Next()
		endOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffsetEnd)}

		if err := addTabixRecord(hfIdx, res.v, bgzf.Chunk{Begin: startOffset, End: endOffset}); err != nil {
			writerErr = fmt.Errorf("tabix add variant at %s:%d: %w", res.v.Chromosome, res.v.Pos, err)
			break
		}

		passedVariants = append(passedVariants, res.v)
		passed++
	}

	// Drain resultCh if writer broke early, so workers can unblock and exit.
	for range resultCh {
	}

	// ── Check errors ──────────────────────────────────────────────────────────
	if ep := readerErr.Load(); ep != nil {
		hfBgzf.Close()
		badBgzf.Close()
		return nil, 0, 0, *ep
	}
	if writerErr != nil {
		hfBgzf.Close()
		badBgzf.Close()
		return nil, 0, 0, writerErr
	}

	_ = bar.Finish()

	// ── Flush and close bgzf streams ──────────────────────────────────────────
	if err = hfBgzf.Close(); err != nil {
		badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("close hard-filtered bgzf: %w", err)
	}
	if err = badBgzf.Close(); err != nil {
		return nil, 0, 0, fmt.Errorf("close rejected bgzf: %w", err)
	}

	// ── Write tabix index ─────────────────────────────────────────────────────
	tbiFile, err := os.Create(hardFilteredVcfPath + ".tbi")
	if err != nil {
		return nil, 0, 0, fmt.Errorf("create tbi file: %w", err)
	}
	defer tbiFile.Close()

	tbiGz := bgzf.NewWriter(tbiFile, 1)
	if err := tabix.WriteTo(tbiGz, hfIdx); err != nil {
		tbiGz.Close()
		return nil, 0, 0, fmt.Errorf("write tabix index: %w", err)
	}
	if err := tbiGz.Close(); err != nil {
		return nil, 0, 0, fmt.Errorf("close tbi bgzf: %w", err)
	}

	original := int(originalCount.Load())
	skipped := int(skippedCount.Load())

	fmt.Printf(
		"Hard filtering complete: %d variant records read (%d gVCF ref blocks skipped) → %d passed, %d rejected\n",
		original, skipped, passed, original-passed,
	)

	return passedVariants, original, passed, nil
}
