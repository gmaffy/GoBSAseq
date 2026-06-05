package filter

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/tabix"
	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/utils"
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

func isHomozygous(gt []int) bool {
	if len(gt) == 0 {
		return false
	}
	first := gt[0]
	if first < 0 {
		return false
	}
	for _, a := range gt[1:] {
		if a < 0 || a != first {
			return false
		}
	}
	return true
}

func getFloat(v *vcfgo.Variant, key string) (float64, bool) {
	raw, err := v.Info().Get(key)
	if err != nil || raw == nil {
		return 0, false
	}
	switch val := raw.(type) {
	case float32:
		return float64(val), true
	case float64:
		return val, true
	case int:
		return float64(val), true
	default:
		return 0, false
	}
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

func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
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
		return false
	}
}

func BsaSeqFilter(v *vcfgo.Variant, cfg utils.AnalysisConfig, bsaType string) bool {
	targetAlt := singleRealAlt(v)
	if targetAlt < 0 {
		return false
	}

	switch bsaType {
	case "2b":
		// bulks only
		for _, idx := range []int{cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		return v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "2p2b":
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
	case "2plb":
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
	case "2phb":
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
	case "hp2b":
		// high parent 2 bulks filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	case "lp2b":
		// low parent  2 bulks filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "hphb":
		// high parent high bulk filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth
	case "hplb":
		// high parent low bulk filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	case "lphb":
		// low parent high bulk filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth
	case "lplb":
		// low parent low bulk filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	default:
		return false
	}

}

type countingWriter struct {
	io.Writer
	n int64
}

func newTabixIndex() *tabix.Index {
	idx := tabix.New()
	idx.Format = 2 // VCF
	idx.NameColumn = 1
	idx.BeginColumn = 2
	idx.EndColumn = 0
	idx.MetaChar = '#'
	return idx
}

type vcfRecord struct {
	chrom string
	start int // 0-based
	end   int // 0-based half-open
}

func (r vcfRecord) RefName() string { return r.chrom }
func (r vcfRecord) Start() int      { return r.start }
func (r vcfRecord) End() int        { return r.end }

func addTabixRecord(idx *tabix.Index, v *vcfgo.Variant, chunk bgzf.Chunk) error {
	rec := vcfRecord{
		chrom: v.Chromosome,
		start: int(v.Pos) - 1,
		end:   int(v.Pos) - 1 + len(v.Ref()),
	}
	return idx.Add(rec, chunk, true, true)
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

func HardFilterVcf(cfg utils.AnalysisConfig, hfcfg utils.HardFilterConfig, bsaseqType string, keepIndices []int) ([]*vcfgo.Variant, int, int, error) {

	err := os.MkdirAll(filepath.Join(cfg.OutputDir, "stats"), 0775)
	if err != nil {
		return nil, 0, 0, err
	}

	hardFilteredVcfPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.hardfiltered.vcf.gz", bsaseqType))
	badVcfPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.lowqaul.vcf.gz", bsaseqType))

	rdr := cfg.Rdr
	origSampleNames := rdr.Header.SampleNames

	var newSampleNames []string
	for _, idx := range keepIndices {
		if idx >= 0 && idx < len(origSampleNames) {
			newSampleNames = append(newSampleNames, origSampleNames[idx])
		}
	}

	writerHeader := vcfgo.Header{
		SampleNames:   newSampleNames,
		Infos:         rdr.Header.Infos,
		SampleFormats: make(map[string]*vcfgo.SampleFormat, len(rdr.Header.SampleFormats)),
		Filters:       rdr.Header.Filters,
		Extras:        rdr.Header.Extras,
		FileFormat:    rdr.Header.FileFormat,
		Contigs:       rdr.Header.Contigs,
		Samples:       rdr.Header.Samples,
		Pedigrees:     rdr.Header.Pedigrees,
	}
	for id, format := range rdr.Header.SampleFormats {
		writerHeader.SampleFormats[id] = format
	}

	for _, id := range []string{"PGT", "PID"} {
		delete(writerHeader.SampleFormats, id)
	}

	// ── Open output files ─────────────────────────────────────────────────────
	hfFile, err := os.Create(hardFilteredVcfPath)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("create hard-filtered VCF: %w", err)
	}
	defer hfFile.Close()

	hfCounting := &countingWriter{Writer: hfFile}
	hfBgzf := bgzf.NewWriter(hfCounting, 1)
	hfWriter, err := vcfgo.NewWriter(hfBgzf, &writerHeader)
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
	badWriter, err := vcfgo.NewWriter(badBgzf, &writerHeader)
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
