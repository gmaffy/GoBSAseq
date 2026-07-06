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
	"github.com/fatih/color"
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

// sampleAllelesOKLenient checks that every non-missing allele in the
// sample's genotype is either ref (0) or targetAlt, ignoring missing calls
// ('.'). Unlike sampleHasOnlyRefOrAlt, a partially-missing genotype is not
// an automatic failure — only a fully-missing one is, since there's no
// signal left to use at all.
func sampleAllelesOKLenient(v *vcfgo.Variant, idx, targetAlt int) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]
	if s == nil || len(s.GT) == 0 {
		return false
	}
	seenCall := false
	for _, allele := range s.GT {
		if allele < 0 {
			continue // missing — tolerated
		}
		seenCall = true
		if allele != 0 && allele != targetAlt {
			return false
		}
	}
	return seenCall
}

// parentAllele resolves which single allele (ref=0 or targetAlt) a parent
// carries, tolerating missing calls within its own genotype (e.g. a diploid
// call dropping to effectively haploid under low coverage). It fails if the
// parent has no usable genotype at all, if its non-missing alleles disagree
// with each other (i.e. it's genuinely heterozygous, not just partially
// missing), or if a called allele is neither ref nor targetAlt.
func parentAllele(v *vcfgo.Variant, idx, targetAlt int) (allele int, ok bool) {
	if idx < 0 || idx >= len(v.Samples) || v.Samples[idx] == nil {
		return -1, false
	}
	resolved := -1
	for _, a := range v.Samples[idx].GT {
		if a < 0 {
			continue // missing — tolerated
		}
		if a != 0 && a != targetAlt {
			return -1, false
		}
		if resolved == -1 {
			resolved = a
		} else if a != resolved {
			return -1, false // genuinely heterozygous
		}
	}
	if resolved == -1 {
		return -1, false // fully missing — can't determine
	}
	return resolved, true
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
	if hfcfg.LightFilter {
		return true
	}
	switch {
	case isSNP:
		if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
			return false
		}
		if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.SNP_QD_Min {
			return false
		}
		if fs, ok := utils.GetFloat(v, "FS"); ok && fs > hfcfg.SNP_FS_Max {
			return false
		}
		if sor, ok := utils.GetFloat(v, "SOR"); ok && sor > hfcfg.SNP_SOR_Max {
			return false
		}
		if mq, ok := utils.GetFloat(v, "MQ"); ok && mq < hfcfg.SNP_MQ_Min {
			return false
		}
		if mqrs, ok := utils.GetFloat(v, "MQRankSum"); ok && mqrs < hfcfg.SNP_MQRankSum_Min {
			return false
		}
		if rprs, ok := utils.GetFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.SNP_ReadPosRankSum_Min {
			return false
		}
		return true

	case isIndel:
		if float64(v.Quality) < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.INDEL_QD_Min {
			return false
		}
		if fs, ok := utils.GetFloat(v, "FS"); ok && fs > hfcfg.INDEL_FS_Max {
			return false
		}
		if sor, ok := utils.GetFloat(v, "SOR"); ok && sor > hfcfg.INDEL_SOR_Max {
			return false
		}
		if rprs, ok := utils.GetFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.INDEL_ReadPosRankSum_Min {
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
		// bulks only — check depths only, genotype calls in bulks are ignored
		if cfg.HighBulkIdx < 0 || cfg.HighBulkIdx >= len(v.Samples) || v.Samples[cfg.HighBulkIdx] == nil {
			return false
		}
		if cfg.LowBulkIdx < 0 || cfg.LowBulkIdx >= len(v.Samples) || v.Samples[cfg.LowBulkIdx] == nil {
			return false
		}
		return v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "2p2b":
		// parents may have partial missing data. If one parent is completely missing,
		// we infer/rescue its allele as the opposite of the resolved parent.
		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		if !hpOK && !lpOK {
			return false // both parents missing — cannot resolve
		}

		if hpOK && lpOK {
			if hpAllele == lpAllele {
				return false // same allele in both parents — not informative
			}
		} else {
			// One parent is resolved, the other is missing.
			// Infer the missing parent has the opposite allele since the site is biallelic.
			if hpOK {
				if hpAllele == 0 {
					lpAllele = targetAlt
				} else {
					lpAllele = 0
				}
			} else {
				if lpAllele == 0 {
					hpAllele = targetAlt
				} else {
					hpAllele = 0
				}
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]
		hb := v.Samples[cfg.HighBulkIdx]
		lb := v.Samples[cfg.LowBulkIdx]

		if hp == nil || lp == nil || hb == nil || lb == nil {
			return false
		}

		// Enforce depth checks only on the parent that was actually resolved
		hpDepthPass := !hpOK || hp.DP >= cfg.HighParentDepth
		lpDepthPass := !lpOK || lp.DP >= cfg.LowParentDepth

		return hpDepthPass && lpDepthPass && hb.DP >= cfg.HighBulkDepth && lb.DP >= cfg.LowBulkDepth

	case "2plb":
		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		if !hpOK && !lpOK {
			return false
		}

		if hpOK && lpOK {
			if hpAllele == lpAllele {
				return false
			}
		} else {
			if hpOK {
				if hpAllele == 0 {
					lpAllele = targetAlt
				} else {
					lpAllele = 0
				}
			} else {
				if lpAllele == 0 {
					hpAllele = targetAlt
				} else {
					hpAllele = 0
				}
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]
		lb := v.Samples[cfg.LowBulkIdx]

		if hp == nil || lp == nil || lb == nil {
			return false
		}

		hpDepthPass := !hpOK || hp.DP >= cfg.HighParentDepth
		lpDepthPass := !lpOK || lp.DP >= cfg.LowParentDepth

		return hpDepthPass && lpDepthPass && lb.DP >= cfg.LowBulkDepth

	case "2phb":
		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		if !hpOK && !lpOK {
			return false
		}

		if hpOK && lpOK {
			if hpAllele == lpAllele {
				return false
			}
		} else {
			if hpOK {
				if hpAllele == 0 {
					lpAllele = targetAlt
				} else {
					lpAllele = 0
				}
			} else {
				if lpAllele == 0 {
					hpAllele = targetAlt
				} else {
					hpAllele = 0
				}
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]
		hb := v.Samples[cfg.HighBulkIdx]

		if hp == nil || lp == nil || hb == nil {
			return false
		}

		hpDepthPass := !hpOK || hp.DP >= cfg.HighParentDepth
		lpDepthPass := !lpOK || lp.DP >= cfg.LowParentDepth

		return hpDepthPass && lpDepthPass && hb.DP >= cfg.HighBulkDepth

	case "hp2b":
		// parent may have partial missing data
		if !sampleAllelesOKLenient(v, cfg.HighParentIdx, targetAlt) {
			return false
		}

		hp := v.Samples[cfg.HighParentIdx]
		hb := v.Samples[cfg.HighBulkIdx]
		lb := v.Samples[cfg.LowBulkIdx]
		if hp == nil || hb == nil || lb == nil {
			return false
		}
		return hp.DP >= cfg.HighParentDepth && hb.DP >= cfg.HighBulkDepth && lb.DP >= cfg.LowBulkDepth

	case "lp2b":
		// parent may have partial missing data
		if !sampleAllelesOKLenient(v, cfg.LowParentIdx, targetAlt) {
			return false
		}
		lp := v.Samples[cfg.LowParentIdx]
		hb := v.Samples[cfg.HighBulkIdx]
		lb := v.Samples[cfg.LowBulkIdx]
		if lp == nil || hb == nil || lb == nil {
			return false
		}
		return lp.DP >= cfg.LowParentDepth && hb.DP >= cfg.HighBulkDepth && lb.DP >= cfg.LowBulkDepth

	case "hphb":
		// parent may have partial missing data
		if !sampleAllelesOKLenient(v, cfg.HighParentIdx, targetAlt) {
			return false
		}
		hp := v.Samples[cfg.HighParentIdx]
		hb := v.Samples[cfg.HighBulkIdx]
		if hp == nil || hb == nil {
			return false
		}
		return hp.DP >= cfg.HighParentDepth && hb.DP >= cfg.HighBulkDepth

	case "hplb":
		// parent may have partial missing data
		if !sampleAllelesOKLenient(v, cfg.HighParentIdx, targetAlt) {
			return false
		}
		hp := v.Samples[cfg.HighParentIdx]
		lb := v.Samples[cfg.LowBulkIdx]
		if hp == nil || lb == nil {
			return false
		}
		return hp.DP >= cfg.HighParentDepth && lb.DP >= cfg.LowBulkDepth

	case "lphb":
		// parent may have partial missing data
		if !sampleAllelesOKLenient(v, cfg.LowParentIdx, targetAlt) {
			return false
		}
		lp := v.Samples[cfg.LowParentIdx]
		hb := v.Samples[cfg.HighBulkIdx]
		if lp == nil || hb == nil {
			return false
		}
		return lp.DP >= cfg.LowParentDepth && hb.DP >= cfg.HighBulkDepth

	case "lplb":
		// parent may have partial missing data
		if !sampleAllelesOKLenient(v, cfg.LowParentIdx, targetAlt) {
			return false
		}
		lp := v.Samples[cfg.LowParentIdx]
		lb := v.Samples[cfg.LowBulkIdx]
		if lp == nil || lb == nil {
			return false
		}
		return lp.DP >= cfg.LowParentDepth && lb.DP >= cfg.LowBulkDepth

	default:
		return false
	}

}

type countingWriter struct {
	io.Writer
	n int64
}

// Write wraps the embedded Writer and counts bytes written so callers can
// determine file offsets. Use atomic to be goroutine-safe.
func (cw *countingWriter) Write(p []byte) (int, error) {
	n, err := cw.Writer.Write(p)
	if n > 0 {
		atomic.AddInt64(&cw.n, int64(n))
	}
	return n, err
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
	for i := 0; i < numWorkers; i++ {
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
	//skipped := int(skippedCount.Load())

	color.Blue("Hard filtering complete: %d variant records read  → %d passed, %d rejected\n",
		original, passed, original-passed,
	)

	return passedVariants, original, passed, nil
}
