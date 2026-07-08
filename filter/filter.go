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

// ─────────────────────────────────────────────────────────────────────────
// Variant classification
// ─────────────────────────────────────────────────────────────────────────

// VariantClass is the structural category of a variant record, derived from
// comparing REF length against each real (non-'.', non-'*', non-symbolic)
// ALT allele.
type VariantClass int

const (
	ClassOther VariantClass = iota
	ClassSNP
	ClassMNP
	ClassIndel
)

func (c VariantClass) String() string {
	switch c {
	case ClassSNP:
		return "SNP"
	case ClassMNP:
		return "MNP"
	case ClassIndel:
		return "INDEL"
	default:
		return "OTHER"
	}
}

func isRealAllele(alt string) bool {
	return alt != "." && alt != "*" && !(len(alt) > 0 && alt[0] == '<')
}

// classifyVariant compares every real ALT allele's length against REF:
//   - any length mismatch   -> INDEL (a mixed SNP+indel multiallelic counts
//     as INDEL, since that's the more conservative/informative bucket)
//   - all match, REF len 1  -> SNP
//   - all match, REF len >1 -> MNP
//   - no real ALT present   -> OTHER (gVCF <NON_REF>-only, spanning deletion
//     only, etc.)
//
// Previously, an MNP (equal-length ref/alt, both >1bp) produced isSNP=false,
// isIndel=false and fell into PassesHardFilter's `default: return false`,
// silently discarding otherwise-good MNP calls. That's fixed here by giving
// MNP its own first-class branch all the way through.
func classifyVariant(v *vcfgo.Variant) VariantClass {
	refLen := len(v.Ref())
	sawReal := false
	sameLenAsRef := true

	for _, alt := range v.Alt() {
		if !isRealAllele(alt) {
			continue
		}
		sawReal = true
		if len(alt) != refLen {
			sameLenAsRef = false
		}
	}

	switch {
	case !sawReal:
		return ClassOther
	case !sameLenAsRef:
		return ClassIndel
	case refLen == 1:
		return ClassSNP
	default:
		return ClassMNP
	}
}

// singleRealAlt returns the 1-based ALT index (matching VCF/vcfgo GT allele
// numbering) of the lone real allele, or -1 if there isn't exactly one.
func singleRealAlt(v *vcfgo.Variant) int {
	realIdx := -1
	count := 0
	for i, alt := range v.Alt() {
		if isRealAllele(alt) {
			realIdx = i
			count++
		}
	}
	if count != 1 {
		return -1
	}
	return realIdx + 1 // vcfgo uses 1-based allele numbering
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

// genotypeIsOnlyRefOrAlt reports whether every allele in the sample's GT is
// either ref (0) or the target ALT allele.
func genotypeIsOnlyRefOrAlt(s *vcfgo.SampleGenotype, targetAlt int) bool {
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

// ─────────────────────────────────────────────────────────────────────────
// AD (allelic depth) filtering
// ─────────────────────────────────────────────────────────────────────────

// ADFilterConfig adds a read-support layer on top of the GT-based checks.
// GT says what was *called*; AD says how much evidence that call actually
// had. This is defined standalone so it drops in without needing changes to
// utils.AnalysisConfig — wire it into that struct later if convenient.
//
// Enabled defaults to false (zero value), so passing ADFilterConfig{}
// reproduces the previous GT-only behavior exactly.
type ADFilterConfig struct {
	Enabled bool

	// MinCalledAlleleDepth is the minimum AD value required for every
	// allele that appears in a sample's GT (e.g. a het call needs read
	// support for both the ref and alt allele it names).
	MinCalledAlleleDepth int

	// MaxOtherAlleleFraction caps how much of a sample's total AD may come
	// from alleles that are neither ref nor the variant's target ALT —
	// guards against contamination/mismapping producing an accepted GT
	// that's actually sitting on noisy multi-allelic support.
	MaxOtherAlleleFraction float64
}

// sampleADSupportsCall is a no-op pass (returns true) when cfg is disabled,
// so callers can unconditionally chain it after genotypeIsOnlyRefOrAlt.
func sampleADSupportsCall(s *vcfgo.SampleGenotype, targetAlt int, cfg ADFilterConfig) bool {
	if !cfg.Enabled {
		return true
	}
	if s == nil || len(s.AD) == 0 {
		return false
	}

	total := 0
	for _, ad := range s.AD {
		if ad < 0 {
			return false // malformed AD
		}
		total += ad
	}
	if total == 0 {
		return false
	}

	other := 0
	for i, ad := range s.AD {
		if i == 0 || i == targetAlt {
			continue
		}
		other += ad
	}
	if float64(other)/float64(total) > cfg.MaxOtherAlleleFraction {
		return false
	}

	seen := map[int]bool{}
	for _, a := range s.GT {
		if a >= 0 {
			seen[a] = true
		}
	}
	for allele := range seen {
		if allele >= len(s.AD) || s.AD[allele] < cfg.MinCalledAlleleDepth {
			return false
		}
	}
	return true
}

// sampleIsClean is the single gate every sample must pass in BsaSeqFilter:
// GT purity first (cheap), then AD support (if enabled).
func sampleIsClean(v *vcfgo.Variant, idx, targetAlt int, adCfg ADFilterConfig) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]
	return genotypeIsOnlyRefOrAlt(s, targetAlt) && sampleADSupportsCall(s, targetAlt, adCfg)
}

// ─────────────────────────────────────────────────────────────────────────
// Hard filter
// ─────────────────────────────────────────────────────────────────────────

// PassesHardFilter applies GATK-style hard-filter thresholds and also
// returns the variant's class, so callers (and the pipeline stats below)
// don't need to re-derive it.
func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) (bool, VariantClass) {
	class := classifyVariant(v)

	switch class {
	case ClassSNP, ClassMNP:
		// GATK publishes no dedicated MNP hard-filter recommendation; MNPs
		// are scored with the SNP thresholds since they come from the same
		// substitution-based annotation math (QD, FS, MQ, ...) as SNPs.
		// Split this out into its own case later if you want MNP-specific
		// thresholds.
		if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
			return false, class
		}
		if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.SNP_QD_Min {
			return false, class
		}
		if fs, ok := utils.GetFloat(v, "FS"); ok && fs > hfcfg.SNP_FS_Max {
			return false, class
		}
		if sor, ok := utils.GetFloat(v, "SOR"); ok && sor > hfcfg.SNP_SOR_Max {
			return false, class
		}
		if mq, ok := utils.GetFloat(v, "MQ"); ok && mq < hfcfg.SNP_MQ_Min {
			return false, class
		}
		if mqrs, ok := utils.GetFloat(v, "MQRankSum"); ok && mqrs < hfcfg.SNP_MQRankSum_Min {
			return false, class
		}
		if rprs, ok := utils.GetFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.SNP_ReadPosRankSum_Min {
			return false, class
		}
		return true, class

	case ClassIndel:
		if float64(v.Quality) < hfcfg.INDEL_QUAL_Min {
			return false, class
		}
		if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.INDEL_QD_Min {
			return false, class
		}
		if fs, ok := utils.GetFloat(v, "FS"); ok && fs > hfcfg.INDEL_FS_Max {
			return false, class
		}
		if sor, ok := utils.GetFloat(v, "SOR"); ok && sor > hfcfg.INDEL_SOR_Max {
			return false, class
		}
		if rprs, ok := utils.GetFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.INDEL_ReadPosRankSum_Min {
			return false, class
		}
		return true, class

	default:
		return false, class
	}
}

// ─────────────────────────────────────────────────────────────────────────
// BSA-seq genotype-pattern filter
// ─────────────────────────────────────────────────────────────────────────

func BsaSeqFilter(v *vcfgo.Variant, cfg utils.AnalysisConfig, adCfg ADFilterConfig, bsaType string) bool {
	targetAlt := singleRealAlt(v)
	if targetAlt < 0 {
		return false
	}
	clean := func(idx int) bool { return sampleIsClean(v, idx, targetAlt, adCfg) }

	switch bsaType {
	case "2b":
		for _, idx := range []int{cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		return v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth &&
			v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "2p2b":
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		hp, lp := v.Samples[cfg.HighParentIdx], v.Samples[cfg.LowParentIdx]
		if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) || hp.GT[0] == lp.GT[0] {
			return false
		}
		return hp.DP >= cfg.HighParentDepth &&
			lp.DP >= cfg.LowParentDepth &&
			v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth &&
			v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "2plb":
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		hp, lp := v.Samples[cfg.HighParentIdx], v.Samples[cfg.LowParentIdx]
		if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) || hp.GT[0] == lp.GT[0] {
			return false
		}
		return hp.DP >= cfg.HighParentDepth && lp.DP >= cfg.LowParentDepth &&
			v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "2phb":
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		hp, lp := v.Samples[cfg.HighParentIdx], v.Samples[cfg.LowParentIdx]
		if !isHomozygous(hp.GT) || !isHomozygous(lp.GT) || hp.GT[0] == lp.GT[0] {
			return false
		}
		return hp.DP >= cfg.HighParentDepth && lp.DP >= cfg.LowParentDepth &&
			v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth

	case "hp2b":
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth &&
			v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth &&
			v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "lp2b":
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth &&
			v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth &&
			v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "hphb":
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth

	case "hplb":
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	case "lphb":
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth

	case "lplb":
		for _, idx := range []int{cfg.LowParentIdx, cfg.LowBulkIdx} {
			if !clean(idx) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

	default:
		return false
	}
}

// ─────────────────────────────────────────────────────────────────────────
// Output plumbing (bgzf / tabix) — unchanged from the original
// ─────────────────────────────────────────────────────────────────────────

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
		if idx >= 0 && idx < len(v.Samples) {
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

// ─────────────────────────────────────────────────────────────────────────
// Per-class stats
// ─────────────────────────────────────────────────────────────────────────

// classStats tracks pass/fail counts per VariantClass, safe for concurrent
// use from filter workers. Having this broken out makes it immediately
// visible if, e.g., MNPs are being rejected at a suspicious rate — instead
// of everything disappearing into a single opaque "rejected" number.
type classStats struct {
	passed [4]atomic.Int64
	failed [4]atomic.Int64
}

func (cs *classStats) record(class VariantClass, passed bool) {
	if passed {
		cs.passed[class].Add(1)
	} else {
		cs.failed[class].Add(1)
	}
}

func (cs *classStats) report() string {
	var b strings.Builder
	for _, c := range []VariantClass{ClassSNP, ClassMNP, ClassIndel, ClassOther} {
		fmt.Fprintf(&b, "  %-6s pass=%-8d fail=%d\n", c, cs.passed[c].Load(), cs.failed[c].Load())
	}
	return b.String()
}

// ─────────────────────────────────────────────────────────────────────────
// Main pipeline
// ─────────────────────────────────────────────────────────────────────────

// HardFilterVcf streams the input VCF through hard-filter + BSA-seq
// genotype-pattern filtering, writing passing variants (with index) and
// rejects to separate bgzf VCFs.
//
// Filtering is parallelized across a worker pool (CPU-bound, no shared
// state), but results are reassembled into original read order before
// being written — required both for VCF's implicit sort-order convention
// and, critically, for tabix indexing, which assumes monotonically
// increasing coordinates. (The previous version wrote results in whatever
// order workers happened to finish, silently producing an unsorted VCF and
// a broken .tbi index.)
func HardFilterVcf(cfg utils.AnalysisConfig, hfcfg utils.HardFilterConfig, adCfg ADFilterConfig, bsaseqType string, keepIndices []int) ([]*vcfgo.Variant, int, int, error) {

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

	// ── Open output files ─────────────────────────────────────────────────
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

	// ── Pipeline types ────────────────────────────────────────────────────
	type job struct {
		seq int
		v   *vcfgo.Variant
	}
	type result struct {
		seq    int
		v      *vcfgo.Variant
		class  VariantClass
		passed bool
	}

	const chanBuf = 512
	jobCh := make(chan job, chanBuf)
	resultCh := make(chan result, chanBuf)
	var readerErr atomic.Pointer[error]

	var (
		originalCount atomic.Int64
		skippedCount  atomic.Int64
	)
	stats := &classStats{}

	// ── Stage 1: Reader goroutine ─────────────────────────────────────────
	go func() {
		defer close(jobCh)
		seq := 0
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
			jobCh <- job{seq: seq, v: v}
			seq++
		}
		if err := rdr.Error(); err != nil && !strings.Contains(err.Error(), "bad sample string") {
			e := fmt.Errorf("VCF read error: %w", err)
			readerErr.Store(&e)
		}
	}()

	// ── Stage 2: Filter worker pool ────────────────────────────────────────
	// CPU-bound, no shared state, one worker per core. Each variant's filter
	// call is wrapped in a recover() so a bug on one malformed record (nil
	// AD, unexpected allele count, etc.) can't panic the whole process and
	// lose every variant already in flight — it's simply routed to rejects.
	numWorkers := runtime.GOMAXPROCS(0)
	var wg sync.WaitGroup
	wg.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go func() {
			defer wg.Done()
			for j := range jobCh {
				res := result{seq: j.seq, v: j.v}
				func() {
					defer func() {
						if r := recover(); r != nil {
							res.passed = false
							res.class = ClassOther
							fmt.Printf("warning: recovered panic filtering variant at %s:%d: %v\n",
								j.v.Chromosome, j.v.Pos, r)
						}
					}()
					hfPass, class := PassesHardFilter(j.v, hfcfg)
					res.class = class
					res.passed = hfPass && BsaSeqFilter(j.v, cfg, adCfg, bsaseqType)
				}()
				stats.record(res.class, res.passed)
				resultCh <- res
			}
		}()
	}
	go func() {
		wg.Wait()
		close(resultCh)
	}()

	// ── Stage 3: Writer goroutine ───────────────────────────────────────────
	// Workers finish out of order; results are buffered in `pending` and
	// only written once they arrive in original read order (tracked by
	// `nextSeq`). Memory use is bounded by how far out-of-order workers get,
	// which in practice stays small and roughly proportional to chanBuf
	// since producers block once resultCh is full.
	var (
		passedVariants []*vcfgo.Variant
		passed         int
		writerErr      error
	)
	pending := make(map[int]result)
	nextSeq := 0

	flush := func(res result) bool { // false => stop processing further results
		sanitizeVariant(res.v, keepIndices)

		if !res.passed {
			badWriter.WriteVariant(res.v)
			return true
		}

		blockOffset, _ := hfBgzf.Next()
		startOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffset)}

		hfWriter.WriteVariant(res.v)

		blockOffsetEnd, _ := hfBgzf.Next()
		endOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffsetEnd)}

		if err := addTabixRecord(hfIdx, res.v, bgzf.Chunk{Begin: startOffset, End: endOffset}); err != nil {
			writerErr = fmt.Errorf("tabix add variant at %s:%d: %w", res.v.Chromosome, res.v.Pos, err)
			return false
		}
		passedVariants = append(passedVariants, res.v)
		passed++
		return true
	}

resultLoop:
	for res := range resultCh {
		pending[res.seq] = res
		for r, ok := pending[nextSeq]; ok; r, ok = pending[nextSeq] {
			delete(pending, nextSeq)
			nextSeq++
			if !flush(r) {
				break resultLoop
			}
		}
	}

	// Drain resultCh if the writer broke early, so blocked workers can
	// unblock and exit.
	for range resultCh {
	}

	// ── Check errors ─────────────────────────────────────────────────────
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

	// ── Flush and close bgzf streams ──────────────────────────────────────
	if err = hfBgzf.Close(); err != nil {
		badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("close hard-filtered bgzf: %w", err)
	}
	if err = badBgzf.Close(); err != nil {
		return nil, 0, 0, fmt.Errorf("close rejected bgzf: %w", err)
	}

	// ── Write tabix index ────────────────────────────────────────────────
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
	fmt.Print("By variant class:\n" + stats.report())

	return passedVariants, original, passed, nil
}
