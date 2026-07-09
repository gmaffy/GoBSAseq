package filter

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/tabix"
	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
)

// realAltIndices returns the 1-based allele indices of every "real" ALT at
// a site (excluding spanning-deletion '*', missing '.', and symbolic <...>
// alleles). Previously only sites with exactly one real ALT were kept and
// everything else was discarded outright; but a BSA-seq cross typically only
// segregates for one of several ALTs at a genuinely multi-allelic site, so
// callers should try each candidate rather than rejecting the whole record.
func realAltIndices(v *vcfgo.Variant) []int {
	var idxs []int
	for i, alt := range v.Alt() {
		if !(alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<')) {
			idxs = append(idxs, i+1) // vcfgo uses 1-based allele numbering
		}
	}
	return idxs
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

// safeSampleRefDepth parses REF depth for a sample without calling vcfgo's
// RefDepth() directly. vcfgo can panic when AD exists but is malformed
// (e.g. missing the comma separator). On parse failure it returns an error.
func safeSampleRefDepth(s *vcfgo.SampleGenotype) (int, error) {
	if s == nil {
		return 0, fmt.Errorf("nil sample")
	}

	// Prefer AD when it is well-formed; if AD exists but lacks the comma,
	// vcfgo would panic, so fall back to RO when available.
	if ad, ok := s.Fields["AD"]; ok {
		if comma := strings.Index(ad, ","); comma >= 0 {
			refStr := ad[:comma]
			if refStr == "" || refStr == "." {
				return 0, fmt.Errorf("invalid AD ref depth %q", refStr)
			}
			return strconv.Atoi(refStr)
		}
	}

	if ro, ok := s.Fields["RO"]; ok {
		if ro == "" || ro == "." {
			return 0, fmt.Errorf("invalid RO ref depth %q", ro)
		}
		return strconv.Atoi(ro)
	}

	return 0, fmt.Errorf("no ref depth field (AD/RO)")
}

// safeSampleAltDepths parses ALT depths for a sample without calling
// vcfgo's AltDepths() directly. It avoids panics/incorrect parsing when AD is
// malformed (e.g. missing the comma separator). On parse failure it returns
// an error.
func safeSampleAltDepths(s *vcfgo.SampleGenotype) ([]int, error) {
	if s == nil {
		return []int{}, fmt.Errorf("nil sample")
	}

	if ad, ok := s.Fields["AD"]; ok {
		if comma := strings.Index(ad, ","); comma >= 0 {
			altStr := ad[comma+1:]
			parts := strings.Split(altStr, ",")
			vals := make([]int, len(parts))
			for i, p := range parts {
				if p == "" || p == "." {
					return []int{}, fmt.Errorf("invalid AD alt depth %q", p)
				}
				v, err := strconv.Atoi(p)
				if err != nil {
					return []int{}, err
				}
				vals[i] = v
			}
			return vals, nil
		}
		// If AD is present but malformed (no comma), fall through to AO.
	}

	if ao, ok := s.Fields["AO"]; ok {
		parts := strings.Split(ao, ",")
		vals := make([]int, len(parts))
		for i, p := range parts {
			if p == "" || p == "." {
				return []int{}, fmt.Errorf("invalid AO alt depth %q", p)
			}
			v, err := strconv.Atoi(p)
			if err != nil {
				return []int{}, err
			}
			vals[i] = v
		}
		return vals, nil
	}

	return []int{}, fmt.Errorf("no alt depth field (AD/AO)")
}

// effectiveDP returns a sample's total depth for coverage-threshold checks,
// falling back to AD-derived depth (ref + all alts) when FORMAT DP is
// missing or zero. Without this, a sample whose caller only populates AD
// (no DP) could pass the AD-based allele-signal check yet still be rejected
// by a depth-threshold comparison that only looks at raw DP.
func effectiveDP(s *vcfgo.SampleGenotype) int {
	if s == nil {
		return 0
	}
	if s.DP > 0 {
		return s.DP
	}
	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR != nil || errA != nil {
		return 0
	}
	total := refDepth
	for _, d := range altDepths {
		total += d
	}
	return total
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

const (
	// maxOtherAlleleFrac bounds how much of a sample's total depth may come
	// from alleles other than REF/target-ALT before the site is considered
	// too messy to trust (contamination, paralogs, a genuinely
	// multi-allelic site). This is intentionally loose: BSA-seq bulks are
	// pools of many individuals, so some background from other alleles is
	// expected and, on its own, is not a reason to discard an otherwise
	// informative variant.
	maxOtherAlleleFrac = 0.20

	// parentAllelePurity is the minimum fraction of a parent's reads that
	// must support one allele before that allele is treated as "the
	// parent's allele". BSA-seq parents are ideally fully homozygous inbred
	// lines, but residual heterozygosity, low-level contamination, and
	// genotyping noise are common in practice. Requiring a strict
	// homozygous GT call throws away many genuinely informative markers;
	// tolerating up to (1 - parentAllelePurity) reads from the other allele
	// keeps those markers while still rejecting parents that are
	// genuinely segregating/heterozygous at a site (which would be
	// uninformative for BSA-seq regardless of filtering strategy).
	parentAllelePurity = 0.85
)

// hasUsableAlleleSignal reports whether a sample's allele-depth data at this
// site is trustworthy enough to use, without requiring a confident discrete
// genotype call. Pooled bulk samples are frequently emitted with a missing
// or low-confidence GT by callers designed around single-diploid-sample
// genotype likelihoods, even when read support itself is perfectly good —
// gating on GT there discards real BSA-seq signal along with noise. Instead,
// the sample is only rejected if there's no usable depth at all, or if a
// large share of its reads support some allele other than REF/target-ALT.
func hasUsableAlleleSignal(v *vcfgo.Variant, idx, targetAlt int) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]

	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR != nil || errA != nil || targetAlt-1 < 0 || targetAlt-1 >= len(altDepths) {
		// No usable AD — fall back to the strict GT check rather than
		// blindly accepting a site we can't actually interpret.
		return sampleHasOnlyRefOrAlt(v, idx, targetAlt)
	}

	altDepth := altDepths[targetAlt-1]
	other := 0
	for i, d := range altDepths {
		if i != targetAlt-1 {
			other += d
		}
	}

	// Use AD-derived total (ref + all alts) rather than raw FORMAT DP: DP
	// often includes reads that never resolved to an allele call (low-qual,
	// soft-clipped, multi-mapping), which would inflate "other" and fail
	// otherwise-clean samples. Fall back to DP only if it accounts for more
	// depth than AD does, and only for the denominator — it never
	// contributes to the "other allele" numerator.
	adTotal := refDepth + altDepth + other
	total := adTotal
	if s.DP > adTotal {
		total = s.DP
	}
	if total <= 0 {
		return false
	}
	return float64(other)/float64(total) <= maxOtherAlleleFrac
}

// parentAllele determines which allele (0 = ref, targetAlt = alt) predominates
// in a parent sample using allele read depths, falling back to the genotype
// call only when depth data isn't usable. confident is false when the
// parent's reads are too evenly split between alleles to trust — i.e. the
// parent looks genuinely heterozygous/segregating at this site, which makes
// it uninformative for BSA-seq regardless of how strict the filter is.
func parentAllele(v *vcfgo.Variant, idx, targetAlt int) (allele int, confident bool) {
	if idx < 0 || idx >= len(v.Samples) {
		return 0, false
	}
	s := v.Samples[idx]

	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR == nil && errA == nil && targetAlt-1 >= 0 && targetAlt-1 < len(altDepths) {
		// Total includes reads supporting *any* allele, not just
		// ref-vs-target-alt, so a parent with real off-target/contaminating
		// reads is judged against its full read pool — matching how
		// hasUsableAlleleSignal treats bulks.
		total := refDepth
		for _, d := range altDepths {
			total += d
		}
		if total > 0 {
			refFrac := float64(refDepth) / float64(total)
			nonRef := total - refDepth
			switch {
			// A single discordant read is tolerated regardless of the
			// fraction it implies — but only when the other allele still
			// clearly outnumbers it (guards against e.g. a lone alt read
			// with zero ref support being misread as "ref-confident").
			// At low parent depth (e.g. the default 5x), this stops one
			// sequencing error from pushing a genuinely homozygous parent
			// below the flat purity threshold.
			case refFrac >= parentAllelePurity || (nonRef <= 1 && refDepth > nonRef):
				return 0, true
			case refFrac <= 1-parentAllelePurity || (refDepth <= 1 && nonRef > refDepth):
				return targetAlt, true
			default:
				return 0, false
			}
		}
	}

	// No usable AD — fall back to the genotype call.
	if len(s.GT) > 0 && s.GT[0] >= 0 && isHomozygous(s.GT) {
		return s.GT[0], true
	}
	return 0, false
}

func classifyVariant(v *vcfgo.Variant) (isSNP, isIndel bool) {
	refLen := len(v.Ref())
	sameLength := true
	for _, alt := range v.Alt() {
		if alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<') {
			continue
		}
		if len(alt) != refLen {
			sameLength = false
			isIndel = true
		}
	}
	// A same-length substitution is a classic SNP when refLen == 1, and an
	// MNP (multi-nucleotide substitution, e.g. AT>GC) when refLen > 1. MNPs
	// carry the same annotations as SNPs and aren't indel alignment
	// artifacts, so they should be filtered with SNP thresholds rather than
	// falling into neither bucket, which previously caused every MNP to be
	// silently discarded regardless of quality.
	isSNP = sameLength
	return
}

func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
	isSNP, isIndel := classifyVariant(v)

	if hfcfg.LightFilter {
		// Light filtering: reject only variants with genuinely poor support.
		// QUAL, QD (SNP/INDEL) and MQ (SNP) reflect absolute call confidence
		// and read-mapping quality, and stay meaningful regardless of how
		// the sample was pooled. FS, SOR, MQRankSum and ReadPosRankSum, by
		// contrast, test whether ref- and alt-supporting reads look
		// different from each other (strand, position, mapping quality) —
		// in a bulk that is a pool of many individuals, ref/alt reads are
		// *expected* to be present in a skewed, non-50/50 ratio by design,
		// which these tests can easily mistake for bias/artifact. Skipping
		// them here keeps real segregating BSA-seq variants that GATK's
		// single-sample-tuned best-practice filters would otherwise discard.
		switch {
		case isSNP:
			if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
				return false
			}
			if qd, ok := utils.GetFloat(v, "QD"); ok && qd < hfcfg.SNP_QD_Min {
				return false
			}
			if mq, ok := utils.GetFloat(v, "MQ"); ok && mq < hfcfg.SNP_MQ_Min {
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
			return true
		default:
			return false
		}
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
	// Try each real ALT as a candidate target and keep the variant if any
	// one of them shows a valid BSA-seq segregation pattern, instead of
	// discarding multi-allelic sites outright.
	for _, targetAlt := range realAltIndices(v) {
		if bsaSeqFilterAllele(v, cfg, bsaType, targetAlt) {
			return true
		}
	}
	return false
}

func bsaSeqFilterAllele(v *vcfgo.Variant, cfg utils.AnalysisConfig, bsaType string, targetAlt int) bool {
	switch bsaType {
	case "2b":
		// bulks only
		for _, idx := range []int{cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}
		return effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth && effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth

	case "2p2b":
		// 2 parents 2 bulks filter
		for _, idx := range []int{cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]

		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		if !hpOK || !lpOK {
			return false // parent allele too ambiguous/mixed to be informative
		}
		if hpAllele == lpAllele {
			return false // same predominant allele in both parents — not informative
		}

		return effectiveDP(hp) >= cfg.HighParentDepth &&
			effectiveDP(lp) >= cfg.LowParentDepth &&
			effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth &&
			effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth
	case "2plb":
		// 2 parents low bulk filter
		if !hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt) {
			return false
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]

		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		if !hpOK || !lpOK {
			return false
		}
		if hpAllele == lpAllele {
			return false // same predominant allele in both parents — not informative
		}

		return effectiveDP(hp) >= cfg.HighParentDepth && effectiveDP(lp) >= cfg.LowParentDepth && effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth
	case "2phb":
		// 2 parents high bulk filter
		if !hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt) {
			return false
		}

		hp := v.Samples[cfg.HighParentIdx]
		lp := v.Samples[cfg.LowParentIdx]

		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		if !hpOK || !lpOK {
			return false
		}
		if hpAllele == lpAllele {
			return false // same predominant allele in both parents — not informative
		}

		return effectiveDP(hp) >= cfg.HighParentDepth && effectiveDP(lp) >= cfg.LowParentDepth && effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth
	case "hp2b":
		// high parent 2 bulks filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		return effectiveDP(hp) >= cfg.HighParentDepth && effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth && effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth
	case "lp2b":
		// low parent  2 bulks filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return effectiveDP(lp) >= cfg.LowParentDepth && effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth && effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth

	case "hphb":
		// high parent high bulk filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return effectiveDP(hp) >= cfg.HighParentDepth && effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth
	case "hplb":
		// high parent low bulk filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.LowBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}
		hp := v.Samples[cfg.HighParentIdx]
		return effectiveDP(hp) >= cfg.HighParentDepth && effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth
	case "lphb":
		// low parent high bulk filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return effectiveDP(lp) >= cfg.LowParentDepth && effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth
	case "lplb":
		// low parent low bulk filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.LowBulkIdx} {
			if !hasUsableAlleleSignal(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return effectiveDP(lp) >= cfg.LowParentDepth && effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth

	default:
		return false
	}

}

type countingWriter struct {
	io.Writer
	n int64
}

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
				if strings.Contains(err.Error(), "bad sample string") ||
					strings.Contains(err.Error(), "strconv.ParseFloat") {
					// vcfgo fails to treat "." (a valid VCF missing-value
					// marker) as a numeric placeholder in some INFO/FORMAT
					// fields. This is non-fatal - the field is simply left
					// unset for this variant - so clear the error and continue.
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

		if err := rdr.Error(); err != nil &&
			!strings.Contains(err.Error(), "bad sample string") &&
			!strings.Contains(err.Error(), "strconv.ParseFloat") {
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
	skipped := int(skippedCount.Load())

	fmt.Printf(
		"Hard filtering complete: %d variant records read (%d gVCF ref blocks skipped) → %d passed, %d rejected\n",
		original, skipped, passed, original-passed,
	)

	return passedVariants, original, passed, nil
}
