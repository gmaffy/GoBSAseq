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
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
)

// realAltIndices returns the one-based indices of the non-symbolic ALT alleles.
func realAltIndices(v *vcfgo.Variant) []int {
	var idxs []int
	for i, alt := range v.Alt() {
		if !(alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<')) {
			idxs = append(idxs, i+1)
		}
	}
	return idxs
}

func safeSampleRefDepth(s *vcfgo.SampleGenotype) (int, error) {
	if s == nil {
		return 0, fmt.Errorf("nil sample")
	}
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

func safeSampleAltDepths(s *vcfgo.SampleGenotype) ([]int, error) {
	if s == nil {
		return []int{}, fmt.Errorf("nil sample")
	}
	if ad, ok := s.Fields["AD"]; ok {
		if comma := strings.Index(ad, ","); comma >= 0 {
			parts := strings.Split(ad[comma+1:], ",")
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

// effectiveDP is FORMAT/DP when present, otherwise the summed AD/RO+AO depths.
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

const parentAllelePurity = 0.85

// secondAltMaxFraction caps the share of a sample's reads that may support a
// second (non-target, non-ref) ALT. A small tolerance keeps clean sites over a
// stray third-allele read without touching the target-ALT frequency, so it does
// not bias the SNP-index.
const secondAltMaxFraction = 0.05

func hasUsableAlleleSignal(v *vcfgo.Variant, idx, targetAlt int) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]
	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR != nil || errA != nil || targetAlt-1 < 0 || targetAlt-1 >= len(altDepths) {
		// AD unusable: fall back to the genotype, accepting only REF/target-ALT.
		if s == nil || len(s.GT) == 0 {
			return false
		}
		for _, allele := range s.GT {
			if allele < 0 || (allele != 0 && allele != targetAlt) {
				return false
			}
		}
		return true
	}

	altDepth := altDepths[targetAlt-1]
	other := 0
	for i, d := range altDepths {
		if i != targetAlt-1 {
			other += d
		}
	}

	total := refDepth + altDepth + other
	if s.DP > total {
		total = s.DP
	}
	if total <= 0 {
		return false
	}
	// Biallelic REF-vs-target-ALT; tolerate only a small fraction of
	// contaminating second-ALT reads (REF/target are renormalised downstream).
	return float64(other) <= secondAltMaxFraction*float64(total)
}

// sampleGQAtLeast reports whether FORMAT/GQ meets minGQ; a missing or
// unparseable GQ is not penalised.
func sampleGQAtLeast(s *vcfgo.SampleGenotype, minGQ int) bool {
	if s == nil {
		return false
	}
	if raw, ok := s.Fields["GQ"]; ok {
		if raw == "" || raw == "." {
			return true
		}
		if g, err := strconv.Atoi(raw); err == nil {
			return g >= minGQ
		}
	}
	if s.GQ > 0 {
		return s.GQ >= minGQ
	}
	return true
}

// sampleUsable checks that depth is within [minDepth, maxDepth] (maxDepth <= 0
// disables the upper bound) and, when minGQ > 0 and GQ is present, that GQ meets
// the floor. The maximum-depth cap drops coverage outliers (collapsed repeats,
// paralogues, CNVs) whose allele ratios are mapping artifacts, not segregation.
func sampleUsable(v *vcfgo.Variant, idx, minDepth, maxDepth, minGQ int) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]
	dp := effectiveDP(s)
	if dp < minDepth {
		return false
	}
	if maxDepth > 0 && dp > maxDepth {
		return false
	}
	if minGQ > 0 && !sampleGQAtLeast(s, minGQ) {
		return false
	}
	return true
}

// parentAllele returns the allele a parent predominantly carries (0 = ref,
// targetAlt = alt) and whether that call is confident. It prefers allele depth
// over GT so polarity is robust to imperfectly-homozygous parents.
func parentAllele(v *vcfgo.Variant, idx, targetAlt int) (allele int, confident bool) {
	if idx < 0 || idx >= len(v.Samples) {
		return 0, false
	}
	s := v.Samples[idx]

	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR == nil && errA == nil && targetAlt-1 >= 0 && targetAlt-1 < len(altDepths) {
		total := refDepth
		for _, d := range altDepths {
			total += d
		}
		if total > 0 {
			refFrac := float64(refDepth) / float64(total)
			nonRef := total - refDepth
			switch {
			case refFrac >= parentAllelePurity || (nonRef <= 1 && refDepth > nonRef):
				return 0, true
			case refFrac <= 1-parentAllelePurity || (refDepth <= 1 && nonRef > refDepth):
				return targetAlt, true
			default:
				return 0, false
			}
		}
	}

	// Depth unusable: fall back to a homozygous genotype call.
	if len(s.GT) > 0 && s.GT[0] >= 0 {
		hom := true
		for _, a := range s.GT[1:] {
			if a < 0 || a != s.GT[0] {
				hom = false
				break
			}
		}
		if hom {
			return s.GT[0], true
		}
	}
	return 0, false
}

// classifyVariant reports whether a record is a SNP/MNP or an indel. Same-length
// substitutions (SNPs and multi-nucleotide substitutions) are treated as SNPs;
// any length-changing ALT makes it an indel.
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
	isSNP = sameLength
	return
}

// FilterProfile records which quality annotations the caller wrote into the VCF,
// derived once from the header. It lets the pipeline report the active filter set
// honestly and lean on GQ when GATK annotations are absent (e.g. DeepVariant).
type FilterProfile struct {
	Caller            string
	HasQD             bool
	HasFS             bool
	HasSOR            bool
	HasMQ             bool
	HasMQRankSum      bool
	HasReadPosRankSum bool
	HasGQ             bool
}

// HasGATKAnnotations reports whether the core GATK hard-filter annotations exist.
func (p FilterProfile) HasGATKAnnotations() bool {
	return p.HasQD || p.HasFS || p.HasSOR || p.HasMQ
}

// DetectFilterProfile inspects the header for available annotations and makes a
// best-effort guess at the calling tool (for reporting).
func DetectFilterProfile(hdr *vcfgo.Header) FilterProfile {
	p := FilterProfile{Caller: "unknown"}
	if hdr == nil {
		return p
	}
	hasInfo := func(id string) bool { _, ok := hdr.Infos[id]; return ok }
	p.HasQD = hasInfo("QD")
	p.HasFS = hasInfo("FS")
	p.HasSOR = hasInfo("SOR")
	p.HasMQ = hasInfo("MQ")
	p.HasMQRankSum = hasInfo("MQRankSum")
	p.HasReadPosRankSum = hasInfo("ReadPosRankSum")
	if _, ok := hdr.SampleFormats["GQ"]; ok {
		p.HasGQ = true
	}

	meta := strings.ToLower(strings.Join(hdr.Extras, " "))
	switch {
	case strings.Contains(meta, "deepvariant"):
		p.Caller = "deepvariant"
	case strings.Contains(meta, "haplotypecaller") || strings.Contains(meta, "gatk"):
		p.Caller = "gatk"
	case p.HasGATKAnnotations():
		p.Caller = "gatk-like"
	case p.HasGQ:
		p.Caller = "deepvariant-like"
	}
	return p
}

// reportActiveFilters prints which hard filters will run given the annotations
// present, plus the BSA-specific gates in effect.
func reportActiveFilters(p FilterProfile, hfcfg utils.HardFilterConfig, cfg utils.AnalysisConfig) {
	mode := "GATK best-practice (strict)"
	if hfcfg.LightFilter {
		mode = "light (BSA-tuned: strand kept, rank-sums dropped)"
	}
	color.Blue("Detected caller: %s | filter mode: %s", p.Caller, mode)

	active := []string{"FILTER=PASS", "QUAL"}
	add := func(present bool, name string) {
		if present {
			active = append(active, name)
		}
	}
	add(p.HasQD, "QD")
	add(p.HasMQ, "MQ")
	add(p.HasFS, "FS")
	add(p.HasSOR, "SOR")
	if !hfcfg.LightFilter {
		add(p.HasMQRankSum, "MQRankSum")
		add(p.HasReadPosRankSum, "ReadPosRankSum")
	}
	color.Blue("Active annotation filters: %s", strings.Join(active, ", "))

	if !p.HasGATKAnnotations() {
		color.Yellow("No GATK annotations found (e.g. DeepVariant): relying on FILTER=PASS, QUAL, depth%s.",
			map[bool]string{true: ", GQ (parents)", false: ""}[p.HasGQ && cfg.MinGQ > 0])
	}
	if cfg.MinGQ > 0 {
		if p.HasGQ {
			color.Blue("BSA gate: min GQ = %d (parent samples only; bulks are depth-only to avoid AF bias)", cfg.MinGQ)
		} else {
			color.Yellow("min-GQ=%d requested but GQ is not present in this VCF; skipping GQ filter.", cfg.MinGQ)
		}
	}
	if cfg.SplitMultiallelic {
		color.Blue("Multi-allelic decomposition: enabled (records split to biallelic before filtering)")
	} else {
		color.Yellow("Multi-allelic decomposition: disabled")
	}
	color.Blue("BSA gate: biallelic second-ALT tolerance = %.0f%% of depth", secondAltMaxFraction*100)
	if cfg.HighBulkMaxDepth > 0 || cfg.LowBulkMaxDepth > 0 || cfg.HighParentMaxDepth > 0 || cfg.LowParentMaxDepth > 0 {
		color.Blue("BSA gate: max-depth caps — HB:%d LB:%d HP:%d LP:%d (0 = off)",
			cfg.HighBulkMaxDepth, cfg.LowBulkMaxDepth, cfg.HighParentMaxDepth, cfg.LowParentMaxDepth)
	}
}

// PassesHardFilter applies caller-agnostic hard filtering. FILTER != PASS and a
// QUAL floor are always enforced; each annotation threshold applies only when the
// annotation is present, so the same path is correct for GATK and DeepVariant.
// FS/SOR stay valid in a pool and run in both modes; MQRankSum/ReadPosRankSum are
// unreliable near allele fixation (at QTLs), so light mode drops them.
func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
	if v.Filter != "" && v.Filter != "." && v.Filter != "PASS" {
		return false
	}
	isSNP, isIndel := classifyVariant(v)
	strict := !hfcfg.LightFilter

	// atLeast/atMost pass when the annotation is absent, else enforce the bound.
	atLeast := func(key string, min float64) bool {
		val, ok := utils.GetFloat(v, key)
		return !ok || val >= min
	}
	atMost := func(key string, max float64) bool {
		val, ok := utils.GetFloat(v, key)
		return !ok || val <= max
	}

	switch {
	case isSNP:
		if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
			return false
		}
		if !atLeast("QD", hfcfg.SNP_QD_Min) ||
			!atLeast("MQ", hfcfg.SNP_MQ_Min) ||
			!atMost("FS", hfcfg.SNP_FS_Max) ||
			!atMost("SOR", hfcfg.SNP_SOR_Max) {
			return false
		}
		if strict {
			if !atLeast("MQRankSum", hfcfg.SNP_MQRankSum_Min) ||
				!atLeast("ReadPosRankSum", hfcfg.SNP_ReadPosRankSum_Min) {
				return false
			}
		}
		return true

	case isIndel:
		if float64(v.Quality) < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if !atLeast("QD", hfcfg.INDEL_QD_Min) ||
			!atMost("FS", hfcfg.INDEL_FS_Max) ||
			!atMost("SOR", hfcfg.INDEL_SOR_Max) {
			return false
		}
		if strict {
			if !atLeast("ReadPosRankSum", hfcfg.INDEL_ReadPosRankSum_Min) {
				return false
			}
		}
		return true

	default:
		return false
	}
}

// BsaSeqTargetAlt returns the one-based ALT index that satisfies the BSA-seq
// segregation filter, or zero when none does. Trying each real ALT keeps
// multi-allelic sites that pass through an ALT other than the first.
func BsaSeqTargetAlt(v *vcfgo.Variant, cfg utils.AnalysisConfig, bsaType string) int {
	for _, targetAlt := range realAltIndices(v) {
		if bsaSeqFilterAllele(v, cfg, bsaType, targetAlt) {
			return targetAlt
		}
	}
	return 0
}

func bsaSeqFilterAllele(v *vcfgo.Variant, cfg utils.AnalysisConfig, bsaType string, targetAlt int) bool {
	minGQ := cfg.MinGQ

	// Bulks: usable biallelic signal + in-range depth, but NO GQ floor. GQ in a
	// pool is AF-correlated, so gating bulks on it would bias the SNP-index.
	bulkOK := func(idx, minDepth, maxDepth int) bool {
		return hasUsableAlleleSignal(v, idx, targetAlt) && sampleUsable(v, idx, minDepth, maxDepth, 0)
	}
	// Lone parent in single-parent modes: allele signal + depth + GQ (a parent is
	// a real genotype, so GQ is meaningful).
	parentSigOK := func(idx, minDepth, maxDepth int) bool {
		return hasUsableAlleleSignal(v, idx, targetAlt) && sampleUsable(v, idx, minDepth, maxDepth, minGQ)
	}
	// Parents in two-parent modes: depth + GQ; informativeness comes from parentAllele.
	parentDepthOK := func(idx, minDepth, maxDepth int) bool {
		return sampleUsable(v, idx, minDepth, maxDepth, minGQ)
	}
	// Both parents must call confidently and disagree, else the site is uninformative.
	parentsInformative := func() bool {
		hpAllele, hpOK := parentAllele(v, cfg.HighParentIdx, targetAlt)
		lpAllele, lpOK := parentAllele(v, cfg.LowParentIdx, targetAlt)
		return hpOK && lpOK && hpAllele != lpAllele
	}

	switch bsaType {
	case "2b":
		return bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth) &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	case "2p2b":
		return parentDepthOK(cfg.HighParentIdx, cfg.HighParentDepth, cfg.HighParentMaxDepth) &&
			parentDepthOK(cfg.LowParentIdx, cfg.LowParentDepth, cfg.LowParentMaxDepth) &&
			parentsInformative() &&
			bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth) &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	case "2plb":
		return parentDepthOK(cfg.HighParentIdx, cfg.HighParentDepth, cfg.HighParentMaxDepth) &&
			parentDepthOK(cfg.LowParentIdx, cfg.LowParentDepth, cfg.LowParentMaxDepth) &&
			parentsInformative() &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	case "2phb":
		return parentDepthOK(cfg.HighParentIdx, cfg.HighParentDepth, cfg.HighParentMaxDepth) &&
			parentDepthOK(cfg.LowParentIdx, cfg.LowParentDepth, cfg.LowParentMaxDepth) &&
			parentsInformative() &&
			bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth)

	case "hp2b":
		return parentSigOK(cfg.HighParentIdx, cfg.HighParentDepth, cfg.HighParentMaxDepth) &&
			bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth) &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	case "lp2b":
		return parentSigOK(cfg.LowParentIdx, cfg.LowParentDepth, cfg.LowParentMaxDepth) &&
			bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth) &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	case "hphb":
		return parentSigOK(cfg.HighParentIdx, cfg.HighParentDepth, cfg.HighParentMaxDepth) &&
			bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth)

	case "hplb":
		return parentSigOK(cfg.HighParentIdx, cfg.HighParentDepth, cfg.HighParentMaxDepth) &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	case "lphb":
		return parentSigOK(cfg.LowParentIdx, cfg.LowParentDepth, cfg.LowParentMaxDepth) &&
			bulkOK(cfg.HighBulkIdx, cfg.HighBulkDepth, cfg.HighBulkMaxDepth)

	case "lplb":
		return parentSigOK(cfg.LowParentIdx, cfg.LowParentDepth, cfg.LowParentMaxDepth) &&
			bulkOK(cfg.LowBulkIdx, cfg.LowBulkDepth, cfg.LowBulkMaxDepth)

	default:
		return false
	}
}

// countingWriter tracks bytes written so tabix chunk offsets can be recorded.
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

// vcfRecord implements the tabix record interface for index building.
type vcfRecord struct {
	chrom string
	start int // 0-based
	end   int // 0-based half-open
}

func (r vcfRecord) RefName() string { return r.chrom }
func (r vcfRecord) Start() int      { return r.start }
func (r vcfRecord) End() int        { return r.end }

// FilteredVariant associates a passing record with the ALT allele that satisfied
// the BSA-seq filter. TargetAlt is one-based VCF allele numbering.
type FilteredVariant struct {
	Variant   *vcfgo.Variant
	TargetAlt int
}

func HardFilterVcf(cfg utils.AnalysisConfig, hfcfg utils.HardFilterConfig, bsaseqType string, keepIndices []int) ([]FilteredVariant, int, int, error) {

	if err := os.MkdirAll(filepath.Join(cfg.OutputDir, "stats"), 0775); err != nil {
		return nil, 0, 0, err
	}

	hardFilteredVcfPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.hardfiltered.vcf.gz", bsaseqType))
	badVcfPath := filepath.Join(cfg.OutputDir, "stats", fmt.Sprintf("GoBSAseq.%s.lowqual.vcf.gz", bsaseqType))

	rdr := cfg.Rdr
	origSampleNames := rdr.Header.SampleNames

	profile := DetectFilterProfile(rdr.Header)
	reportActiveFilters(profile, hfcfg, cfg)

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

	// Output files.
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

	hfIdx := tabix.New()
	hfIdx.Format = 2 // VCF
	hfIdx.NameColumn = 1
	hfIdx.BeginColumn = 2
	hfIdx.EndColumn = 0
	hfIdx.MetaChar = '#'

	bar := progressbar.Default(-1, "Hard filtering variants")

	type variantResult struct {
		seq       int
		v         *vcfgo.Variant
		targetAlt int
		passed    bool
	}
	type queuedVariant struct {
		seq int
		v   *vcfgo.Variant
	}

	const chanBuf = 512
	filterCh := make(chan queuedVariant, chanBuf) // reader  → filter workers
	resultCh := make(chan variantResult, chanBuf) // workers → writer
	var readerErr atomic.Pointer[error]
	var originalCount atomic.Int64

	// Stage 1: reader. Preserves source order via seq so the writer can emit a
	// coordinate-sorted VCF despite concurrent filtering.
	go func() {
		defer close(filterCh)
		seq := 0
		for {
			v := rdr.Read()
			if v == nil {
				break
			}
			if err := rdr.Error(); err != nil {
				// vcfgo mis-handles "." (a valid missing-value marker) in some
				// numeric INFO/FORMAT fields; that is non-fatal, so clear and go on.
				if strings.Contains(err.Error(), "bad sample string") ||
					strings.Contains(err.Error(), "strconv.ParseFloat") {
					rdr.Clear()
				} else {
					e := fmt.Errorf("VCF parse error at line %d: %w", v.LineNumber, err)
					readerErr.Store(&e)
					return
				}
			}
			alts := v.Alt()
			if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
				continue
			}

			// Decompose multi-allelics so every record entering the filter is
			// biallelic; already-biallelic records pass through unchanged. Split
			// records share a POS, so the writer's order buffer stays sorted.
			records := []*vcfgo.Variant{v}
			if cfg.SplitMultiallelic {
				records = SplitMultiallelic(v)
			}
			for _, rec := range records {
				originalCount.Add(1)
				_ = bar.Add(1)
				filterCh <- queuedVariant{seq: seq, v: rec}
				seq++
			}
		}

		if err := rdr.Error(); err != nil &&
			!strings.Contains(err.Error(), "bad sample string") &&
			!strings.Contains(err.Error(), "strconv.ParseFloat") {
			e := fmt.Errorf("VCF read error: %w", err)
			readerErr.Store(&e)
		}
	}()

	// Stage 2: filter workers. Filtering is CPU-bound and reads only immutable
	// variant fields, so no locking is needed; run one worker per core.
	numWorkers := runtime.GOMAXPROCS(0)
	var wg sync.WaitGroup
	wg.Add(numWorkers)
	for i := 0; i < numWorkers; i++ {
		go func() {
			defer wg.Done()
			for queued := range filterCh {
				v := queued.v
				targetAlt := 0
				passed := PassesHardFilter(v, hfcfg)
				if passed {
					targetAlt = BsaSeqTargetAlt(v, cfg, bsaseqType)
					passed = targetAlt != 0
				}
				resultCh <- variantResult{seq: queued.seq, v: v, targetAlt: targetAlt, passed: passed}
			}
		}()
	}
	go func() {
		wg.Wait()
		close(resultCh)
	}()

	// Stage 3: single writer (bgzf / vcfgo writers are not goroutine-safe).
	var (
		passedVariants []FilteredVariant
		passed         int
		writerErr      error
	)

	writeResult := func(res variantResult) error {
		// Subset samples to kept indices and drop the PGT/PID format fields.
		newSamples := make([]*vcfgo.SampleGenotype, len(keepIndices))
		for i, idx := range keepIndices {
			if idx < len(res.v.Samples) {
				newSamples[i] = res.v.Samples[idx]
			}
		}
		res.v.Samples = newSamples
		newFormat := make([]string, 0, len(res.v.Format))
		for _, f := range res.v.Format {
			if f != "PGT" && f != "PID" {
				newFormat = append(newFormat, f)
			}
		}
		res.v.Format = newFormat

		if !res.passed {
			badWriter.WriteVariant(res.v)
			return nil
		}

		blockOffset, _ := hfBgzf.Next()
		startOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffset)}
		hfWriter.WriteVariant(res.v)
		blockOffsetEnd, _ := hfBgzf.Next()
		endOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffsetEnd)}

		rec := vcfRecord{chrom: res.v.Chromosome, start: int(res.v.Pos) - 1, end: int(res.v.Pos) - 1 + len(res.v.Ref())}
		if err := hfIdx.Add(rec, bgzf.Chunk{Begin: startOffset, End: endOffset}, true, true); err != nil {
			return fmt.Errorf("tabix add variant at %s:%d: %w", res.v.Chromosome, res.v.Pos, err)
		}

		passedVariants = append(passedVariants, FilteredVariant{Variant: res.v, TargetAlt: res.targetAlt})
		passed++
		return nil
	}

	// Workers complete out of order; buffer results until the next source record
	// arrives so VCF output and tabix chunks stay coordinate-sorted.
	pending := make(map[int]variantResult)
	nextSeq := 0
	for res := range resultCh {
		pending[res.seq] = res
		for {
			next, ok := pending[nextSeq]
			if !ok {
				break
			}
			delete(pending, nextSeq)
			if err := writeResult(next); err != nil {
				writerErr = err
				break
			}
			nextSeq++
		}
		if writerErr != nil {
			break
		}
	}
	// Drain if the writer broke early so workers can unblock and exit.
	for range resultCh {
	}

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

	if err = hfBgzf.Close(); err != nil {
		badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("close hard-filtered bgzf: %w", err)
	}
	if err = badBgzf.Close(); err != nil {
		return nil, 0, 0, fmt.Errorf("close rejected bgzf: %w", err)
	}

	// Tabix index.
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
	color.Blue(
		"Hard filtering complete:\n%d variant records read.\n%d passed.\n%d rejected\n",
		original, passed, original-passed,
	)

	return passedVariants, original, passed, nil
}
