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

// hasAnyRealAlt is a fast pre-filter check: returns true if the variant has
// at least one real ALT allele. Used in the reader goroutine to skip
// non-informative records before they enter the filter pipeline.
func hasAnyRealAlt(v *vcfgo.Variant) bool {
	for _, alt := range v.Alt() {
		if alt != "." && alt != "*" && !(len(alt) > 0 && alt[0] == '<') {
			return true
		}
	}
	return false
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
// missing or zero. Unlike the previous version, it uses whatever partial
// depth information is available rather than returning 0 if either ref or
// alt depth parsing fails.
func effectiveDP(s *vcfgo.SampleGenotype) int {
	if s == nil {
		return 0
	}
	if s.DP > 0 {
		return s.DP
	}
	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)

	// Use whatever depth information we have — partial is better than zero.
	total := 0
	if errR == nil {
		total = refDepth
	}
	if errA == nil {
		for _, d := range altDepths {
			total += d
		}
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
	// maxOtherAlleleFracBase is the baseline bound for how much of a
	// sample's total depth may come from alleles other than REF/target-ALT.
	// The actual threshold is depth-dependent (see otherAlleleThreshold).
	maxOtherAlleleFracBase = 0.25

	// parentAllelePurityBase is the baseline minimum fraction of a parent's
	// reads that must support one allele. The actual threshold is
	// depth-dependent (see parentAlleleThreshold).
	parentAllelePurityBase = 0.85
)

// otherAlleleThreshold returns a depth-dependent threshold for the maximum
// allowed fraction of reads from "other" alleles. At low depth, sequencing
// noise can easily create spurious "other" reads, so we're tighter; at high
// depth the estimate is more precise and we can afford to be more lenient
// with genuine population-level multi-allelic variation.
func otherAlleleThreshold(totalDepth int) float64 {
	switch {
	case totalDepth <= 10:
		return 0.15
	case totalDepth <= 30:
		return 0.20
	default:
		return maxOtherAlleleFracBase
	}
}

// parentAlleleThreshold returns a depth-dependent purity threshold for
// determining whether a parent's allele call is confident. At low depths
// (common in BSA-seq parents), a single sequencing error can push a truly
// homozygous sample below a flat threshold — e.g. 4:1 at 5x is 80%, which
// would fail an 85% cutoff despite being strong evidence for homozygosity.
func parentAlleleThreshold(totalDepth int) float64 {
	switch {
	case totalDepth <= 5:
		return 0.65 // 3:1 at 4x (75%) passes; 4:1 at 5x (80%) passes
	case totalDepth <= 10:
		return 0.70 // 7:3 at 10x (70%) passes
	case totalDepth <= 20:
		return 0.75
	case totalDepth <= 30:
		return 0.80
	default:
		return parentAllelePurityBase
	}
}

// hasUsableAlleleSignal reports whether a sample's allele-depth data at this
// site is trustworthy enough to use, without requiring a confident discrete
// genotype call. Pooled bulk samples are frequently emitted with a missing
// or low-confidence GT by callers designed around single-diploid-sample
// genotype likelihoods, even when read support itself is perfectly good —
// gating on GT there discards real BSA-seq signal along with noise.
//
// When AD/AO fields are unavailable, the function no longer falls back to a
// strict GT check (which would reject most pooled samples with GT='./.').
// Instead, it accepts the sample if FORMAT DP indicates any coverage at all,
// trusting that the site has *some* signal even if we can't precisely
// quantify allele fractions.
func hasUsableAlleleSignal(v *vcfgo.Variant, idx, targetAlt int) bool {
	if idx < 0 || idx >= len(v.Samples) {
		return false
	}
	s := v.Samples[idx]

	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR != nil || errA != nil || targetAlt-1 < 0 || targetAlt-1 >= len(altDepths) {
		// No usable AD/AO — but don't reject pooled samples just because
		// they lack allele-depth fields or have GT='./.'. If the sample has
		// any depth at all, assume there's usable signal. Only fall back to
		// the strict GT check if there's no depth information whatsoever.
		if s.DP > 0 {
			return true
		}
		// No DP either — try to salvage from whatever we parsed.
		if errR == nil && refDepth > 0 {
			return true
		}
		if errA == nil {
			for _, d := range altDepths {
				if d > 0 {
					return true
				}
			}
		}
		// Truly no information — fall back to GT as last resort.
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

	threshold := otherAlleleThreshold(total)
	return float64(other)/float64(total) <= threshold
}

// parentAllele determines which allele (0 = ref, targetAlt = alt) predominates
// in a parent sample using allele read depths, falling back to the genotype
// call only when depth data isn't usable. confident is false when the
// parent's reads are too evenly split between alleles to trust — i.e. the
// parent looks genuinely heterozygous/segregating at this site, which makes
// it uninformative for BSA-seq regardless of how strict the filter is.
//
// The purity threshold is depth-dependent: at low coverage (common in BSA-seq
// parents), a flat 85% cutoff would discard samples like 4 REF / 1 ALT
// (80%) that are almost certainly homozygous. The depth-aware threshold
// allows such cases while still rejecting truly ambiguous parents.
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
			threshold := parentAlleleThreshold(total)

			switch {
			// A single discordant read is tolerated regardless of the
			// fraction it implies — but only when the other allele still
			// clearly outnumbers it (guards against e.g. a lone alt read
			// with zero ref support being misread as "ref-confident").
			// At low parent depth (e.g. the default 5x), this stops one
			// sequencing error from pushing a genuinely homozygous parent
			// below the flat purity threshold.
			case refFrac >= threshold || (nonRef <= 1 && refDepth > nonRef):
				return 0, true
			case refFrac <= 1-threshold || (refDepth <= 1 && nonRef > refDepth):
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

// parentAlleleScore returns a continuous score from -1 to +1 indicating how
// strongly the parent supports REF (-1) or ALT (+1). A score of 0 means
// perfectly ambiguous. This is used for the "one parent uncertain" logic
// where we want to know the direction of the uncertain parent's signal
// even if it doesn't meet the confidence threshold.
func parentAlleleScore(v *vcfgo.Variant, idx, targetAlt int) float64 {
	if idx < 0 || idx >= len(v.Samples) {
		return 0
	}
	s := v.Samples[idx]

	refDepth, errR := safeSampleRefDepth(s)
	altDepths, errA := safeSampleAltDepths(s)
	if errR != nil || errA != nil || targetAlt-1 < 0 || targetAlt-1 >= len(altDepths) {
		return 0
	}

	altDepth := altDepths[targetAlt-1]
	total := refDepth + altDepth
	if total <= 0 {
		return 0
	}
	return float64(altDepth)/float64(total) - float64(refDepth)/float64(total)
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

// checkBulkSignalAndDepth validates that a bulk sample has usable allele
// signal and meets the depth threshold. Returns false if the bulk fails
// either check.
func checkBulkSignalAndDepth(v *vcfgo.Variant, idx int, minDepth int) bool {
	if !hasUsableAlleleSignal(v, idx, -1) {
		return false
	}
	return effectiveDP(v.Samples[idx]) >= minDepth
}

// checkParentDepth validates that a parent sample meets the depth threshold.
func checkParentDepth(v *vcfgo.Variant, idx int, minDepth int) bool {
	return effectiveDP(v.Samples[idx]) >= minDepth
}

// parentsAreInformative checks whether two parents show different alleles
// with sufficient confidence. Unlike the original strict check (both parents
// must be confident AND different), this allows one parent to be "uncertain"
// as long as it doesn't clearly conflict with the confident parent's allele.
//
// The rationale: BSA-seq parents often have low coverage (5-10x). If Parent
// 1 clearly shows REF (6:1, 86%) but Parent 2 is borderline (3:1, 75% at
// the depth-adjusted threshold), the variant is still likely informative —
// Parent 2 probably IS homozygous REF but just doesn't meet the strict
// threshold due to low depth + one sequencing error. The key is that Parent
// 2 doesn't clearly show ALT, so we're not misassigning the cross structure.
//
// Returns: (informative, highParentAllele, lowParentAllele)
//   - informative: true if the parents can be used for this variant
//   - highParentAllele / lowParentAllele: the assigned alleles (0=REF, targetAlt=ALT)
//     Only valid when informative=true.
func parentsAreInformative(v *vcfgo.Variant, highParentIdx, lowParentIdx, targetAlt int) (bool, int, int) {
	hpAllele, hpOK := parentAllele(v, highParentIdx, targetAlt)
	lpAllele, lpOK := parentAllele(v, lowParentIdx, targetAlt)

	// Best case: both confident and different — definitely informative.
	if hpOK && lpOK && hpAllele != lpAllele {
		return true, hpAllele, lpAllele
	}

	// Both confident but same allele — not informative for this cross.
	if hpOK && lpOK && hpAllele == lpAllele {
		return false, 0, 0
	}

	// One confident, one uncertain — check if the uncertain parent is at
	// least CONSISTENT with being the opposite allele (i.e., doesn't clearly
	// show the same allele as the confident parent).
	if hpOK && !lpOK {
		lpScore := parentAlleleScore(v, lowParentIdx, targetAlt)
		// If the uncertain parent's score is in the same direction as the
		// confident parent, they likely share the same allele — reject.
		// hpAllele=0 (REF) means we want lpScore < 0 (REF-biased) to reject.
		// hpAllele=targetAlt (ALT) means we want lpScore > 0 to reject.
		if hpAllele == 0 && lpScore < -0.1 {
			return false, 0, 0 // Low parent also looks REF
		}
		if hpAllele == targetAlt && lpScore > 0.1 {
			return false, 0, 0 // Low parent also looks ALT
		}
		// Uncertain parent doesn't clearly conflict — assign it the opposite
		// allele and mark as informative.
		if hpAllele == 0 {
			return true, 0, targetAlt
		}
		return true, targetAlt, 0
	}

	if !hpOK && lpOK {
		hpScore := parentAlleleScore(v, highParentIdx, targetAlt)
		if lpAllele == 0 && hpScore < -0.1 {
			return false, 0, 0 // High parent also looks REF
		}
		if lpAllele == targetAlt && hpScore > 0.1 {
			return false, 0, 0 // High parent also looks ALT
		}
		if lpAllele == 0 {
			return true, targetAlt, 0
		}
		return true, 0, targetAlt
	}

	// Both uncertain — too risky to use. However, if their scores point in
	// opposite directions with reasonable magnitude, we might still trust it.
	hpScore := parentAlleleScore(v, highParentIdx, targetAlt)
	lpScore := parentAlleleScore(v, lowParentIdx, targetAlt)

	// If both point the same way, they're not informative.
	if (hpScore > 0.1 && lpScore > 0.1) || (hpScore < -0.1 && lpScore < -0.1) {
		return false, 0, 0
	}

	// If they point in opposite directions with at least moderate confidence,
	// accept them. This handles cases where both are just below the
	// confidence threshold but clearly support different alleles.
	if hpScore > 0.15 && lpScore < -0.15 {
		return true, targetAlt, 0
	}
	if hpScore < -0.15 && lpScore > 0.15 {
		return true, 0, targetAlt
	}

	return false, 0, 0
}

func bsaSeqFilterAllele(v *vcfgo.Variant, cfg utils.AnalysisConfig, bsaType string, targetAlt int) bool {
	switch bsaType {
	case "2b":
		// Bulks only — require at least one bulk to have good signal,
		// and any bulk with good signal must meet its depth threshold.
		// Previously, BOTH bulks had to pass, which discarded variants
		// where one bulk had a sequencing dropout or paralog issue.
		highOK := hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt)
		lowOK := hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt)

		if !highOK && !lowOK {
			return false // Neither bulk has usable signal
		}

		// Any bulk with usable signal must meet its depth threshold
		if highOK && effectiveDP(v.Samples[cfg.HighBulkIdx]) < cfg.HighBulkDepth {
			return false
		}
		if lowOK && effectiveDP(v.Samples[cfg.LowBulkIdx]) < cfg.LowBulkDepth {
			return false
		}
		return true

	case "2p2b":
		// 2 parents 2 bulks filter
		highOK := hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt)
		lowOK := hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt)

		if !highOK && !lowOK {
			return false
		}
		if highOK && effectiveDP(v.Samples[cfg.HighBulkIdx]) < cfg.HighBulkDepth {
			return false
		}
		if lowOK && effectiveDP(v.Samples[cfg.LowBulkIdx]) < cfg.LowBulkDepth {
			return false
		}

		informative, _, _ := parentsAreInformative(v, cfg.HighParentIdx, cfg.LowParentIdx, targetAlt)
		if !informative {
			return false
		}

		return checkParentDepth(v, cfg.HighParentIdx, cfg.HighParentDepth) &&
			checkParentDepth(v, cfg.LowParentIdx, cfg.LowParentDepth)

	case "2plb":
		// 2 parents low bulk filter
		if !hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt) {
			return false
		}
		if effectiveDP(v.Samples[cfg.LowBulkIdx]) < cfg.LowBulkDepth {
			return false
		}

		informative, _, _ := parentsAreInformative(v, cfg.HighParentIdx, cfg.LowParentIdx, targetAlt)
		if !informative {
			return false
		}

		return checkParentDepth(v, cfg.HighParentIdx, cfg.HighParentDepth) &&
			checkParentDepth(v, cfg.LowParentIdx, cfg.LowParentDepth)

	case "2phb":
		// 2 parents high bulk filter
		if !hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt) {
			return false
		}
		if effectiveDP(v.Samples[cfg.HighBulkIdx]) < cfg.HighBulkDepth {
			return false
		}

		informative, _, _ := parentsAreInformative(v, cfg.HighParentIdx, cfg.LowParentIdx, targetAlt)
		if !informative {
			return false
		}

		return checkParentDepth(v, cfg.HighParentIdx, cfg.HighParentDepth) &&
			checkParentDepth(v, cfg.LowParentIdx, cfg.LowParentDepth)

	case "hp2b":
		// high parent 2 bulks filter
		highOK := hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt)
		lowOK := hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt)

		// High parent must have usable signal
		if !hasUsableAlleleSignal(v, cfg.HighParentIdx, targetAlt) {
			return false
		}

		if !highOK && !lowOK {
			return false
		}
		if highOK && effectiveDP(v.Samples[cfg.HighBulkIdx]) < cfg.HighBulkDepth {
			return false
		}
		if lowOK && effectiveDP(v.Samples[cfg.LowBulkIdx]) < cfg.LowBulkDepth {
			return false
		}

		return checkParentDepth(v, cfg.HighParentIdx, cfg.HighParentDepth)

	case "lp2b":
		// low parent 2 bulks filter
		highOK := hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt)
		lowOK := hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt)

		// Low parent must have usable signal
		if !hasUsableAlleleSignal(v, cfg.LowParentIdx, targetAlt) {
			return false
		}

		if !highOK && !lowOK {
			return false
		}
		if highOK && effectiveDP(v.Samples[cfg.HighBulkIdx]) < cfg.HighBulkDepth {
			return false
		}
		if lowOK && effectiveDP(v.Samples[cfg.LowBulkIdx]) < cfg.LowBulkDepth {
			return false
		}

		return checkParentDepth(v, cfg.LowParentIdx, cfg.LowParentDepth)

	case "hphb":
		// high parent high bulk filter
		if !hasUsableAlleleSignal(v, cfg.HighParentIdx, targetAlt) {
			return false
		}
		if !hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt) {
			return false
		}
		return checkParentDepth(v, cfg.HighParentIdx, cfg.HighParentDepth) &&
			effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth

	case "hplb":
		// high parent low bulk filter
		if !hasUsableAlleleSignal(v, cfg.HighParentIdx, targetAlt) {
			return false
		}
		if !hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt) {
			return false
		}
		return checkParentDepth(v, cfg.HighParentIdx, cfg.HighParentDepth) &&
			effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth

	case "lphb":
		// low parent high bulk filter
		if !hasUsableAlleleSignal(v, cfg.LowParentIdx, targetAlt) {
			return false
		}
		if !hasUsableAlleleSignal(v, cfg.HighBulkIdx, targetAlt) {
			return false
		}
		return checkParentDepth(v, cfg.LowParentIdx, cfg.LowParentDepth) &&
			effectiveDP(v.Samples[cfg.HighBulkIdx]) >= cfg.HighBulkDepth

	case "lplb":
		// low parent low bulk filter
		if !hasUsableAlleleSignal(v, cfg.LowParentIdx, targetAlt) {
			return false
		}
		if !hasUsableAlleleSignal(v, cfg.LowBulkIdx, targetAlt) {
			return false
		}
		return checkParentDepth(v, cfg.LowParentIdx, cfg.LowParentDepth) &&
			effectiveDP(v.Samples[cfg.LowBulkIdx]) >= cfg.LowBulkDepth

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

			// Skip variants with no real ALT alleles (spanning deletions only,
			// symbolic alleles only, missing only). This catches more cases than
			// the original check which only looked at len==1 with specific values.
			if !hasAnyRealAlt(v) {
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
	//skipped := int(skippedCount.Load())

	color.Blue(
		"\nHard filtering complete:\n%d variant records read\n%d variants passed\n%d variants rejected\n",
		original, passed, original-passed,
	)

	return passedVariants, original, passed, nil
}
