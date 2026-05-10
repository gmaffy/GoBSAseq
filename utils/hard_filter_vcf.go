package utils

import (
	"fmt"
	"os"
	"runtime"
	"strings"
	"sync"

	"github.com/biogo/hts/bgzf"
	"github.com/brentp/vcfgo"
	"github.com/schollz/progressbar/v3"
)

// ---------------------------------------------------------------------------
// GATK hard-filtering thresholds (germline short variants best practices).
//
// SNPs:
//   QD < 2.0  |  QUAL < 30.0  |  FS > 60.0  |  SOR > 3.0
//   MQ < 40.0  |  MQRankSum < -12.5  |  ReadPosRankSum < -8.0
//
// INDELs:
//   QD < 2.0  |  QUAL < 30.0  |  FS > 200.0  |  ReadPosRankSum < -20.0
//
// Callers of HardFilterVcf are responsible for populating HardFilterConfig
// with these values. Call DefaultHardFilterConfig() to get a correctly
// initialised struct.
// ---------------------------------------------------------------------------

// DefaultHardFilterConfig returns a HardFilterConfig pre-populated with the
// canonical GATK hard-filtering thresholds for germline short variants.
// Override individual fields only when your data distribution justifies it.
func DefaultHardFilterConfig() HardFilterConfig {
	return HardFilterConfig{
		// SNP thresholds
		SNP_QUAL_Min:           30.0,
		SNP_QD_Min:             2.0,
		SNP_FS_Max:             60.0,
		SNP_SOR_Max:            3.0,
		SNP_MQ_Min:             40.0,
		SNP_MQRankSum_Min:      -12.5,
		SNP_ReadPosRankSum_Min: -8.0,
		// INDEL thresholds
		INDEL_QUAL_Min:           30.0,
		INDEL_QD_Min:             2.0,
		INDEL_FS_Max:             200.0,
		INDEL_ReadPosRankSum_Min: -20.0,
	}
}

// ---------------------------------------------------------------------------
// "bad sample string" — root cause and fix
//
// GATK HaplotypeCaller emits FORMAT fields like PGT and PID with '.' for
// non-phased sites, producing sample columns such as:
//
//   GT:AD:DP:GQ:PGT:PID:PL   0/1:18,8:26:99:.:.:186,0,704
//
// vcfgo's parseSample sees the '.' token for PGT (declared as String in the
// header) and emits "bad sample string" because its internal type switch does
// not treat '.' as a universal missing-value sentinel for String fields.
//
// The VCF 4.2 spec §1.4.2 states: "If the field is not present (missing), a
// dot ('.') should be used." A '.' in any typed FORMAT column is valid.
//
// Fix strategy — two complementary steps:
//
//  1. stripProblematicFormatFields: removes PGT and PID from the header's
//     SampleFormats map before any parsing begins. vcfgo then treats those
//     columns as untyped strings stored in SampleGenotype.Fields, which
//     never triggers the type-assertion path that generates the error.
//
//  2. parseSampleTolerant: after rdr.Read(), calls rdr.Error(), inspects
//     the error message, and calls rdr.Clear() for known-benign "bad sample
//     string" errors. For any other error it returns it to the caller.
//     This is the pattern recommended by vcfgo's own documentation:
//     "after every rdr.Read() call, the caller can check rdr.Error() and
//     get feedback on the errors without stopping execution."
// ---------------------------------------------------------------------------

// stripProblematicFormatFields removes FORMAT fields that GATK emits with '.'
// missing values and that vcfgo's typed parser cannot handle gracefully.
// Removing them from the header causes vcfgo to store the raw token in
// SampleGenotype.Fields instead of attempting type conversion.
//
// This must be called before the first rdr.Read().

// getFloat safely retrieves a float-typed INFO field.  vcfgo can return the
// value as float32 or float64 depending on the header declaration, so we
// handle both.
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

func stripProblematicFormatFields(hdr *vcfgo.Header) {
	for _, id := range []string{"PGT", "PID"} {
		delete(hdr.SampleFormats, id)
	}
}

// isBadSampleStringError returns true for the specific vcfgo error that
// arises from GATK's legitimate use of '.' in typed FORMAT columns.
func isBadSampleStringError(err error) bool {
	return err != nil && strings.Contains(err.Error(), "bad sample string")
}

// ---------------------------------------------------------------------------
// Variant-type helpers
// ---------------------------------------------------------------------------

// isHomozygous returns true when every allele in the genotype is the same
// non-missing call. An empty or all-missing GT is treated as non-homozygous.
func isHomozygous(gt []int) bool {
	if len(gt) == 0 {
		return false
	}
	first := gt[0]
	if first < 0 { // missing allele
		return false
	}
	for _, a := range gt[1:] {
		if a < 0 || a != first {
			return false
		}
	}
	return true
}

// IsSNP returns true when the REF and every non-symbolic ALT allele are
// exactly one base long. Symbolic alleles (".", "*", "<...>") are skipped so
// that a mixed record with one symbolic entry is not misclassified.
func IsSNP(v *vcfgo.Variant) bool {
	if len(v.Ref()) != 1 {
		return false
	}
	for _, alt := range v.Alt() {
		if IsSymbolic(alt) {
			continue
		}
		if len(alt) != 1 {
			return false
		}
	}
	return true
}

// IsIndel returns true when any non-symbolic ALT allele differs in length
// from the REF allele (the standard VCF left-aligned indel representation).
func IsIndel(v *vcfgo.Variant) bool {
	refLen := len(v.Ref())
	for _, alt := range v.Alt() {
		if IsSymbolic(alt) {
			continue
		}
		if len(alt) != refLen {
			return true
		}
	}
	return false
}

// IsSymbolic reports whether an ALT token is symbolic (not a sequence allele).
func IsSymbolic(alt string) bool {
	return alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<')
}

// ---------------------------------------------------------------------------
// GATK quality filters
// ---------------------------------------------------------------------------

// isQualitySNP applies GATK SNP hard-filtering thresholds.
// An annotation that is absent from the INFO field is treated as passing
// (consistent with GATK VariantFiltration behaviour: missing → not filtered).
func isQualitySNP(v *vcfgo.Variant, cfg HardFilterConfig) bool {
	// QUAL is a core VCF field; always present (may be '.' → 0 in vcfgo).
	if float64(v.Quality) < cfg.SNP_QUAL_Min {
		return false
	}
	// INFO annotations: skip the test when the field is absent.
	if qd, ok := getFloat(v, "QD"); ok && qd < cfg.SNP_QD_Min {
		return false
	}
	if fs, ok := getFloat(v, "FS"); ok && fs > cfg.SNP_FS_Max {
		return false
	}
	if sor, ok := getFloat(v, "SOR"); ok && sor > cfg.SNP_SOR_Max {
		return false
	}
	if mq, ok := getFloat(v, "MQ"); ok && mq < cfg.SNP_MQ_Min {
		return false
	}
	if mqrs, ok := getFloat(v, "MQRankSum"); ok && mqrs < cfg.SNP_MQRankSum_Min {
		return false
	}
	if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < cfg.SNP_ReadPosRankSum_Min {
		return false
	}
	return true
}

// isQualityINDEL applies GATK INDEL hard-filtering thresholds.
// MQ, MQRankSum, and SOR are deliberately omitted for INDELs — this matches
// the official GATK recommendation, which does not apply those filters to
// indels because their distributions differ significantly from SNPs.
func isQualityINDEL(v *vcfgo.Variant, cfg HardFilterConfig) bool {
	if float64(v.Quality) < cfg.INDEL_QUAL_Min {
		return false
	}
	if qd, ok := getFloat(v, "QD"); ok && qd < cfg.INDEL_QD_Min {
		return false
	}
	if fs, ok := getFloat(v, "FS"); ok && fs > cfg.INDEL_FS_Max {
		return false
	}
	if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < cfg.INDEL_ReadPosRankSum_Min {
		return false
	}
	return true
}

// ---------------------------------------------------------------------------
// BSA-seq filters
// ---------------------------------------------------------------------------

// passesTwoBulkBSAseqFilters checks that:
//  1. The site is biallelic (exactly one ALT allele).
//  2. All four samples (high parent, low parent, high bulk, low bulk) have
//     fully called, non-missing genotypes.
//  3. Both parents are homozygous and differ from each other — i.e. the
//     parents are fixed for alternate alleles, which is the prerequisite for
//     meaningful BSA-seq allele-frequency contrasts.
//  4. Per-sample depth meets the configured minimums.
func passesTwoBulkBSAseqFilters(v *vcfgo.Variant, cfg AnalysisConfig) bool {
	// BSA-seq requires exactly one "real" (non-symbolic) ALT allele for
	// unambiguous allele-frequency calculation between the two bulks.
	// We allow <NON_REF> as an additional symbolic allele (common in gVCFs).
	realAltIdx := -1
	realAltCount := 0
	for i, alt := range v.Alt() {
		if !IsSymbolic(alt) {
			realAltIdx = i
			realAltCount++
		}
	}

	if realAltCount != 1 {
		return false
	}

	// The 'real' allele index in vcfgo's 1-based allele numbering is realAltIdx + 1.
	targetAltAllele := realAltIdx + 1

	indices := []int{cfg.HighParentIdx, cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx}
	for _, idx := range indices {
		if idx < 0 || idx >= len(v.Samples) {
			return false
		}
		s := v.Samples[idx]
		if s == nil || len(s.GT) == 0 {
			return false
		}
		for _, allele := range s.GT {
			if allele < 0 { // missing allele ('.')
				return false
			}
			// Genotypes should only carry the REF (0) or the one real ALT allele.
			// Sites where parents or bulks carry a third allele (symbolic or
			// otherwise) are excluded.
			if allele != 0 && allele != targetAltAllele {
				return false
			}
		}
	}

	hpS := v.Samples[cfg.HighParentIdx]
	lpS := v.Samples[cfg.LowParentIdx]
	hbS := v.Samples[cfg.HighBulkIdx]
	lbS := v.Samples[cfg.LowBulkIdx]

	// Parents must be homozygous and fixed for different alleles.
	if !isHomozygous(hpS.GT) || !isHomozygous(lpS.GT) {
		return false
	}
	if hpS.GT[0] == lpS.GT[0] {
		// Both parents carry the same allele — not informative for BSA-seq.
		return false
	}

	// Enforce minimum sequencing depth per sample.
	return hpS.DP >= cfg.HighParentDepth &&
		lpS.DP >= cfg.LowParentDepth &&
		hbS.DP >= cfg.HighBulkDepth &&
		lbS.DP >= cfg.LowBulkDepth
}

// ---------------------------------------------------------------------------
// filterResult carries a variant and its keep/discard decision through the
// pipeline. The original input index is retained so results can be written in
// the same order as the input (required for bgzf files that will be indexed
// with tabix).
// ---------------------------------------------------------------------------
type filterResult struct {
	v    *vcfgo.Variant
	idx  int // original read order
	keep bool
}

// ---------------------------------------------------------------------------
// HardFilterVcf reads variants from rdr, applies GATK hard-filtering and
// BSA-seq filters, writes passing variants to hardFilteredVcfPath and failing
// variants to badVcfPath, and returns the slice of passing variants together
// with the total and passing counts.
//
// Call stripProblematicFormatFields(rdr.Header) before constructing rdr if
// the VCF was produced by GATK HaplotypeCaller (handles PGT/PID '.' values).
//
// Output files are written in input order (required for downstream tabix
// indexing). Filtering itself is parallelised across all available CPUs.
// ---------------------------------------------------------------------------
func HardFilterVcf(
	rdr *vcfgo.Reader,
	hardFilteredVcfPath string,
	badVcfPath string,
	cfg AnalysisConfig,
	hfcfg HardFilterConfig,
) ([]*vcfgo.Variant, int, int, error) {

	// ------------------------------------------------------------------ //
	// Strip FORMAT fields that GATK emits with '.' and that vcfgo cannot
	// parse without emitting "bad sample string" errors.
	// Must happen before the first Read() call.
	// ------------------------------------------------------------------ //
	stripProblematicFormatFields(rdr.Header)

	// ------------------------------------------------------------------ //
	// Open output files
	// ------------------------------------------------------------------ //

	hfFile, err := os.Create(hardFilteredVcfPath)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("create hard-filtered VCF: %w", err)
	}
	defer hfFile.Close()

	hfBgzf := bgzf.NewWriter(hfFile, 1)
	hfWriter, err := vcfgo.NewWriter(hfBgzf, rdr.Header)
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
	badWriter, err := vcfgo.NewWriter(badBgzf, rdr.Header)
	if err != nil {
		hfBgzf.Close()
		badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("create rejected VCF writer: %w", err)
	}

	// ------------------------------------------------------------------ //
	// Concurrency setup
	//
	// Architecture:
	//   main goroutine  →  variantChan  →  worker pool  →  resultChan
	//                                                     ↓
	//                                              order buffer
	//                                                     ↓
	//                                           writer goroutine
	//
	// The order buffer ensures variants are written in input order so the
	// output VCF can be tabix-indexed without re-sorting.
	// ------------------------------------------------------------------ //

	const chanBuf = 10_000

	variantChan := make(chan filterResult, chanBuf)
	resultChan := make(chan filterResult, chanBuf)

	var workerWG sync.WaitGroup
	var writerWG sync.WaitGroup

	bar := progressbar.Default(-1, "Hard filtering variants")

	// ------------------------------------------------------------------ //
	// Writer goroutine — collects results, reorders them, writes to disk.
	// ------------------------------------------------------------------ //
	var (
		hardFilteredVariants []*vcfgo.Variant
		hardFiltered         int
		writeErr             error
	)

	writerWG.Add(1)
	go func() {
		defer writerWG.Done()

		// pending holds out-of-order results until their turn arrives.
		pending := make(map[int]filterResult)
		nextWrite := 0

		for res := range resultChan {
			pending[res.idx] = res

			// Drain any consecutive results that are now ready.
			for {
				r, ok := pending[nextWrite]
				if !ok {
					break
				}
				delete(pending, nextWrite)
				nextWrite++

				if r.keep {
					hardFilteredVariants = append(hardFilteredVariants, r.v)
					hfWriter.WriteVariant(r.v)
					hardFiltered++
				} else {
					badWriter.WriteVariant(r.v)
				}
				_ = bar.Add(1)
			}
		}
	}()

	// ------------------------------------------------------------------ //
	// Worker pool — applies filters concurrently.
	// ------------------------------------------------------------------ //
	numWorkers := runtime.NumCPU()
	for i := 0; i < numWorkers; i++ {
		workerWG.Add(1)
		go func() {
			defer workerWG.Done()
			for item := range variantChan {
				v := item.v
				keep := false
				switch {
				case IsSNP(v):
					keep = isQualitySNP(v, hfcfg) && passesTwoBulkBSAseqFilters(v, cfg)
				case IsIndel(v):
					keep = isQualityINDEL(v, hfcfg) && passesTwoBulkBSAseqFilters(v, cfg)
					// Mixed / MNP / symbolic-only records: discard.
					// They are neither SNPs nor indels and cannot be meaningfully
					// hard-filtered with the scalar thresholds above.
				}
				resultChan <- filterResult{v: v, idx: item.idx, keep: keep}
			}
		}()
	}

	// ------------------------------------------------------------------ //
	// Main read loop
	//
	// vcfgo accumulates errors internally; rdr.Error() returns them all as
	// a single VCFError after each Read(). We must call rdr.Clear() after
	// inspecting the error or the same errors will be re-reported on every
	// subsequent call to rdr.Error().
	//
	// "bad sample string" errors are benign: vcfgo still populates GT, AD,
	// DP, GQ correctly — the error only signals that a later typed field
	// (e.g. PGT, PID) contained '.' which vcfgo couldn't coerce to the
	// declared type. stripProblematicFormatFields() above eliminates the
	// most common source, but we keep the fallback guard here in case other
	// VCF producers emit similar patterns for fields we didn't pre-strip.
	// ------------------------------------------------------------------ //
	original := 0
	skipped := 0

	for {
		v := rdr.Read()
		if v == nil {
			break
		}

		if err := rdr.Error(); err != nil {
			if isBadSampleStringError(err) {
				// Clear the accumulated errors so they don't compound.
				// The variant itself is usable — GT/AD/DP/GQ are populated.
				rdr.Clear()
			} else {
				// Unexpected parse error — stop processing.
				close(variantChan)
				workerWG.Wait()
				close(resultChan)
				writerWG.Wait()
				_ = hfBgzf.Close()
				_ = badBgzf.Close()
				return nil, 0, 0, fmt.Errorf("VCF parse error at line %d: %w", v.LineNumber, err)
			}
		}

		// Discard gVCF reference-confidence blocks: ALT == '<NON_REF>' or
		// '.' only. These are not variant sites and must not enter filtering.
		alts := v.Alt()
		if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
			skipped++
			continue
		}

		variantChan <- filterResult{v: v, idx: original}
		original++
	}

	// Final check for errors after the loop (in case Read() returned nil due to error).
	if err := rdr.Error(); err != nil && !isBadSampleStringError(err) {
		close(variantChan)
		workerWG.Wait()
		close(resultChan)
		writerWG.Wait()
		_ = hfBgzf.Close()
		_ = badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("VCF read error: %w", err)
	}

	close(variantChan)
	workerWG.Wait()
	close(resultChan)
	writerWG.Wait()
	_ = bar.Finish()

	if writeErr != nil {
		_ = hfBgzf.Close()
		_ = badBgzf.Close()
		return nil, 0, 0, writeErr
	}

	// Flush vcfgo writers before closing the underlying bgzf streams.
	// vcfgo.Writer writes synchronously so no additional flush is needed;
	// closing the bgzf layer finalises the BGZF EOF block.
	if err = hfBgzf.Close(); err != nil {
		_ = badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("close hard-filtered bgzf: %w", err)
	}
	if err = badBgzf.Close(); err != nil {
		return nil, 0, 0, fmt.Errorf("close rejected bgzf: %w", err)
	}

	fmt.Printf(
		"Hard filtering complete: %d variant records read (%d gVCF ref blocks skipped) → %d passed, %d rejected\n",
		original, skipped, hardFiltered, original-hardFiltered,
	)

	return hardFilteredVariants, original, hardFiltered, nil
}
