package utils

import (
	"fmt"
	"io"
	"os"
	"strings"

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

func passesHardFilter(v *vcfgo.Variant, hfcfg HardFilterConfig) bool {
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

func bsaSeqFilter(v *vcfgo.Variant, cfg AnalysisConfig, bsaType int) bool {
	targetAlt := singleRealAlt(v)
	if targetAlt < 0 {
		return false
	}

	switch bsaType {
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
		// 2 parents 1 bulk filter
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

		return hp.DP >= cfg.HighParentDepth && lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth
	case 3:
		// high parent 2 bulks filter
		for _, idx := range []int{cfg.HighParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		hp := v.Samples[cfg.HighParentIdx]
		return hp.DP >= cfg.HighParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	case 4:
		// low parent 1 2 bulks filter
		for _, idx := range []int{cfg.LowParentIdx, cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}
		lp := v.Samples[cfg.LowParentIdx]
		return lp.DP >= cfg.LowParentDepth && v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth
	default:
		// bulks only
		for _, idx := range []int{cfg.HighBulkIdx, cfg.LowBulkIdx} {
			if !sampleHasOnlyRefOrAlt(v, idx, targetAlt) {
				return false
			}
		}

		return v.Samples[cfg.HighBulkIdx].DP >= cfg.HighBulkDepth && v.Samples[cfg.LowBulkIdx].DP >= cfg.LowBulkDepth

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

// ---------------------------------------------------------------------------
// HardFilterVcf
// ---------------------------------------------------------------------------

// HardFilterVcf reads rdr and passes each variant through two sequential gates:
//  1. passesHardFilter  — GATK quality thresholds (SNP / INDEL).
//  2. bsaFilter         — caller-supplied BSA-seq sample filter.
//
// Variants passing both are written (bgzf-compressed) to hardFilteredVcfPath
// with a companion .tbi index. Rejected variants go to badVcfPath.
// Returns the passing variants, total records read, and the passing count.
func HardFilterVcf(rdr *vcfgo.Reader, hardFilteredVcfPath string, badVcfPath string, cfg AnalysisConfig, hfcfg HardFilterConfig, bsaseqType int) ([]*vcfgo.Variant, int, int, error) {

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

	hfIdx := newTabixIndex()

	bar := progressbar.Default(-1, "Hard filtering variants")

	var (
		passedVariants []*vcfgo.Variant
		original       int
		skipped        int
		passed         int
	)

	// ── Main read loop ────────────────────────────────────────────────────────
	//
	// "bad sample string" errors are benign: vcfgo still populates GT, AD, DP,
	// GQ correctly — the error only signals that a later typed field (e.g. PGT,
	// PID) contained '.' which vcfgo couldn't coerce to the declared type.
	// Call rdr.Clear() after inspecting or the same errors compound on every
	// subsequent call.
	for {
		v := rdr.Read()
		if v == nil {
			break
		}

		if err := rdr.Error(); err != nil {
			if strings.Contains(err.Error(), "bad sample string") {
				rdr.Clear() // benign — GT/AD/DP/GQ are still populated
			} else {
				hfBgzf.Close()
				badBgzf.Close()
				return nil, 0, 0, fmt.Errorf("VCF parse error at line %d: %w", v.LineNumber, err)
			}
		}

		// Discard gVCF reference-confidence blocks (ALT == '<NON_REF>' or '.').
		alts := v.Alt()
		if len(alts) == 0 || (len(alts) == 1 && (alts[0] == "<NON_REF>" || alts[0] == ".")) {
			skipped++
			continue
		}

		original++
		_ = bar.Add(1)

		if !passesHardFilter(v, hfcfg) || !bsaSeqFilter(v, cfg, bsaseqType) {
			badWriter.WriteVariant(v)
			continue
		}

		// ── Record tabix offset, write, record end offset ─────────────────────
		blockOffset, _ := hfBgzf.Next()
		startOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffset)}

		hfWriter.WriteVariant(v)

		blockOffsetEnd, _ := hfBgzf.Next()
		endOffset := bgzf.Offset{File: hfCounting.n, Block: uint16(blockOffsetEnd)}

		if err := addTabixRecord(hfIdx, v, bgzf.Chunk{Begin: startOffset, End: endOffset}); err != nil {
			hfBgzf.Close()
			badBgzf.Close()
			return nil, 0, 0, fmt.Errorf("tabix add variant at %s:%d: %w", v.Chromosome, v.Pos, err)
		}

		passedVariants = append(passedVariants, v)
		passed++
	}

	if err := rdr.Error(); err != nil && !strings.Contains(err.Error(), "bad sample string") {
		hfBgzf.Close()
		badBgzf.Close()
		return nil, 0, 0, fmt.Errorf("VCF read error: %w", err)
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

	fmt.Printf(
		"Hard filtering complete: %d variant records read (%d gVCF ref blocks skipped) → %d passed, %d rejected\n",
		original, skipped, passed, original-passed,
	)

	return passedVariants, original, passed, nil
}
