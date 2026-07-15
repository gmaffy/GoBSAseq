package filter

import (
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
)

// SplitMultiallelic decomposes a multi-allelic record into one biallelic record
// per real ALT (the pure-Go equivalent of `bcftools norm -m-`). Splitting up
// front means every downstream step sees clean REF-vs-single-ALT records, which
// fixes two problems with analysing multi-allelic sites in place:
//
//   - Variant classification. A site such as A>C,<indel> is neither purely a SNP
//     nor purely an indel; classified as a whole it lands in one bucket and the
//     other allele is filtered with the wrong thresholds (or dropped). Each
//     decomposed record is classified on its own allele.
//   - Allele bookkeeping. TargetAlt becomes trivially 1 and the per-sample AD is
//     reduced to [ref, alt], so the SNP-index denominator is unambiguous.
//
// Per record the function reindexes, by the header-declared Number:
//
//   - FORMAT/GT   → alleles remapped (the split ALT → 1, ref → 0, other ALTs → 0).
//   - FORMAT R    → [ref, alt]        (e.g. AD).
//   - FORMAT A    → [alt].
//   - FORMAT G    → the three diploid genotypes over {ref, alt} (e.g. PL, GL).
//   - INFO A/R/G  → deleted. These per-allele INFO annotations (AC, AF, MLEAC …)
//     are not consumed anywhere in the pipeline — the hard filter reads only
//     Number=1 fields (QD, FS, SOR, MQ, rank-sums) — and vcfgo exposes no way to
//     give each split record its own INFO object, so dropping them yields a valid
//     biallelic record rather than one carrying stale multi-allelic cardinality.
//
// Symbolic and star ALTs (<NON_REF>, *, .) are not emitted as their own records.
// Records with fewer than two real ALTs are returned unchanged.
func SplitMultiallelic(v *vcfgo.Variant) []*vcfgo.Variant {
	realIdx := realAltIndices(v) // one-based ALT numbers of real (non-symbolic) alts
	if len(realIdx) < 2 {
		return []*vcfgo.Variant{v}
	}

	alts := v.Alt()
	nAlts := len(alts)

	// Drop per-allele INFO once on the shared INFO object; every clone then
	// carries the same, valid (Number=1/flag-only) INFO.
	if v.Info() != nil && v.Header != nil {
		for _, key := range v.Info().Keys() {
			if info, ok := v.Header.Infos[key]; ok {
				switch info.Number {
				case "A", "R", "G":
					v.Info().Delete(key)
				}
			}
		}
	}

	out := make([]*vcfgo.Variant, 0, len(realIdx))
	for _, oneBased := range realIdx {
		altPos := oneBased - 1 // zero-based index into alts

		nv := *v // value copy carries unexported fields (Info_, Id_, Header …)
		nv.Alternate = []string{alts[altPos]}

		nv.Samples = make([]*vcfgo.SampleGenotype, len(v.Samples))
		for i, s := range v.Samples {
			nv.Samples[i] = splitSample(s, v.Format, v.Header, altPos, nAlts)
		}
		vc := nv
		out = append(out, &vc)
	}
	return out
}

// splitSample deep-copies a sample genotype and reindexes its FORMAT fields to
// the biallelic view of ALT index altPos (zero-based). alleleIdx is the original
// allele number of that ALT (altPos+1, since allele 0 is REF).
func splitSample(s *vcfgo.SampleGenotype, format []string, hdr *vcfgo.Header, altPos, nAlts int) *vcfgo.SampleGenotype {
	if s == nil {
		return nil
	}
	alleleIdx := altPos + 1

	fields := make(map[string]string, len(s.Fields))
	for k, val := range s.Fields {
		fields[k] = val
	}

	ns := &vcfgo.SampleGenotype{
		Phased: s.Phased,
		DP:     s.DP,
		GQ:     s.GQ,
		MQ:     s.MQ,
		Fields: fields,
	}
	// Remap the parsed integer genotype for in-memory consumers (parentAllele,
	// stats), mirroring the string remap applied to Fields["GT"].
	ns.GT = remapGTInts(s.GT, alleleIdx)

	for _, key := range format {
		val, ok := fields[key]
		if !ok {
			continue
		}
		if key == "GT" {
			fields[key] = remapGTString(val, alleleIdx)
			continue
		}
		number := ""
		if hdr != nil {
			if f, ok := hdr.SampleFormats[key]; ok {
				number = f.Number
			}
		}
		switch number {
		case "A", "R", "G":
			fields[key] = reindexNumbered(val, number, altPos, nAlts)
		}
	}
	return ns
}

// remapGTInts remaps parsed genotype alleles to the biallelic view: the split
// allele becomes 1, REF stays 0, any other ALT collapses to 0, missing stays
// missing. This matches the standard bcftools decomposition of e.g. 1/2 into two
// 0/1 records.
func remapGTInts(gt []int, alleleIdx int) []int {
	out := make([]int, len(gt))
	for i, a := range gt {
		switch {
		case a < 0:
			out[i] = a // missing
		case a == alleleIdx:
			out[i] = 1
		default:
			out[i] = 0
		}
	}
	return out
}

// remapGTString applies the same remapping to the raw GT string while preserving
// the original phase separators ('/' vs '|').
func remapGTString(gt string, alleleIdx int) string {
	if gt == "" || gt == "." {
		return gt
	}
	var b strings.Builder
	allele := strings.Builder{}
	flush := func() {
		if allele.Len() == 0 {
			return
		}
		tok := allele.String()
		allele.Reset()
		switch {
		case tok == ".":
			b.WriteString(".")
		case tok == strconv.Itoa(alleleIdx):
			b.WriteString("1")
		default:
			b.WriteString("0")
		}
	}
	for _, r := range gt {
		if r == '/' || r == '|' {
			flush()
			b.WriteRune(r)
			continue
		}
		allele.WriteRune(r)
	}
	flush()
	return b.String()
}

// reindexNumbered selects the elements of a comma-separated FORMAT value that
// survive decomposition to ALT index altPos (zero-based), according to the VCF
// Number class. On any cardinality mismatch it returns "." (missing) so the
// output stays valid rather than carrying wrong-length data.
func reindexNumbered(val, number string, altPos, nAlts int) string {
	if val == "" || val == "." {
		return val
	}
	parts := strings.Split(val, ",")

	switch number {
	case "A": // one value per ALT
		if len(parts) != nAlts || altPos >= len(parts) {
			return "."
		}
		return parts[altPos]

	case "R": // one value per allele (REF + ALTs)
		if len(parts) != nAlts+1 || altPos+1 >= len(parts) {
			return "."
		}
		return parts[0] + "," + parts[altPos+1]

	case "G": // one value per diploid genotype
		a := altPos + 1
		orders := []int{
			0,               // (0,0)
			a * (a + 1) / 2, // (0,a)
			a*(a+1)/2 + a,   // (a,a)
		}
		picked := make([]string, len(orders))
		for i, o := range orders {
			if o >= len(parts) {
				return "."
			}
			picked[i] = parts[o]
		}
		return strings.Join(picked, ",")
	}
	return val
}
