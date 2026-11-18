package filter

import (
	"bytes"
	"math"
	"strconv"
	"strings"
)

// SampleAD is one sample's data for one biallelic (REF vs target ALT) view of
// a record. Depths come from AD (or RO/AO); Other is the summed depth of every
// ALT allele other than the target, which multi-allelic splitting would
// otherwise silently discard.
type SampleAD struct {
	Ref, Alt, Other int32
	DP              int32 // FORMAT/DP, -1 if absent
	GQ              int32 // FORMAT/GQ, -1 if absent
	HasAD           bool  // allele depths were present and well-formed
	// Genotype allele counts in the biallelic view: REF, target ALT, any other
	// ALT, and missing ('.') alleles.
	GTRef, GTAlt, GTOther, GTMissing uint8
}

// AlleleDepth is REF + target ALT + other ALT reads.
func (s SampleAD) AlleleDepth() int32 { return s.Ref + s.Alt + s.Other }

// GT returns the biallelic genotype as allele indices (0 = REF, 1 = target
// ALT, 2 = another ALT, -1 = missing), unphased and sorted.
func (s SampleAD) GT() []int {
	gt := make([]int, 0, int(s.GTRef)+int(s.GTAlt)+int(s.GTOther)+int(s.GTMissing))
	for range s.GTMissing {
		gt = append(gt, -1)
	}
	for range s.GTRef {
		gt = append(gt, 0)
	}
	for range s.GTAlt {
		gt = append(gt, 1)
	}
	for range s.GTOther {
		gt = append(gt, 2)
	}
	return gt
}

// rawSample is a parsed sample column (all alleles). Slices point into
// per-worker scratch space and are only valid while the record is processed.
type rawSample struct {
	text  []byte
	ad    []int32 // REF + every ALT; nil if absent or malformed
	dp    int32
	gq    int32
	gt    []int8 // allele indices, -1 missing
	ploid int
}

// splitTabs splits the first n tab-separated fields of line into cols and
// returns the remainder (after the n-th tab) or nil when the line has exactly n
// fields.
func splitTabs(line []byte, cols [][]byte) ([]byte, bool) {
	rest := line
	for i := range cols {
		j := bytes.IndexByte(rest, '\t')
		if j < 0 {
			if i == len(cols)-1 {
				cols[i] = rest
				return nil, true
			}
			return nil, false
		}
		cols[i] = rest[:j]
		rest = rest[j+1:]
	}
	return rest, true
}

// atoi32 parses a non-negative decimal integer; ok=false on anything else.
func atoi32(b []byte) (int32, bool) {
	if len(b) == 0 || len(b) > 10 {
		return 0, false
	}
	var n int64
	for _, c := range b {
		if c < '0' || c > '9' {
			return 0, false
		}
		n = n*10 + int64(c-'0')
	}
	if n > math.MaxInt32 {
		return 0, false
	}
	return int32(n), true
}

// parseFloatField parses an INFO/QUAL number. "." and NaN count as missing,
// matching GATK VariantFiltration, which does not filter on absent values.
func parseFloatField(b []byte) (float64, bool) {
	if len(b) == 0 || (len(b) == 1 && b[0] == '.') {
		return 0, false
	}
	v, err := strconv.ParseFloat(string(b), 64)
	if err != nil || math.IsNaN(v) {
		return 0, false
	}
	return v, true
}

// formatIndex records the position of the FORMAT keys the filter reads.
type formatIndex struct {
	keys                   [][]byte
	gt, ad, dp, gq, ro, ao int
}

func parseFormat(format []byte, keys [][]byte) formatIndex {
	fi := formatIndex{gt: -1, ad: -1, dp: -1, gq: -1, ro: -1, ao: -1}
	keys = keys[:0]
	for i := 0; ; i++ {
		j := bytes.IndexByte(format, ':')
		k := format
		if j >= 0 {
			k = format[:j]
		}
		keys = append(keys, k)
		switch string(k) {
		case "GT":
			fi.gt = i
		case "AD":
			fi.ad = i
		case "DP":
			fi.dp = i
		case "GQ":
			fi.gq = i
		case "RO":
			fi.ro = i
		case "AO":
			fi.ao = i
		}
		if j < 0 {
			break
		}
		format = format[j+1:]
	}
	fi.keys = keys
	return fi
}

// parseIntList parses "a,b,c" into dst; ok=false if any element is missing or
// not an integer.
func parseIntList(b []byte, dst []int32) ([]int32, bool) {
	dst = dst[:0]
	for {
		j := bytes.IndexByte(b, ',')
		e := b
		if j >= 0 {
			e = b[:j]
		}
		v, ok := atoi32(e)
		if !ok {
			return dst, false
		}
		dst = append(dst, v)
		if j < 0 {
			return dst, true
		}
		b = b[j+1:]
	}
}

// parseSample fills rs from one sample column. nAlleles = 1 + number of ALTs.
func parseSample(col []byte, fi *formatIndex, nAlleles int, rs *rawSample, adBuf []int32, gtBuf []int8) ([]int32, []int8) {
	rs.text = col
	rs.ad, rs.gt, rs.dp, rs.gq, rs.ploid = nil, nil, -1, -1, 0
	var ro, ao []byte
	for i := 0; ; i++ {
		j := bytes.IndexByte(col, ':')
		v := col
		if j >= 0 {
			v = col[:j]
		}
		switch i {
		case fi.gt:
			gtBuf = gtBuf[:0]
			if !(len(v) == 0 || (len(v) == 1 && v[0] == '.')) {
				start := 0
				for k := 0; k <= len(v); k++ {
					if k == len(v) || v[k] == '/' || v[k] == '|' {
						tok := v[start:k]
						if a, ok := atoi32(tok); ok && a < 127 {
							gtBuf = append(gtBuf, int8(a))
						} else {
							gtBuf = append(gtBuf, -1)
						}
						start = k + 1
					}
				}
			}
			rs.gt = gtBuf
			rs.ploid = len(gtBuf)
		case fi.ad:
			if lst, ok := parseIntList(v, adBuf); ok && len(lst) == nAlleles {
				adBuf = lst
				rs.ad = lst
			}
		case fi.dp:
			if d, ok := atoi32(v); ok {
				rs.dp = d
			}
		case fi.gq:
			if q, ok := atoi32(v); ok {
				rs.gq = q
			} else if f, ok := parseFloatField(v); ok { // some callers write float GQ
				rs.gq = int32(f)
			}
		case fi.ro:
			ro = v
		case fi.ao:
			ao = v
		}
		if j < 0 {
			break
		}
		col = col[j+1:]
	}
	// FreeBayes-style RO/AO when AD is absent.
	if rs.ad == nil && ro != nil && ao != nil {
		r, okR := atoi32(ro)
		if okR {
			lst, okA := parseIntList(ao, adBuf[:0])
			if okA && len(lst) == nAlleles-1 {
				full := append(lst[:0:0], r)
				full = append(full, lst...)
				adBuf = full
				rs.ad = full
			}
		}
	}
	return adBuf, gtBuf
}

// biallelicView projects a raw sample onto REF vs ALT allele alt (1-based).
func biallelicView(rs *rawSample, alt int) SampleAD {
	s := SampleAD{DP: rs.dp, GQ: rs.gq}
	if rs.ad != nil {
		s.HasAD = true
		s.Ref = rs.ad[0]
		for i := 1; i < len(rs.ad); i++ {
			if i == alt {
				s.Alt = rs.ad[i]
			} else {
				s.Other += rs.ad[i]
			}
		}
	}
	for _, a := range rs.gt {
		switch {
		case a < 0:
			s.GTMissing++
		case a == 0:
			s.GTRef++
		case int(a) == alt:
			s.GTAlt++
		default:
			s.GTOther++
		}
	}
	return s
}

// isRealAlt excludes symbolic, spanning-deletion and missing ALT alleles.
func isRealAlt(a []byte) bool {
	return !(len(a) == 0 || (len(a) == 1 && (a[0] == '.' || a[0] == '*')) || a[0] == '<')
}

// ── split-record rendering (bcftools norm -m- conventions) ───────────────────

// remapGTString maps the split ALT to 1, REF to 0 and any other ALT to 0,
// preserving phase separators (bcftools norm -m- behaviour).
func remapGTString(gt []byte, alleleIdx int) []byte {
	if len(gt) == 0 || (len(gt) == 1 && gt[0] == '.') {
		return gt
	}
	out := make([]byte, 0, len(gt))
	start := 0
	want := strconv.Itoa(alleleIdx)
	for k := 0; k <= len(gt); k++ {
		if k == len(gt) || gt[k] == '/' || gt[k] == '|' {
			tok := string(gt[start:k])
			switch tok {
			case ".":
				out = append(out, '.')
			case want:
				out = append(out, '1')
			default:
				out = append(out, '0')
			}
			if k < len(gt) {
				out = append(out, gt[k])
			}
			start = k + 1
		}
	}
	return out
}

// reindexNumbered keeps the elements of a Number=A/R/G value that belong to
// the biallelic REF/ALT[altPos] view. Wrong cardinality, or Number=G for a
// ploidy other than 1 or 2, yields "." (missing) rather than wrong data.
func reindexNumbered(val []byte, number string, altPos, nAlts int) []byte {
	if len(val) == 0 || (len(val) == 1 && val[0] == '.') {
		return val
	}
	parts := strings.Split(string(val), ",")
	missing := []byte(".")
	a := altPos + 1
	switch number {
	case "A":
		if len(parts) != nAlts {
			return missing
		}
		return []byte(parts[altPos])
	case "R":
		if len(parts) != nAlts+1 {
			return missing
		}
		return []byte(parts[0] + "," + parts[a])
	case "G":
		n := nAlts + 1
		switch len(parts) {
		case n * (n + 1) / 2: // diploid: (0,0) (0,a) (a,a)
			return []byte(parts[0] + "," + parts[a*(a+1)/2] + "," + parts[a*(a+1)/2+a])
		case n: // haploid
			return []byte(parts[0] + "," + parts[a])
		}
		return missing
	}
	return val
}

// appendSplitRecord renders the biallelic record for ALT altPos (0-based):
// per-allele INFO keys (Number A/R/G) are dropped, FORMAT A/R/G values are
// reindexed and GT is remapped.
func appendSplitRecord(dst []byte, cols [][]byte, alts [][]byte, samples []rawSample, fi *formatIndex, h *Header, altPos int) []byte {
	nAlts := len(alts)
	for i := 0; i < 4; i++ { // CHROM POS ID REF
		dst = append(dst, cols[i]...)
		dst = append(dst, '\t')
	}
	dst = append(dst, alts[altPos]...)
	dst = append(dst, '\t')
	dst = append(dst, cols[5]...) // QUAL
	dst = append(dst, '\t')
	dst = append(dst, cols[6]...) // FILTER
	dst = append(dst, '\t')

	// INFO without per-allele keys.
	info := cols[7]
	wrote := false
	if !(len(info) == 1 && info[0] == '.') {
		for len(info) > 0 {
			j := bytes.IndexByte(info, ';')
			kv := info
			if j >= 0 {
				kv = info[:j]
				info = info[j+1:]
			} else {
				info = nil
			}
			key := kv
			if e := bytes.IndexByte(kv, '='); e >= 0 {
				key = kv[:e]
			}
			switch h.Infos[string(key)].Number {
			case "A", "R", "G":
				continue
			}
			if wrote {
				dst = append(dst, ';')
			}
			dst = append(dst, kv...)
			wrote = true
		}
	}
	if !wrote {
		dst = append(dst, '.')
	}

	if len(cols[8]) == 0 {
		return dst
	}
	dst = append(dst, '\t')
	dst = append(dst, cols[8]...)

	for si := range samples {
		dst = append(dst, '\t')
		col := samples[si].text
		for k := 0; ; k++ {
			j := bytes.IndexByte(col, ':')
			v := col
			if j >= 0 {
				v = col[:j]
			}
			if k > 0 {
				dst = append(dst, ':')
			}
			switch {
			case k == fi.gt:
				dst = append(dst, remapGTString(v, altPos+1)...)
			case k < len(fi.keys):
				switch num := h.Formats[string(fi.keys[k])].Number; num {
				case "A", "R", "G":
					dst = append(dst, reindexNumbered(v, num, altPos, nAlts)...)
				default:
					dst = append(dst, v...)
				}
			default:
				dst = append(dst, v...)
			}
			if j < 0 {
				break
			}
			col = col[j+1:]
		}
	}
	return dst
}
