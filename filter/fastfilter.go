package filter

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/schollz/progressbar/v3"
)

// ── BSA-seq constants ─────────────────────────────────────────────────────────

// parentAllelePurity is the fraction of a parent's reads that must support one
// allele for the parent to count as homozygous for it.
const parentAllelePurity = 0.85

// secondAltMaxFraction caps the share of a sample's reads on ALT alleles other
// than the target. Reads on other ALTs are excluded from the SNP-index, so a
// large share would make the biallelic view unrepresentative.
const secondAltMaxFraction = 0.05

// minNonParentReads is the minimum number of reads of the allele the lone
// parent does NOT carry that a single bulk must show in one-parent one-bulk
// modes when --require-nonparent-allele is set (see bsaCheck).
const minNonParentReads = 2

// pooledMAFAutoFactor sets the automatic segregation floor: the pooled minor
// allele frequency of both bulks must be at least this fraction of the
// population's expected minor allele frequency (0.2 for F2).
const pooledMAFAutoFactor = 0.4

// Site is one biallelic site that passed filtering, in compact form. Samples
// are in keepIndices order (the order RawStats expects).
type Site struct {
	Chrom     string
	Pos       int64
	Ref, Alt  string
	HighIsAlt bool // the high-parent allele is the ALT allele (false: REF). In 2b mode it is set only by --ref-parent low
	Samples   [4]SampleAD
	NSamples  uint8
}

// roles holds each role's position within keepIndices (-1 = absent).
type roles struct{ hp, lp, hb, lb int }

func modeRoles(bsaType string) (hp, lp, hb, lb bool, err error) {
	switch bsaType {
	case "2p2b":
		return true, true, true, true, nil
	case "2phb":
		return true, true, true, false, nil
	case "2plb":
		return true, true, false, true, nil
	case "hp2b":
		return true, false, true, true, nil
	case "lp2b":
		return false, true, true, true, nil
	case "hphb":
		return true, false, true, false, nil
	case "hplb":
		return true, false, false, true, nil
	case "lphb":
		return false, true, true, false, nil
	case "lplb":
		return false, true, false, true, nil
	case "2b":
		return false, false, true, true, nil
	}
	return false, false, false, false, fmt.Errorf("unsupported bsaseq type %q", bsaType)
}

// Filterer applies GATK hard filters and the BSA-seq filters to VCF lines. It
// is immutable after construction and safe for concurrent use.
type Filterer struct {
	h            *Header
	cfg          utils.AnalysisConfig
	hf           utils.HardFilterConfig
	bsaType      string
	keep         []int
	r            roles
	colToKeep    map[int]int // VCF sample index → position in keep
	maxCol       int
	minPooledMAF float64
	split        bool
}

// NewFilterer validates the mode and sample indices against the header.
func NewFilterer(h *Header, cfg utils.AnalysisConfig, hf utils.HardFilterConfig, bsaType string, keepIndices []int) (*Filterer, error) {
	hp, lp, hb, lb, err := modeRoles(bsaType)
	if err != nil {
		return nil, err
	}
	if len(keepIndices) == 0 || len(keepIndices) > 4 {
		return nil, fmt.Errorf("need 1–4 sample indices, got %d", len(keepIndices))
	}
	f := &Filterer{h: h, cfg: cfg, hf: hf, bsaType: bsaType, keep: keepIndices,
		colToKeep: map[int]int{}, maxCol: -1, split: cfg.SplitMultiallelic}
	for p, idx := range keepIndices {
		if idx < 0 || idx >= len(h.Samples) {
			return nil, fmt.Errorf("sample index %d out of range (VCF has %d samples)", idx, len(h.Samples))
		}
		if _, dup := f.colToKeep[idx]; dup {
			return nil, fmt.Errorf("sample index %d used for two roles", idx)
		}
		f.colToKeep[idx] = p
		if idx > f.maxCol {
			f.maxCol = idx
		}
	}
	pos := func(present bool, idx int, role string) (int, error) {
		if !present {
			return -1, nil
		}
		if p, ok := f.colToKeep[idx]; ok {
			return p, nil
		}
		return -1, fmt.Errorf("%s sample index %d is not among the kept samples %v", role, idx, keepIndices)
	}
	if f.r.hp, err = pos(hp, cfg.HighParentIdx, "high parent"); err != nil {
		return nil, err
	}
	if f.r.lp, err = pos(lp, cfg.LowParentIdx, "low parent"); err != nil {
		return nil, err
	}
	if f.r.hb, err = pos(hb, cfg.HighBulkIdx, "high bulk"); err != nil {
		return nil, err
	}
	if f.r.lb, err = pos(lb, cfg.LowBulkIdx, "low bulk"); err != nil {
		return nil, err
	}

	// Segregation floor for bulk modes without two informative parents.
	switch {
	case cfg.MinPooledMAF >= 0:
		f.minPooledMAF = cfg.MinPooledMAF
	default: // < 0: automatic, from the population's expected minor allele frequency
		p0, err := utils.ExpectedAF(cfg.Population)
		if err != nil {
			p0 = 0.5
		}
		f.minPooledMAF = pooledMAFAutoFactor * min(p0, 1-p0)
	}
	return f, nil
}

// ── GATK hard filter ──────────────────────────────────────────────────────────

// siteInfo holds the site-level annotations used by the GATK hard filter.
type siteInfo struct {
	qd, fs, sor, mq, mqrs, rprs                   float64
	hasQD, hasFS, hasSOR, hasMQ, hasMQRS, hasRPRS bool
}

func parseSiteInfo(info []byte) siteInfo {
	var si siteInfo
	if len(info) == 1 && info[0] == '.' {
		return si
	}
	for len(info) > 0 {
		j := bytes.IndexByte(info, ';')
		kv := info
		if j >= 0 {
			kv, info = info[:j], info[j+1:]
		} else {
			info = nil
		}
		e := bytes.IndexByte(kv, '=')
		if e < 0 {
			continue
		}
		v := kv[e+1:]
		switch string(kv[:e]) {
		case "QD":
			si.qd, si.hasQD = parseFloatField(v)
		case "FS":
			si.fs, si.hasFS = parseFloatField(v)
		case "SOR":
			si.sor, si.hasSOR = parseFloatField(v)
		case "MQ":
			si.mq, si.hasMQ = parseFloatField(v)
		case "MQRankSum":
			si.mqrs, si.hasMQRS = parseFloatField(v)
		case "ReadPosRankSum":
			si.rprs, si.hasRPRS = parseFloatField(v)
		}
	}
	return si
}

// hardPass applies the GATK germline hard-filter recommendations, with GATK
// VariantFiltration semantics: a record fails when an annotation is present and
// beyond its threshold; absent annotations pass (e.g. rank-sums are not
// computed when no sample is heterozygous). SNPs (and same-length MNPs) use the
// SNP set; any length change uses the INDEL set. Light mode keeps FS/SOR but
// drops the rank-sum tests, which are unreliable near allele fixation in pools.
func (f *Filterer) hardPass(qual float64, hasQual bool, si *siteInfo, isSNP bool) bool {
	hf := &f.hf
	strict := !hf.LightFilter
	if isSNP {
		switch {
		case hasQual && qual < hf.SNP_QUAL_Min,
			si.hasQD && si.qd < hf.SNP_QD_Min,
			si.hasMQ && si.mq < hf.SNP_MQ_Min,
			si.hasFS && si.fs > hf.SNP_FS_Max,
			si.hasSOR && si.sor > hf.SNP_SOR_Max,
			strict && si.hasMQRS && si.mqrs < hf.SNP_MQRankSum_Min,
			strict && si.hasRPRS && si.rprs < hf.SNP_ReadPosRankSum_Min:
			return false
		}
		return true
	}
	switch {
	case hasQual && qual < hf.INDEL_QUAL_Min,
		si.hasQD && si.qd < hf.INDEL_QD_Min,
		si.hasFS && si.fs > hf.INDEL_FS_Max,
		si.hasSOR && si.sor > hf.INDEL_SOR_Max,
		strict && si.hasRPRS && si.rprs < hf.INDEL_ReadPosRankSum_Min:
		return false
	}
	return true
}

// ── BSA-seq filter ────────────────────────────────────────────────────────────

// drop is the reason a biallelic record was rejected (dropNone = kept). The
// filter counts each reason so the cost of every rule is visible.
type drop int

const (
	dropNone drop = iota
	dropMalformed
	dropFilterField
	dropGATK
	dropNoAD
	dropMinDepth
	dropMaxDepth
	dropParentGQ
	dropOtherAlt
	dropParentNotHom
	dropParentsSame
	dropSegregation
	dropNonParent
	nDrops
)

var dropNames = [nDrops]string{
	dropMalformed:    "malformed record or missing sample column",
	dropFilterField:  "FILTER not PASS/.",
	dropGATK:         "GATK hard filter (QUAL/QD/FS/SOR/MQ/rank-sums)",
	dropNoAD:         "bulk without usable allele depths (AD or RO/AO)",
	dropMinDepth:     "below minimum depth (-p/-b)",
	dropMaxDepth:     "above maximum depth cap",
	dropParentGQ:     "parent GQ below --min-GQ",
	dropOtherAlt:     "other-ALT reads > 5% of allele depth",
	dropParentNotHom: "parent not homozygous for REF or this ALT",
	dropParentsSame:  "parents carry the same allele",
	dropSegregation:  "segregation floor (pooled bulk minor AF, --min-pooled-maf)",
	dropNonParent:    "bulk lacks the non-parent allele (--require-nonparent-allele)",
}

// bulkCheck: allele depths present, informative depth (REF + target ALT, the
// reads the SNP-index uses) within bounds, and little support for other ALTs.
// No GQ gate: in a pool, GQ tracks allele frequency and would bias the SNP-index.
func bulkCheck(s *SampleAD, minDepth, maxDepth int) drop {
	switch d := s.Ref + s.Alt; {
	case !s.HasAD:
		return dropNoAD
	case d == 0 || int(d) < minDepth:
		return dropMinDepth
	case maxDepth > 0 && int(max(s.AlleleDepth(), s.DP)) > maxDepth:
		return dropMaxDepth
	case float64(s.Other) > secondAltMaxFraction*float64(s.AlleleDepth()):
		return dropOtherAlt
	}
	return dropNone
}

// parentCheck: depth (allele depth, or DP when AD is absent) within bounds, GQ
// above the floor when present (a parent is a real genotype), and — for a lone
// parent, whose reads are not otherwise cross-checked — little other-ALT
// support. Returns the parent's homozygous allele (0 = REF, 1 = target ALT).
func parentCheck(s *SampleAD, minDepth, maxDepth, minGQ int, lone bool) (int, drop) {
	d := s.AlleleDepth()
	if !s.HasAD {
		d = max(s.DP, 0)
	}
	switch {
	case int(d) < minDepth:
		return 0, dropMinDepth
	case maxDepth > 0 && int(max(s.AlleleDepth(), s.DP)) > maxDepth:
		return 0, dropMaxDepth
	case minGQ > 0 && s.GQ >= 0 && int(s.GQ) < minGQ:
		return 0, dropParentGQ
	case lone && s.HasAD && float64(s.Other) > secondAltMaxFraction*float64(s.AlleleDepth()):
		return 0, dropOtherAlt
	}
	if a, ok := parentAllele(s); ok {
		return a, dropNone
	}
	return 0, dropParentNotHom
}

// parentAllele returns the allele a parent is homozygous for (0 = REF,
// 1 = target ALT) and whether the call is confident. With allele depths, the
// TARGET allele's own share must reach parentAllelePurity (reads on other ALTs
// count against both). A parent with AD present but zero informative reads is
// not confident. Without AD, a homozygous REF or target-ALT genotype is used.
func parentAllele(s *SampleAD) (int, bool) {
	if s.HasAD {
		total := s.AlleleDepth()
		if total == 0 {
			return 0, false
		}
		ref, alt := float64(s.Ref)/float64(total), float64(s.Alt)/float64(total)
		switch {
		case ref >= parentAllelePurity || (s.Alt+s.Other <= 1 && s.Ref > s.Alt+s.Other):
			return 0, true
		case alt >= parentAllelePurity || (s.Ref+s.Other <= 1 && s.Alt > s.Ref+s.Other):
			return 1, true
		}
		return 0, false
	}
	if s.GTMissing > 0 || s.GTOther > 0 {
		return 0, false
	}
	switch {
	case s.GTRef > 0 && s.GTAlt == 0:
		return 0, true
	case s.GTAlt > 0 && s.GTRef == 0:
		return 1, true
	}
	return 0, false
}

// bsaCheck applies the mode-specific BSA-seq filter to one biallelic view and
// returns the polarity (is the high-parent allele the ALT?) or a drop reason.
//
//   - Bulks: bulkCheck.
//   - Two-parent modes: both parents pass parentCheck and differ.
//   - One-parent modes: the lone parent passes parentCheck (incl. other-ALT).
//   - 2b / hp2b / lp2b: without two informative parents nothing proves the site
//     segregates in the cross, so the pooled minor allele frequency of the two
//     bulks must reach the segregation floor (removes fixed differences from the
//     reference and sites polymorphic only in other samples of the VCF).
//   - One-parent one-bulk modes, only with --require-nonparent-allele: the bulk
//     must show the allele the parent does NOT carry at least minNonParentReads
//     times. Off by default: it removes non-segregating sites but also true QTL
//     markers where the bulk is fixed for the parent's allele.
func (f *Filterer) bsaCheck(sa []SampleAD) (highIsAlt bool, why drop) {
	c := &f.cfg
	if f.r.hb >= 0 {
		if why = bulkCheck(&sa[f.r.hb], c.HighBulkDepth, c.HighBulkMaxDepth); why != dropNone {
			return
		}
	}
	if f.r.lb >= 0 {
		if why = bulkCheck(&sa[f.r.lb], c.LowBulkDepth, c.LowBulkMaxDepth); why != dropNone {
			return
		}
	}

	twoBulks := f.r.hb >= 0 && f.r.lb >= 0
	twoParents := f.r.hp >= 0 && f.r.lp >= 0
	ha, la := -1, -1
	if f.r.hp >= 0 {
		if ha, why = parentCheck(&sa[f.r.hp], c.HighParentDepth, c.HighParentMaxDepth, c.MinGQ, !twoParents); why != dropNone {
			return
		}
	}
	if f.r.lp >= 0 {
		if la, why = parentCheck(&sa[f.r.lp], c.LowParentDepth, c.LowParentMaxDepth, c.MinGQ, !twoParents); why != dropNone {
			return
		}
	}
	if twoParents && ha == la {
		return false, dropParentsSame
	}
	// The high-parent allele: the high parent's, else the one the low parent lacks.
	highIsAlt = ha == 1 || (ha < 0 && la == 0)
	// No parent sample (2b): polarity is known only if the reference genome is a
	// parent. --ref-parent low means REF is the low parent's allele, so ALT is
	// the high one; high and none leave REF as the reference allele for SI.
	if ha < 0 && la < 0 && c.RefParent == "low" {
		highIsAlt = true
	}

	if twoBulks && !twoParents && f.minPooledMAF > 0 {
		h, l := &sa[f.r.hb], &sa[f.r.lb]
		af := float64(h.Alt+l.Alt) / float64(h.Ref+h.Alt+l.Ref+l.Alt)
		if af < f.minPooledMAF || af > 1-f.minPooledMAF {
			return false, dropSegregation
		}
	}
	if !twoBulks && !twoParents && c.RequireNonParentAllele {
		b := max(f.r.hb, f.r.lb)
		nonParent := sa[b].Alt
		if max(ha, la) == 1 { // the lone parent carries ALT
			nonParent = sa[b].Ref
		}
		if nonParent < minNonParentReads {
			return false, dropNonParent
		}
	}
	return highIsAlt, dropNone
}

// ── per-record processing ─────────────────────────────────────────────────────

type passRec struct {
	chrom    string
	beg, end int // 0-based half-open span of REF, for tabix
	n        int // bytes of the rendered record (including '\n')
}

type batchResult struct {
	seq      int
	pass     []byte // rendered passing records
	recs     []passRec
	rej      []byte // rendered rejected records
	sites    []Site
	entering int // biallelic records entering the filter (after splitting)
	drops    [nDrops]int
	nLines   int
}

// workerState holds per-goroutine scratch buffers and string interning.
type workerState struct {
	cols    [][]byte
	alts    [][]byte
	keys    [][]byte
	raw     [4]rawSample
	adBufs  [4][]int32
	gtBufs  [4][]int8
	views   [4]SampleAD
	chroms  map[string]string
	strings map[string]string
}

func newWorkerState() *workerState {
	return &workerState{cols: make([][]byte, 9), chroms: map[string]string{}, strings: map[string]string{}}
}

func (ws *workerState) intern(m map[string]string, b []byte) string {
	if s, ok := m[string(b)]; ok {
		return s
	}
	s := string(b)
	if len(m) < 1<<16 {
		m[s] = s
	}
	return s
}

// appendSubset renders the record with only the kept sample columns.
func (f *Filterer) appendSubset(dst []byte, cols [][]byte, ws *workerState) []byte {
	for i := 0; i < 9; i++ {
		if i > 0 {
			dst = append(dst, '\t')
		}
		dst = append(dst, cols[i]...)
	}
	for p := range f.keep {
		dst = append(dst, '\t')
		dst = append(dst, ws.raw[p].text...)
	}
	return append(dst, '\n')
}

// processLine filters one VCF data line and appends its outputs to res.
// Malformed lines are rejected (written to the rejected VCF unchanged).
func (f *Filterer) processLine(line []byte, ws *workerState, res *batchResult) {
	cols := ws.cols
	rest, ok := splitTabs(line, cols)
	var pos int64
	if ok {
		pos, _ = strconv.ParseInt(string(cols[1]), 10, 64)
	}
	if !ok || rest == nil || pos < 1 {
		res.rej = append(append(res.rej, line...), '\n')
		res.entering++
		res.drops[dropMalformed]++
		return
	}

	// Real ALT alleles. Records without one (".", <NON_REF>, "*") carry no
	// variant to analyse and are skipped, as before.
	ws.alts = ws.alts[:0]
	altField := cols[4]
	for {
		j := bytes.IndexByte(altField, ',')
		a := altField
		if j >= 0 {
			a = altField[:j]
		}
		ws.alts = append(ws.alts, a)
		if j < 0 {
			break
		}
		altField = altField[j+1:]
	}
	nReal := 0
	for _, a := range ws.alts {
		if isRealAlt(a) {
			nReal++
		}
	}
	if nReal == 0 {
		return
	}
	nAlleles := len(ws.alts) + 1

	// Kept sample columns only; stop after the last one needed.
	fi := parseFormat(cols[8], ws.keys)
	ws.keys = fi.keys
	found := 0
	for col := 0; col <= f.maxCol; col++ {
		j := bytes.IndexByte(rest, '\t')
		field := rest
		if j >= 0 {
			field = rest[:j]
		}
		if p, ok := f.colToKeep[col]; ok {
			ws.adBufs[p], ws.gtBufs[p] = parseSample(field, &fi, nAlleles, &ws.raw[p], ws.adBufs[p], ws.gtBufs[p])
			found++
		}
		if j < 0 {
			break
		}
		rest = rest[j+1:]
	}
	if found != len(f.keep) {
		res.rej = append(append(res.rej, line...), '\n')
		res.entering++
		res.drops[dropMalformed]++
		return
	}

	// Site-level checks shared by every ALT.
	filterOK := len(cols[6]) == 0 || bytes.Equal(cols[6], []byte("PASS")) || bytes.Equal(cols[6], []byte("."))
	qual, hasQual := parseFloatField(cols[5])
	var si siteInfo
	if filterOK {
		si = parseSiteInfo(cols[7])
	}
	refLen := len(cols[3])

	var chrom string
	anyPass := false
	views := ws.views[:len(f.keep)]
	for altPos, a := range ws.alts {
		if !isRealAlt(a) {
			continue
		}
		res.entering++
		if !filterOK {
			res.drops[dropFilterField]++
			continue
		}
		if !f.hardPass(qual, hasQual, &si, len(a) == refLen) {
			res.drops[dropGATK]++
			continue
		}
		for p := range views {
			views[p] = biallelicView(&ws.raw[p], altPos+1)
		}
		highIsAlt, why := f.bsaCheck(views)
		if why != dropNone {
			res.drops[why]++
			continue
		}

		if chrom == "" {
			chrom = ws.intern(ws.chroms, cols[0])
		}
		site := Site{Chrom: chrom, Pos: pos, Ref: ws.intern(ws.strings, cols[3]), Alt: ws.intern(ws.strings, a),
			HighIsAlt: highIsAlt, NSamples: uint8(len(views))}
		copy(site.Samples[:], views)
		res.sites = append(res.sites, site)

		start := len(res.pass)
		if f.split && len(ws.alts) > 1 {
			res.pass = appendSplitRecord(res.pass, cols, ws.alts, ws.raw[:len(f.keep)], &fi, f.h, altPos)
			res.pass = append(res.pass, '\n')
		} else {
			res.pass = f.appendSubset(res.pass, cols, ws)
		}
		res.recs = append(res.recs, passRec{chrom: chrom, beg: int(pos) - 1, end: int(pos) - 1 + refLen, n: len(res.pass) - start})
		anyPass = true
		if !f.split {
			break // without splitting, the first passing ALT represents the record
		}
	}
	if !anyPass {
		res.rej = f.appendSubset(res.rej, cols, ws)
	}
}

// ── filter profile reporting ──────────────────────────────────────────────────

// FilterProfile records which quality annotations the caller wrote into the VCF.
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
func (p FilterProfile) HasGATKAnnotations() bool { return p.HasQD || p.HasFS || p.HasSOR || p.HasMQ }

// DetectFilterProfile inspects the header for available annotations and makes a
// best-effort guess at the calling tool (for reporting).
func DetectFilterProfile(h *Header) FilterProfile {
	p := FilterProfile{Caller: "unknown"}
	if h == nil {
		return p
	}
	p.HasQD, p.HasFS, p.HasSOR, p.HasMQ = h.HasInfo("QD"), h.HasInfo("FS"), h.HasInfo("SOR"), h.HasInfo("MQ")
	p.HasMQRankSum, p.HasReadPosRankSum = h.HasInfo("MQRankSum"), h.HasInfo("ReadPosRankSum")
	_, p.HasGQ = h.Formats["GQ"]

	var meta strings.Builder
	for _, l := range h.Meta {
		if strings.HasPrefix(l, "##contig") || strings.HasPrefix(l, "##INFO") || strings.HasPrefix(l, "##FORMAT") {
			continue
		}
		meta.WriteString(strings.ToLower(l))
	}
	m := meta.String()
	switch {
	case strings.Contains(m, "deepvariant"):
		p.Caller = "deepvariant"
	case strings.Contains(m, "haplotypecaller") || strings.Contains(m, "gatk"):
		p.Caller = "gatk"
	case p.HasGATKAnnotations():
		p.Caller = "gatk-like"
	case p.HasGQ:
		p.Caller = "deepvariant-like"
	}
	return p
}

func (f *Filterer) reportActiveFilters(p FilterProfile) {
	hf, cfg := f.hf, f.cfg
	mode := "GATK best-practice (strict)"
	if hf.LightFilter {
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
	if !hf.LightFilter {
		add(p.HasMQRankSum, "MQRankSum")
		add(p.HasReadPosRankSum, "ReadPosRankSum")
	}
	color.Blue("Active annotation filters: %s", strings.Join(active, ", "))
	if !p.HasGATKAnnotations() {
		color.Yellow("No GATK annotations found (e.g. DeepVariant): relying on FILTER=PASS, QUAL and depth.")
	}
	if cfg.MinGQ > 0 {
		if p.HasGQ {
			color.Blue("BSA gate: min GQ = %d (parent samples only; bulks are depth-only to avoid AF bias)", cfg.MinGQ)
		} else {
			color.Yellow("min-GQ=%d requested but GQ is not present in this VCF; skipping GQ filter.", cfg.MinGQ)
		}
	}
	if f.split {
		color.Blue("Multi-allelic decomposition: enabled (other-ALT reads are kept for the tolerance check)")
	} else {
		color.Yellow("Multi-allelic decomposition: disabled (first passing ALT represents each record)")
	}
	color.Blue("BSA gate: bulk depth = REF+ALT allele depth; other-ALT tolerance = %.0f%% of allele depth", secondAltMaxFraction*100)
	if f.r.hb >= 0 && f.r.lb >= 0 && !(f.r.hp >= 0 && f.r.lp >= 0) {
		if f.minPooledMAF > 0 {
			color.Blue("BSA gate: segregation — pooled bulk minor allele frequency ≥ %.3f", f.minPooledMAF)
		} else {
			color.Yellow("BSA gate: segregation filter disabled (--min-pooled-maf 0)")
		}
	}
	if f.r.hp < 0 && f.r.lp < 0 {
		switch cfg.RefParent {
		case "high", "low":
			color.Blue("Polarity (2b): reference genome is the %s parent (--ref-parent), so ΔSI is signed", cfg.RefParent)
		default:
			color.Yellow("Polarity (2b): unknown (--ref-parent none). SI is the REF fraction, whose sign vs the parents " +
				"varies by SNP, so 2b uses |ΔSI|: QTLs can be found but not which bulk carries which parent's allele. " +
				"BRM uses the polarity-free ΔSI² form; no HMM LOD.")
		}
	} else if cfg.RefParent == "high" || cfg.RefParent == "low" {
		color.Yellow("--ref-parent %s ignored: polarity comes from the parent sample(s)", cfg.RefParent)
	}
	if !(f.r.hb >= 0 && f.r.lb >= 0) && !(f.r.hp >= 0 && f.r.lp >= 0) {
		if cfg.RequireNonParentAllele {
			color.Blue("BSA gate: one parent + one bulk — the bulk must carry ≥%d reads of the non-parent allele", minNonParentReads)
		} else {
			color.Yellow("BSA note: one parent + one bulk cannot tell a bulk fixed for the parent's allele from a " +
				"non-segregating site; use --require-nonparent-allele to drop such sites")
		}
	}
	if cfg.HighBulkMaxDepth > 0 || cfg.LowBulkMaxDepth > 0 || cfg.HighParentMaxDepth > 0 || cfg.LowParentMaxDepth > 0 {
		color.Blue("BSA gate: max-depth caps — HB:%d LB:%d HP:%d LP:%d (0 = off)",
			cfg.HighBulkMaxDepth, cfg.LowBulkMaxDepth, cfg.HighParentMaxDepth, cfg.LowParentMaxDepth)
	}
}

// ── pipeline ──────────────────────────────────────────────────────────────────

type lineBatch struct {
	seq    int
	arena  []byte
	starts []int // line i is arena[starts[i]:starts[i+1]]
}

const (
	batchMaxLines = 8192
	batchMaxBytes = 8 << 20
)

// writeHeader writes the original "##" lines and a #CHROM line naming only the
// kept samples.
func writeHeader(w io.Writer, h *Header, keep []int) error {
	var b bytes.Buffer
	for _, l := range h.Meta {
		b.WriteString(l)
		b.WriteByte('\n')
	}
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for _, idx := range keep {
		b.WriteByte('\t')
		b.WriteString(h.Samples[idx])
	}
	b.WriteByte('\n')
	_, err := w.Write(b.Bytes())
	return err
}

type indexEntry struct {
	chrom      string
	beg, end   int
	begV, endV vpos
}

// HardFilterVcf streams cfg.VCF through the GATK hard filter and the BSA-seq
// filter. It writes GoBSAseq.<mode>.hardfiltered.vcf.gz (+ .tbi) with the kept
// samples, and GoBSAseq.<mode>.lowqual.vcf.gz with records where no ALT passed.
// It returns the passing sites in input order, the number of biallelic records
// that entered the filter (after splitting), and the number that passed.
//
// Throughput comes from: multi-threaded BGZF decompression; parsing in parallel
// workers; parsing only QUAL/FILTER/the needed INFO keys and the kept sample
// columns; compact per-site output instead of parsed VCF objects; and parallel
// BGZF compression whose block offsets are known exactly, so the tabix index is
// built in the same pass.
func HardFilterVcf(cfg utils.AnalysisConfig, hfcfg utils.HardFilterConfig, bsaType string, keepIndices []int) ([]Site, int, int, error) {
	br, closeIn, err := openVCFStream(cfg.VCF)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("open VCF %q: %w", cfg.VCF, err)
	}
	defer closeIn()
	h, err := ParseHeader(br)
	if err != nil {
		return nil, 0, 0, err
	}
	f, err := NewFilterer(h, cfg, hfcfg, bsaType, keepIndices)
	if err != nil {
		return nil, 0, 0, err
	}
	f.reportActiveFilters(DetectFilterProfile(h))
	return f.run(br)
}

func (f *Filterer) run(br *bufio.Reader) ([]Site, int, int, error) {
	cfg := f.cfg
	statsDir := filepath.Join(cfg.OutputDir, "stats")
	if err := os.MkdirAll(statsDir, 0775); err != nil {
		return nil, 0, 0, err
	}
	hfPath := filepath.Join(statsDir, fmt.Sprintf("GoBSAseq.%s.hardfiltered.vcf.gz", f.bsaType))
	badPath := filepath.Join(statsDir, fmt.Sprintf("GoBSAseq.%s.lowqual.vcf.gz", f.bsaType))

	hfFile, err := os.Create(hfPath)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("create hard-filtered VCF: %w", err)
	}
	defer hfFile.Close()
	badFile, err := os.Create(badPath)
	if err != nil {
		return nil, 0, 0, fmt.Errorf("create rejected VCF: %w", err)
	}
	defer badFile.Close()

	hfOut := bufio.NewWriterSize(hfFile, 1<<20)
	badOut := bufio.NewWriterSize(badFile, 1<<20)
	hfW := newBGZFWriter(hfOut, -1)
	badW := newBGZFWriter(badOut, -1)
	closeWriters := func() error {
		e1 := hfW.Close()
		e2 := badW.Close()
		e3 := hfOut.Flush()
		e4 := badOut.Flush()
		for _, e := range []error{e1, e2, e3, e4} {
			if e != nil {
				return e
			}
		}
		return nil
	}
	if err := writeHeader(hfW, f.h, f.keep); err != nil {
		closeWriters()
		return nil, 0, 0, err
	}
	if err := writeHeader(badW, f.h, f.keep); err != nil {
		closeWriters()
		return nil, 0, 0, err
	}

	workers := runtime.GOMAXPROCS(0)
	batches := make(chan *lineBatch, workers*2)
	results := make(chan *batchResult, workers*2)
	arenas := sync.Pool{New: func() any { return make([]byte, 0, batchMaxBytes+(1<<16)) }}

	// Reader: cut the stream into line batches.
	var readErr error
	go func() {
		defer close(batches)
		seq := 0
		for {
			b := &lineBatch{seq: seq, arena: arenas.Get().([]byte)[:0], starts: make([]int, 0, batchMaxLines+1)}
			var err error
			for len(b.starts) < batchMaxLines && len(b.arena) < batchMaxBytes {
				var start int
				b.arena, start, err = readLineInto(br, b.arena)
				if len(b.arena) > start && b.arena[start] != '#' {
					b.starts = append(b.starts, start)
				} else {
					b.arena = b.arena[:start]
				}
				if err != nil {
					break
				}
			}
			if len(b.starts) > 0 {
				b.starts = append(b.starts, len(b.arena))
				batches <- b
				seq++
			}
			if err != nil {
				if err != io.EOF {
					readErr = err
				}
				return
			}
		}
	}()

	// Workers.
	var wg sync.WaitGroup
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			ws := newWorkerState()
			for b := range batches {
				res := &batchResult{seq: b.seq, nLines: len(b.starts) - 1}
				for i := 0; i+1 < len(b.starts); i++ {
					f.processLine(b.arena[b.starts[i]:b.starts[i+1]], ws, res)
				}
				arenas.Put(b.arena[:0])
				results <- res
			}
		}()
	}
	go func() { wg.Wait(); close(results) }()

	// Writer: restore input order, write outputs, record exact index positions.
	bar := progressbar.Default(-1, "Filtering variants")
	var (
		sites    []Site
		entries  []indexEntry
		entering int
		drops    [nDrops]int
		pending  = map[int]*batchResult{}
		next     = 0
		writeErr error
	)
	emit := func(res *batchResult) {
		sites = append(sites, res.sites...)
		entering += res.entering
		for i, n := range res.drops {
			drops[i] += n
		}
		off := 0
		for _, r := range res.recs {
			rec := res.pass[off : off+r.n]
			off += r.n
			hfW.Reserve(len(rec))
			beg := hfW.Pos()
			hfW.Write(rec)
			entries = append(entries, indexEntry{chrom: r.chrom, beg: r.beg, end: r.end, begV: beg, endV: hfW.Pos()})
		}
		if len(res.rej) > 0 {
			badW.Write(res.rej)
		}
		_ = bar.Add(res.nLines)
	}
	for res := range results {
		pending[res.seq] = res
		for {
			r, ok := pending[next]
			if !ok {
				break
			}
			delete(pending, next)
			emit(r)
			next++
		}
	}
	_ = bar.Finish()
	fmt.Println()

	if err := closeWriters(); err != nil {
		writeErr = err
	}
	if readErr != nil {
		return nil, 0, 0, fmt.Errorf("VCF read error: %w", readErr)
	}
	if writeErr != nil {
		return nil, 0, 0, fmt.Errorf("write filtered VCF: %w", writeErr)
	}

	if err := writeTabix(hfPath+".tbi", hfW, entries); err != nil {
		os.Remove(hfPath + ".tbi")
		color.Yellow("Tabix index not written: %v", err)
	}

	color.Blue("Hard filtering complete:\n%d variant records entered the filter (after multi-allelic splitting).\n%d passed.\n%d rejected:",
		entering, len(sites), entering-len(sites))
	for i, n := range drops {
		if n > 0 {
			color.Blue("  %10d  %s", n, dropNames[i])
		}
	}
	return sites, entering, len(sites), nil
}
