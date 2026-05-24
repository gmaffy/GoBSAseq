package run

import (
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/brentp/vcfgo"
	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/best_of_both/utils"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
	"gonum.org/v1/gonum/stat/distuv"
)

const (
	eps              = 1e-12
	minSNPsPerWindow = 3
	consensusMin     = 3
	defaultBRMAlpha  = 0.05
	afpFloor         = 0.05
)

var statNames = []string{"high_si", "low_si", "delta_si", "g_statistic", "ed4", "lod", "beta_binomial_bf", "brm"}

type mode int

const (
	modeTwoParentsTwoBulks mode = iota
	modeBulksOnly
	modeTwoParentsOneBulk
	modeOneParentOneBulk
)

type marker struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string

	HighDepth sampleDepth
	LowDepth  sampleDepth
	OneDepth  sampleDepth

	HighSI float64
	LowSI  float64
	Delta  float64
	G      float64
	ED4    float64
	LOD    float64
	BF     float64
	BRM    float64
	Depth  int
}

type sampleDepth struct {
	Track int
	Other int
	DP    int
}

type smoothPoint struct {
	Chrom      string
	Start      int64
	End        int64
	Center     int64
	MarkerN    int
	MeanHighDP int
	MeanLowDP  int
	Stats      map[string]float64
	Thresholds map[string]thresholdPair
	Z          map[string]float64
	BRMHit     bool
}

type thresholdPair struct {
	Upper map[float64]float64
	Lower map[float64]float64
}

type qtl struct {
	Chrom     string
	Start     int64
	End       int64
	Stat      string
	Method    string
	Threshold string
	PeakPos   int64
	PeakValue float64
	Markers   int
}

type brmBlock struct {
	Chrom     string
	Start     int64
	End       int64
	PeakPos   int64
	Peak      float64
	Threshold float64
}

type runData struct {
	Mode     mode
	Markers  []marker
	Rejected int
}

func Run(cfg utils.AnalysisConfig, hf utils.HardFilterConfig) error {
	start := time.Now()
	if cfg.VCF == "" {
		return errors.New("variant file is required")
	}
	if cfg.WindowSize <= 0 || cfg.StepSize <= 0 {
		return errors.New("window-size and step-size must be positive")
	}
	if cfg.Rep <= 0 {
		cfg.Rep = 1000
	}
	if len(cfg.Alphas) == 0 {
		cfg.Alphas = []float64{0.05, 0.01}
	}
	if cfg.OutputDir == "" {
		cfg.OutputDir = "."
	}
	if err := os.MkdirAll(filepath.Join(cfg.OutputDir, "plots"), 0o755); err != nil {
		return err
	}

	logStage("Starting best-of-both GoBSAseq")
	color.White("  VCF: %s", cfg.VCF)
	color.White("  Output: %s", cfg.OutputDir)
	color.White("  Window/step: %d/%d bp", cfg.WindowSize, cfg.StepSize)
	color.White("  Simulations: %d; alphas: %s", cfg.Rep, formatAlphas(cfg.Alphas))

	data, err := readMarkers(cfg, hf)
	if err != nil {
		return err
	}
	if len(data.Markers) == 0 {
		return errors.New("no informative markers passed filtering")
	}
	color.White("  Mode: %s", modeName(data.Mode))
	color.White("  Informative markers: %d; rejected/skipped: %d", len(data.Markers), data.Rejected)

	logStage("Calculating BRM seed values")
	sortMarkers(data.Markers)
	assignBRM(data.Markers, cfg.WindowSize)

	logStage("Smoothing with tricube spatial and depth weights")
	points := smoothMarkers(data.Markers, cfg.WindowSize, cfg.StepSize)
	if len(points) == 0 {
		return errors.New("no smoothed windows produced; reduce --window-size/--step-size or check marker density")
	}
	color.White("  Smoothed windows: %d", len(points))

	logStage("Computing depth-adaptive simulated thresholds")
	attachAdaptiveThresholds(points, cfg, data.Mode)

	logStage("Normalising robust z-scores")
	attachRobustZ(points)

	logStage("Detecting BRM blocks and QTLs")
	blocks := detectBRMBlocks(points, cfg)
	markBRMHits(points, blocks)
	permQTLs := detectPermutationQTLs(points, cfg)
	zQTLs := detectZQTLs(points, cfg)
	consensusQTLs := detectConsensusQTLs(points, cfg)
	maxZQTLs := detectMaxZQTLs(points, cfg)
	highConfidence := intersectWithBRM(append(append([]qtl{}, permQTLs...), zQTLs...), blocks, "HighConfidence")
	allQTLs := append(append(append(append(permQTLs, zQTLs...), consensusQTLs...), maxZQTLs...), highConfidence...)
	color.White("  BRM blocks: %d", len(blocks))
	color.White("  QTLs: permutation=%d z=%d consensus=%d maxZ=%d high-confidence=%d", len(permQTLs), len(zQTLs), len(consensusQTLs), len(maxZQTLs), len(highConfidence))

	logStage("Writing tables")
	if err := writeMarkers(filepath.Join(cfg.OutputDir, "markers.tsv"), data.Markers); err != nil {
		return err
	}
	if err := writeSmoothed(filepath.Join(cfg.OutputDir, "smoothed.tsv"), points, cfg.Alphas); err != nil {
		return err
	}
	if err := writeBRM(filepath.Join(cfg.OutputDir, "brm_blocks.tsv"), blocks); err != nil {
		return err
	}
	if err := writeQTLs(filepath.Join(cfg.OutputDir, "qtls.tsv"), allQTLs); err != nil {
		return err
	}

	logStage("Rendering go-echarts plots")
	if err := plotAll(filepath.Join(cfg.OutputDir, "plots"), points, blocks, cfg.Alphas); err != nil {
		return err
	}

	if cfg.SnpEffDB != "" || cfg.Gff != "" || cfg.GeneDesc != "" || cfg.Prg != "" {
		logStage("Gene-space analysis")
		color.Yellow("  Gene-space parameters were supplied. This best_of_both scaffold records QTLs for downstream annotation; direct SnpEff/GeneSpace execution can be enabled after validating local databases.")
	}

	color.Green("Complete in %s", time.Since(start).Round(time.Millisecond))
	return nil
}

func readMarkers(cfg utils.AnalysisConfig, hf utils.HardFilterConfig) (runData, error) {
	logStage("Opening VCF and resolving samples")
	r, closeFn, err := openVCF(cfg.VCF)
	if err != nil {
		return runData{}, err
	}
	defer closeFn()
	rdr, err := vcfgo.NewReader(r, false)
	if err != nil {
		return runData{}, err
	}
	sampleIdx := map[string]int{}
	for i, name := range rdr.Header.SampleNames {
		sampleIdx[name] = i
	}
	m, err := inferMode(cfg)
	if err != nil {
		return runData{}, err
	}
	if err := validateSamples(cfg, m, sampleIdx); err != nil {
		return runData{}, err
	}

	var markers []marker
	rejected := 0
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		if err := rdr.Error(); err != nil {
			if strings.Contains(err.Error(), "bad sample string") {
				rdr.Clear()
				rejected++
				continue
			}
			return runData{}, err
		}
		item, ok := markerFromVariant(v, cfg, hf, m, sampleIdx)
		if !ok {
			rejected++
			continue
		}
		markers = append(markers, item)
	}
	return runData{Mode: m, Markers: markers, Rejected: rejected}, nil
}

func inferMode(cfg utils.AnalysisConfig) (mode, error) {
	hasTwoParents := cfg.HighParentName != "" && cfg.LowParentName != ""
	hasOneParent := cfg.OneParentName != ""
	hasTwoBulks := cfg.HighBulkName != "" && cfg.LowBulkName != ""
	hasOneBulk := cfg.OneBulkName != "" || (cfg.HighBulkName != "" && cfg.LowBulkName == "") || (cfg.LowBulkName != "" && cfg.HighBulkName == "")
	switch {
	case hasTwoParents && hasTwoBulks:
		return modeTwoParentsTwoBulks, nil
	case !hasTwoParents && !hasOneParent && hasTwoBulks:
		return modeBulksOnly, nil
	case hasTwoParents && hasOneBulk:
		return modeTwoParentsOneBulk, nil
	case hasOneParent && hasOneBulk:
		return modeOneParentOneBulk, nil
	default:
		return modeBulksOnly, fmt.Errorf("unsupported setup: select at least one bulk and optionally one/two parents")
	}
}

func validateSamples(cfg utils.AnalysisConfig, m mode, idx map[string]int) error {
	names := requiredNames(cfg, m)
	for _, name := range names {
		if _, ok := idx[name]; !ok {
			return fmt.Errorf("sample %q is not present in VCF header", name)
		}
	}
	return nil
}

func requiredNames(cfg utils.AnalysisConfig, m mode) []string {
	switch m {
	case modeTwoParentsTwoBulks:
		return []string{cfg.HighParentName, cfg.LowParentName, cfg.HighBulkName, cfg.LowBulkName}
	case modeBulksOnly:
		return []string{cfg.HighBulkName, cfg.LowBulkName}
	case modeTwoParentsOneBulk:
		name, _ := oneBulk(cfg)
		return []string{cfg.HighParentName, cfg.LowParentName, name}
	default:
		name, _ := oneBulk(cfg)
		return []string{cfg.OneParentName, name}
	}
}

func markerFromVariant(v *vcfgo.Variant, cfg utils.AnalysisConfig, hf utils.HardFilterConfig, m mode, sampleIdx map[string]int) (marker, bool) {
	realAlt := singleRealAlt(v)
	if realAlt < 0 || !passesHardFilter(v, hf) {
		return marker{}, false
	}
	trackAllele := realAlt
	if m == modeTwoParentsTwoBulks || m == modeTwoParentsOneBulk {
		hp := v.Samples[sampleIdx[cfg.HighParentName]]
		lp := v.Samples[sampleIdx[cfg.LowParentName]]
		if hp.DP < cfg.HighParentDepth || lp.DP < cfg.LowParentDepth || !isHomozygous(hp.GT) || !isHomozygous(lp.GT) || hp.GT[0] == lp.GT[0] {
			return marker{}, false
		}
		if !sampleHasOnlyRefOrAlt(hp, realAlt) || !sampleHasOnlyRefOrAlt(lp, realAlt) {
			return marker{}, false
		}
		trackAllele = hp.GT[0]
	}

	out := marker{Chrom: v.Chromosome, Pos: int64(v.Pos), Ref: v.Ref(), Alt: v.Alt()[realAlt-1]}
	switch m {
	case modeTwoParentsTwoBulks, modeBulksOnly:
		hb := v.Samples[sampleIdx[cfg.HighBulkName]]
		lb := v.Samples[sampleIdx[cfg.LowBulkName]]
		if hb.DP < cfg.HighBulkDepth || lb.DP < cfg.LowBulkDepth || !sampleHasOnlyRefOrAlt(hb, realAlt) || !sampleHasOnlyRefOrAlt(lb, realAlt) {
			return marker{}, false
		}
		var ok bool
		out.HighDepth, ok = orientedDepth(hb, realAlt, trackAllele)
		if !ok {
			return marker{}, false
		}
		out.LowDepth, ok = orientedDepth(lb, realAlt, trackAllele)
		if !ok {
			return marker{}, false
		}
		calcTwoBulkStats(&out)
	default:
		bulkName, minDepth := oneBulk(cfg)
		b := v.Samples[sampleIdx[bulkName]]
		if b.DP < minDepth || !sampleHasOnlyRefOrAlt(b, realAlt) {
			return marker{}, false
		}
		var ok bool
		out.OneDepth, ok = orientedDepth(b, realAlt, trackAllele)
		if !ok {
			return marker{}, false
		}
		calcOneBulkStats(&out)
	}
	return out, true
}

func oneBulk(cfg utils.AnalysisConfig) (string, int) {
	switch {
	case cfg.OneBulkName != "":
		return cfg.OneBulkName, cfg.OneBulkDepth
	case cfg.HighBulkName != "":
		return cfg.HighBulkName, cfg.HighBulkDepth
	default:
		return cfg.LowBulkName, cfg.LowBulkDepth
	}
}

func openVCF(path string) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, func() {}, err
	}
	if strings.HasSuffix(strings.ToLower(path), ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			_ = f.Close()
			return nil, func() {}, err
		}
		return gz, func() { _ = gz.Close(); _ = f.Close() }, nil
	}
	return f, func() { _ = f.Close() }, nil
}

func singleRealAlt(v *vcfgo.Variant) int {
	real := -1
	count := 0
	for i, alt := range v.Alt() {
		if alt != "." && alt != "*" && !(len(alt) > 0 && alt[0] == '<') {
			real = i + 1
			count++
		}
	}
	if count != 1 {
		return -1
	}
	return real
}

func sampleHasOnlyRefOrAlt(s *vcfgo.SampleGenotype, targetAlt int) bool {
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

func isHomozygous(gt []int) bool {
	if len(gt) == 0 || gt[0] < 0 {
		return false
	}
	for _, allele := range gt[1:] {
		if allele != gt[0] {
			return false
		}
	}
	return true
}

func orientedDepth(s *vcfgo.SampleGenotype, realAlt, trackAllele int) (sampleDepth, bool) {
	ref, err := s.RefDepth()
	if err != nil {
		return sampleDepth{}, false
	}
	alts, err := s.AltDepths()
	if err != nil || len(alts) < realAlt {
		return sampleDepth{}, false
	}
	alt := alts[realAlt-1]
	if trackAllele == 0 {
		return sampleDepth{Track: ref, Other: alt, DP: ref + alt}, true
	}
	return sampleDepth{Track: alt, Other: ref, DP: ref + alt}, true
}

func passesHardFilter(v *vcfgo.Variant, hf utils.HardFilterConfig) bool {
	isSNP, isIndel := classifyVariant(v)
	switch {
	case isSNP:
		return float64(v.Quality) >= hf.SNP_QUAL_Min &&
			infoMin(v, "QD", hf.SNP_QD_Min) && infoMax(v, "SOR", hf.SNP_SOR_Max) &&
			infoMax(v, "FS", hf.SNP_FS_Max) && infoMin(v, "MQ", hf.SNP_MQ_Min) &&
			infoMin(v, "MQRankSum", hf.SNP_MQRankSum_Min) && infoMin(v, "ReadPosRankSum", hf.SNP_ReadPosRankSum_Min)
	case isIndel:
		return float64(v.Quality) >= hf.INDEL_QUAL_Min &&
			infoMin(v, "QD", hf.INDEL_QD_Min) && infoMax(v, "FS", hf.INDEL_FS_Max) &&
			infoMax(v, "SOR", hf.INDEL_SOR_Max) && infoMin(v, "ReadPosRankSum", hf.INDEL_ReadPosRankSum_Min)
	default:
		return false
	}
}

func classifyVariant(v *vcfgo.Variant) (bool, bool) {
	refLen := len(v.Ref())
	isSNP := refLen == 1
	isIndel := false
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
	return isSNP, isIndel
}

func infoFloat(v *vcfgo.Variant, key string) (float64, bool) {
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

func infoMin(v *vcfgo.Variant, key string, min float64) bool {
	val, ok := infoFloat(v, key)
	return !ok || val >= min
}

func infoMax(v *vcfgo.Variant, key string, max float64) bool {
	val, ok := infoFloat(v, key)
	return !ok || val <= max
}

func calcTwoBulkStats(m *marker) {
	m.HighSI = ratio(m.HighDepth.Track, m.HighDepth.DP)
	m.LowSI = ratio(m.LowDepth.Track, m.LowDepth.DP)
	m.Delta = m.HighSI - m.LowSI
	m.G = gStatistic(m.HighDepth.Track, m.HighDepth.Other, m.LowDepth.Track, m.LowDepth.Other)
	m.ED4 = math.Pow(m.Delta, 4)
	m.LOD = lod(m.HighDepth.Other, m.HighDepth.Track, m.LowDepth.Other, m.LowDepth.Track)
	m.BF = betaBinomialBF(m.HighDepth.Track, m.HighDepth.Other, m.LowDepth.Track, m.LowDepth.Other)
	m.Depth = min(m.HighDepth.DP, m.LowDepth.DP)
	m.BRM = m.Delta
}

func calcOneBulkStats(m *marker) {
	m.HighDepth = m.OneDepth
	m.HighSI = ratio(m.OneDepth.Track, m.OneDepth.DP)
	m.LowSI = math.NaN()
	m.Delta = m.HighSI - 0.5
	m.G = gOneBulk(m.OneDepth.Track, m.OneDepth.Other)
	m.ED4 = math.Pow(m.Delta, 4)
	m.LOD = lodOneBulk(m.OneDepth.Track, m.OneDepth.Other)
	m.BF = betaBinomialOneBF(m.OneDepth.Track, m.OneDepth.Other)
	m.Depth = m.OneDepth.DP
	m.BRM = m.Delta
}

func ratio(a, b int) float64 {
	if b <= 0 {
		return math.NaN()
	}
	return float64(a) / float64(b)
}

func gStatistic(a1, r1, a2, r2 int) float64 {
	obs := []float64{float64(a1) + 0.5, float64(r1) + 0.5, float64(a2) + 0.5, float64(r2) + 0.5}
	row1, row2 := obs[0]+obs[1], obs[2]+obs[3]
	col1, col2 := obs[0]+obs[2], obs[1]+obs[3]
	total := row1 + row2
	exp := []float64{row1 * col1 / total, row1 * col2 / total, row2 * col1 / total, row2 * col2 / total}
	g := 0.0
	for i := range obs {
		g += obs[i] * math.Log(obs[i]/exp[i])
	}
	return 2 * g
}

func gOneBulk(alt, ref int) float64 {
	a, r := float64(alt)+0.5, float64(ref)+0.5
	e := (a + r) / 2
	return 2 * (a*math.Log(a/e) + r*math.Log(r/e))
}

func lod(ref1, alt1, ref2, alt2 int) float64 {
	n1, n2 := float64(ref1+alt1), float64(ref2+alt2)
	if n1 == 0 || n2 == 0 {
		return 0
	}
	p1 := clamp(float64(alt1) / n1)
	p2 := clamp(float64(alt2) / n2)
	p0 := clamp(float64(alt1+alt2) / (n1 + n2))
	ll := func(k, n, p float64) float64 { return k*math.Log(p) + (n-k)*math.Log(1-p) }
	return (ll(float64(alt1), n1, p1) + ll(float64(alt2), n2, p2) - ll(float64(alt1), n1, p0) - ll(float64(alt2), n2, p0)) / math.Ln10
}

func lodOneBulk(alt, ref int) float64 {
	n := float64(alt + ref)
	if n == 0 {
		return 0
	}
	p := clamp(float64(alt) / n)
	return n * (p*math.Log(p/0.5) + (1-p)*math.Log((1-p)/0.5)) / math.Ln10
}

func clamp(p float64) float64 {
	if p < 1e-10 {
		return 1e-10
	}
	if p > 1-1e-10 {
		return 1 - 1e-10
	}
	return p
}

func betaBinomialBF(a1, r1, a2, r2 int) float64 {
	sep := logBeta(float64(a1)+0.5, float64(r1)+0.5) - logBeta(0.5, 0.5)
	sep += logBeta(float64(a2)+0.5, float64(r2)+0.5) - logBeta(0.5, 0.5)
	shared := logBeta(float64(a1+a2)+0.5, float64(r1+r2)+0.5) - logBeta(0.5, 0.5)
	return (sep - shared) / math.Ln10
}

func betaBinomialOneBF(alt, ref int) float64 {
	total := float64(alt + ref)
	if total == 0 {
		return 0
	}
	flexible := logBeta(float64(alt)+1, float64(ref)+1) - logBeta(1, 1)
	fixed := total * math.Log(0.5)
	return (flexible - fixed) / math.Ln10
}

func logBeta(a, b float64) float64 {
	la, _ := math.Lgamma(a)
	lb, _ := math.Lgamma(b)
	lab, _ := math.Lgamma(a + b)
	return la + lb - lab
}

func sortMarkers(markers []marker) {
	sort.Slice(markers, func(i, j int) bool {
		if markers[i].Chrom == markers[j].Chrom {
			return markers[i].Pos < markers[j].Pos
		}
		return markers[i].Chrom < markers[j].Chrom
	})
}

func splitByChrom(markers []marker) map[string][]int {
	out := map[string][]int{}
	for i, m := range markers {
		out[m.Chrom] = append(out[m.Chrom], i)
	}
	return out
}

func assignBRM(markers []marker, window int) {
	byChrom := splitByChrom(markers)
	for _, idxs := range byChrom {
		for _, idx := range idxs {
			markers[idx].BRM = localLinear(markers, idxs, markers[idx].Pos, float64(window)/2, func(m marker) float64 { return m.Delta })
		}
	}
}

func smoothMarkers(markers []marker, window, step int) []smoothPoint {
	byChrom := splitByChrom(markers)
	bandwidth := float64(window) / 2
	var out []smoothPoint
	for chrom, idxs := range byChrom {
		minPos := markers[idxs[0]].Pos
		maxPos := markers[idxs[len(idxs)-1]].Pos
		for center := minPos; center <= maxPos; center += int64(step) {
			stats := map[string]float64{}
			n, meanH, meanL := 0, 0, 0
			for _, stat := range statNames {
				if stat == "brm" {
					stats[stat] = localLinear(markers, idxs, center, bandwidth, getter(stat))
				} else {
					v, count, hdp, ldp := weightedMean(markers, idxs, center, bandwidth, getter(stat))
					stats[stat] = v
					if count > n {
						n, meanH, meanL = count, hdp, ldp
					}
				}
			}
			if n >= minSNPsPerWindow {
				start := maxInt64(1, center-int64(window)/2)
				out = append(out, smoothPoint{Chrom: chrom, Start: start, End: center + int64(window)/2, Center: center, MarkerN: n, MeanHighDP: meanH, MeanLowDP: meanL, Stats: stats})
			}
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Chrom == out[j].Chrom {
			return out[i].Center < out[j].Center
		}
		return out[i].Chrom < out[j].Chrom
	})
	return out
}

func weightedMean(markers []marker, idxs []int, center int64, bandwidth float64, get func(marker) float64) (float64, int, int, int) {
	sumW, sumY, sumH, sumL := 0.0, 0.0, 0.0, 0.0
	n := 0
	for _, idx := range idxs {
		m := markers[idx]
		d := math.Abs(float64(m.Pos - center))
		w := tricubeWeight(d, bandwidth) * math.Sqrt(float64(max(1, m.Depth)))
		if w <= 0 {
			continue
		}
		y := get(m)
		if math.IsNaN(y) || math.IsInf(y, 0) {
			continue
		}
		sumW += w
		sumY += w * y
		sumH += float64(m.HighDepth.DP)
		sumL += float64(m.LowDepth.DP)
		n++
	}
	if sumW == 0 {
		return math.NaN(), n, 0, 0
	}
	return sumY / sumW, n, int(sumH / float64(n)), int(sumL / float64(max(1, n)))
}

func localLinear(markers []marker, idxs []int, center int64, bandwidth float64, get func(marker) float64) float64 {
	var sw, sx, sy, sxx, sxy float64
	for _, idx := range idxs {
		m := markers[idx]
		x := float64(m.Pos - center)
		w := tricubeWeight(math.Abs(x), bandwidth) * math.Sqrt(float64(max(1, m.Depth)))
		if w <= 0 {
			continue
		}
		y := get(m)
		if math.IsNaN(y) || math.IsInf(y, 0) {
			continue
		}
		sw += w
		sx += w * x
		sy += w * y
		sxx += w * x * x
		sxy += w * x * y
	}
	den := sw*sxx - sx*sx
	if sw == 0 {
		return math.NaN()
	}
	if math.Abs(den) < eps {
		return sy / sw
	}
	return (sxx*sy - sx*sxy) / den
}

func tricubeWeight(d, D float64) float64 {
	if D <= 0 || d >= D {
		return 0
	}
	x := 1 - math.Pow(d/D, 3)
	return x * x * x
}

func getter(name string) func(marker) float64 {
	return func(m marker) float64 {
		switch name {
		case "high_si":
			return m.HighSI
		case "low_si":
			return m.LowSI
		case "delta_si":
			return m.Delta
		case "g_statistic":
			return m.G
		case "ed4":
			return m.ED4
		case "lod":
			return m.LOD
		case "beta_binomial_bf":
			return m.BF
		case "brm":
			return m.BRM
		default:
			return math.NaN()
		}
	}
}

var thresholdCache sync.Map

func attachAdaptiveThresholds(points []smoothPoint, cfg utils.AnalysisConfig, m mode) {
	highAF := simulateAF(cfg.Population, cfg.HighBulkSize, cfg.Rep)
	lowAF := simulateAF(cfg.Population, cfg.LowBulkSize, cfg.Rep)
	progressEvery := max(1, len(points)/10)
	for i := range points {
		points[i].Thresholds = calcThresholdsCached(points[i].MeanHighDP, points[i].MeanLowDP, highAF, lowAF, cfg, m)
		if (i+1)%progressEvery == 0 || i+1 == len(points) {
			color.White("  Threshold windows completed: %d/%d", i+1, len(points))
		}
	}
}

func calcThresholdsCached(highDP, lowDP int, highAF, lowAF float64, cfg utils.AnalysisConfig, m mode) map[string]thresholdPair {
	key := fmt.Sprintf("%d_%d_%.4f_%.4f_%d_%d", highDP, lowDP, highAF, lowAF, cfg.Rep, m)
	if val, ok := thresholdCache.Load(key); ok {
		return val.(map[string]thresholdPair)
	}
	t := calcThresholds(highDP, lowDP, highAF, lowAF, cfg, m)
	actual, _ := thresholdCache.LoadOrStore(key, t)
	return actual.(map[string]thresholdPair)
}

func calcThresholds(highDP, lowDP int, highAF, lowAF float64, cfg utils.AnalysisConfig, m mode) map[string]thresholdPair {
	out := map[string]thresholdPair{}
	vals := map[string][]float64{}
	for _, stat := range statNames {
		vals[stat] = []float64{}
	}
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	if m == modeTwoParentsOneBulk || m == modeOneParentOneBulk || lowDP <= 0 {
		n := max(1, highDP)
		for i := 0; i < cfg.Rep; i++ {
			alt := binomial(rng, n, 0.5)
			mk := marker{OneDepth: sampleDepth{Track: alt, Other: n - alt, DP: n}}
			calcOneBulkStats(&mk)
			for _, stat := range statNames {
				vals[stat] = append(vals[stat], getter(stat)(mk))
			}
		}
	} else {
		highDP, lowDP = max(1, highDP), max(1, lowDP)
		for i := 0; i < cfg.Rep; i++ {
			ha := binomial(rng, highDP, highAF)
			la := binomial(rng, lowDP, lowAF)
			mk := marker{HighDepth: sampleDepth{Track: ha, Other: highDP - ha, DP: highDP}, LowDepth: sampleDepth{Track: la, Other: lowDP - la, DP: lowDP}}
			calcTwoBulkStats(&mk)
			for _, stat := range statNames {
				vals[stat] = append(vals[stat], getter(stat)(mk))
			}
		}
	}
	for _, stat := range statNames {
		sort.Float64s(vals[stat])
		t := thresholdPair{Upper: map[float64]float64{}, Lower: map[float64]float64{}}
		for _, alpha := range cfg.Alphas {
			t.Upper[alpha] = percentile(vals[stat], 1-alpha)
			t.Lower[alpha] = percentile(vals[stat], alpha)
		}
		out[stat] = t
	}
	return out
}

func simulateAF(pop string, bulkSize int, rep int) float64 {
	if bulkSize <= 0 {
		return 0.5
	}
	probs := []float64{0.25, 0.5, 0.25}
	switch strings.ToUpper(pop) {
	case "RIL":
		probs = []float64{0.5, 0, 0.5}
	case "BC":
		probs = []float64{0.5, 0.5, 0}
	}
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	total := 0.0
	for i := 0; i < max(1, rep); i++ {
		sum := 0.0
		for j := 0; j < bulkSize; j++ {
			x := rng.Float64()
			switch {
			case x < probs[0]:
				sum += 0
			case x < probs[0]+probs[1]:
				sum += 0.5
			default:
				sum += 1
			}
		}
		total += sum / float64(bulkSize)
	}
	return total / float64(max(1, rep))
}

func binomial(rng *rand.Rand, n int, p float64) int {
	x := 0
	for i := 0; i < n; i++ {
		if rng.Float64() < p {
			x++
		}
	}
	return x
}

func attachRobustZ(points []smoothPoint) {
	for _, stat := range statNames {
		values := make([]float64, 0, len(points))
		for _, p := range points {
			v := p.Stats[stat]
			if !math.IsNaN(v) && !math.IsInf(v, 0) {
				values = append(values, v)
			}
		}
		med, mad := robustBackground(values, 0.01)
		scale := 1.4826 * mad
		if scale == 0 {
			scale = 1
		}
		for i := range points {
			if points[i].Z == nil {
				points[i].Z = map[string]float64{}
			}
			points[i].Z[stat] = (points[i].Stats[stat] - med) / scale
		}
	}
}

func robustBackground(values []float64, trim float64) (float64, float64) {
	if len(values) == 0 {
		return 0, 1
	}
	cp := append([]float64(nil), values...)
	sort.Float64s(cp)
	cut := int(math.Round(float64(len(cp)) * (1 - trim)))
	if cut < 1 {
		cut = 1
	}
	cp = cp[:cut]
	med := percentile(cp, 0.5)
	dev := make([]float64, len(cp))
	for i, v := range cp {
		dev[i] = math.Abs(v - med)
	}
	sort.Float64s(dev)
	return med, percentile(dev, 0.5)
}

func detectPermutationQTLs(points []smoothPoint, cfg utils.AnalysisConfig) []qtl {
	var out []qtl
	for _, stat := range statNames {
		if stat == "high_si" || stat == "low_si" {
			continue
		}
		for _, alpha := range cfg.Alphas {
			label := fmt.Sprintf("p%.2g", alpha)
			out = append(out, detectIntervals(points, stat, "Permutation", label, cfg, func(p smoothPoint) (float64, bool) {
				v := p.Stats[stat]
				t := p.Thresholds[stat]
				if signedStat(stat) {
					return math.Abs(v), v >= t.Upper[alpha] || v <= t.Lower[alpha]
				}
				return v, v >= t.Upper[alpha]
			})...)
		}
	}
	return out
}

func detectZQTLs(points []smoothPoint, cfg utils.AnalysisConfig) []qtl {
	var out []qtl
	for _, stat := range statNames {
		if stat == "high_si" || stat == "low_si" {
			continue
		}
		for _, z := range []float64{2, 3} {
			label := fmt.Sprintf("z%.0f", z)
			out = append(out, detectIntervals(points, stat, "ZScore", label, cfg, func(p smoothPoint) (float64, bool) {
				v := math.Abs(p.Z[stat])
				return v, v >= z
			})...)
		}
	}
	return out
}

func detectConsensusQTLs(points []smoothPoint, cfg utils.AnalysisConfig) []qtl {
	alpha := 0.01
	if len(cfg.Alphas) > 0 {
		alpha = cfg.Alphas[len(cfg.Alphas)-1]
	}
	return detectIntervals(points, "consensus", "Consensus", fmt.Sprintf("p%.2g", alpha), cfg, func(p smoothPoint) (float64, bool) {
		count := 0
		for _, stat := range []string{"delta_si", "g_statistic", "ed4", "lod", "beta_binomial_bf", "brm"} {
			v := p.Stats[stat]
			t := p.Thresholds[stat]
			if signedStat(stat) {
				if v >= t.Upper[alpha] || v <= t.Lower[alpha] {
					count++
				}
			} else if v >= t.Upper[alpha] {
				count++
			}
		}
		return float64(count), count >= consensusMin
	})
}

func detectMaxZQTLs(points []smoothPoint, cfg utils.AnalysisConfig) []qtl {
	return detectIntervals(points, "max_abs_z", "MaxZ", "z3", cfg, func(p smoothPoint) (float64, bool) {
		maxZ := 0.0
		for _, stat := range []string{"delta_si", "g_statistic", "ed4", "lod", "beta_binomial_bf", "brm"} {
			maxZ = math.Max(maxZ, math.Abs(p.Z[stat]))
		}
		return maxZ, maxZ >= 3
	})
}

func detectIntervals(points []smoothPoint, stat, method, threshold string, cfg utils.AnalysisConfig, pass func(smoothPoint) (float64, bool)) []qtl {
	var out []qtl
	var cur *qtl
	for _, p := range points {
		v, ok := pass(p)
		if !ok {
			if cur != nil {
				out = appendQTL(out, *cur, cfg.MinQTLWidth)
				cur = nil
			}
			continue
		}
		if cur == nil || cur.Chrom != p.Chrom || p.Start-cur.End > cfg.MergeDistance {
			if cur != nil {
				out = appendQTL(out, *cur, cfg.MinQTLWidth)
			}
			cur = &qtl{Chrom: p.Chrom, Start: p.Start, End: p.End, Stat: stat, Method: method, Threshold: threshold, PeakPos: p.Center, PeakValue: v, Markers: p.MarkerN}
			continue
		}
		cur.End = maxInt64(cur.End, p.End)
		cur.Markers += p.MarkerN
		if v > cur.PeakValue {
			cur.PeakValue = v
			cur.PeakPos = p.Center
		}
	}
	if cur != nil {
		out = appendQTL(out, *cur, cfg.MinQTLWidth)
	}
	return out
}

func appendQTL(out []qtl, q qtl, minWidth int64) []qtl {
	if q.End-q.Start+1 >= minWidth {
		return append(out, q)
	}
	return out
}

func detectBRMBlocks(points []smoothPoint, cfg utils.AnalysisConfig) []brmBlock {
	popLevel := 1
	if strings.ToUpper(cfg.Population) == "RIL" {
		popLevel = 0
	}
	u := distuv.UnitNormal.Quantile(1 - defaultBRMAlpha/2)
	n1, n2 := float64(max(1, cfg.HighBulkSize)), float64(max(1, cfg.LowBulkSize))
	scale := (n1 + n2) / (math.Pow(2, float64(popLevel)) * n1 * n2)
	return detectBRMIntervals(points, func(p smoothPoint) (float64, float64, bool) {
		afp := math.Max(afpFloor, math.Min(1-afpFloor, (p.Stats["high_si"]+p.Stats["low_si"])/2))
		th := u * math.Sqrt(scale*afp*(1-afp))
		v := math.Abs(p.Stats["brm"])
		return v, th, v >= th
	})
}

func detectBRMIntervals(points []smoothPoint, pass func(smoothPoint) (float64, float64, bool)) []brmBlock {
	var out []brmBlock
	var cur *brmBlock
	for _, p := range points {
		v, th, ok := pass(p)
		if !ok {
			if cur != nil {
				out = append(out, *cur)
				cur = nil
			}
			continue
		}
		if cur == nil || cur.Chrom != p.Chrom {
			if cur != nil {
				out = append(out, *cur)
			}
			cur = &brmBlock{Chrom: p.Chrom, Start: p.Start, End: p.End, PeakPos: p.Center, Peak: v, Threshold: th}
			continue
		}
		cur.End = p.End
		if v > cur.Peak {
			cur.Peak, cur.PeakPos, cur.Threshold = v, p.Center, th
		}
	}
	if cur != nil {
		out = append(out, *cur)
	}
	return out
}

func markBRMHits(points []smoothPoint, blocks []brmBlock) {
	for i := range points {
		for _, b := range blocks {
			if points[i].Chrom == b.Chrom && points[i].Start <= b.End && points[i].End >= b.Start {
				points[i].BRMHit = true
				break
			}
		}
	}
}

func intersectWithBRM(qtls []qtl, blocks []brmBlock, method string) []qtl {
	var out []qtl
	for _, q := range qtls {
		for _, b := range blocks {
			if q.Chrom == b.Chrom && q.Start <= b.End && q.End >= b.Start {
				q.Method = method
				out = append(out, q)
				break
			}
		}
	}
	return out
}

func writeMarkers(path string, markers []marker) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()
	_ = w.Write([]string{"chrom", "pos", "ref", "alt", "high_track_dp", "high_other_dp", "low_track_dp", "low_other_dp", "high_si", "low_si", "delta_si", "g_statistic", "ed4", "lod", "beta_binomial_bf", "brm"})
	for _, m := range markers {
		_ = w.Write([]string{m.Chrom, i64(m.Pos), m.Ref, m.Alt, itoa(m.HighDepth.Track), itoa(m.HighDepth.Other), itoa(m.LowDepth.Track), itoa(m.LowDepth.Other), ftoa(m.HighSI), ftoa(m.LowSI), ftoa(m.Delta), ftoa(m.G), ftoa(m.ED4), ftoa(m.LOD), ftoa(m.BF), ftoa(m.BRM)})
	}
	return w.Error()
}

func writeSmoothed(path string, points []smoothPoint, alphas []float64) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()
	header := []string{"chrom", "start", "end", "center", "marker_count", "mean_high_dp", "mean_low_dp", "brm_block"}
	for _, stat := range statNames {
		header = append(header, stat, stat+"_z")
		for _, a := range alphas {
			header = append(header, stat+"_upper_p"+ftoa(a), stat+"_lower_p"+ftoa(a))
		}
	}
	_ = w.Write(header)
	for _, p := range points {
		row := []string{p.Chrom, i64(p.Start), i64(p.End), i64(p.Center), itoa(p.MarkerN), itoa(p.MeanHighDP), itoa(p.MeanLowDP), strconv.FormatBool(p.BRMHit)}
		for _, stat := range statNames {
			row = append(row, ftoa(p.Stats[stat]), ftoa(p.Z[stat]))
			for _, a := range alphas {
				row = append(row, ftoa(p.Thresholds[stat].Upper[a]), ftoa(p.Thresholds[stat].Lower[a]))
			}
		}
		_ = w.Write(row)
	}
	return w.Error()
}

func writeBRM(path string, blocks []brmBlock) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()
	_ = w.Write([]string{"chrom", "start", "end", "peak_pos", "peak_abs_brm", "threshold"})
	for _, b := range blocks {
		_ = w.Write([]string{b.Chrom, i64(b.Start), i64(b.End), i64(b.PeakPos), ftoa(b.Peak), ftoa(b.Threshold)})
	}
	return w.Error()
}

func writeQTLs(path string, qtls []qtl) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()
	_ = w.Write([]string{"chrom", "start", "end", "stat", "method", "threshold", "peak_pos", "peak_value", "supporting_markers"})
	for _, q := range qtls {
		_ = w.Write([]string{q.Chrom, i64(q.Start), i64(q.End), q.Stat, q.Method, q.Threshold, i64(q.PeakPos), ftoa(q.PeakValue), itoa(q.Markers)})
	}
	return w.Error()
}

func plotAll(outDir string, points []smoothPoint, blocks []brmBlock, alphas []float64) error {
	if err := os.MkdirAll(outDir, 0o755); err != nil {
		return err
	}
	page := components.NewPage()
	page.PageTitle = "GoBSAseq Best Of Both"
	page.SetLayout(components.PageFlexLayout)
	for _, stat := range statNames {
		page.AddCharts(plotStat(stat, points, alphas, false))
		page.AddCharts(plotStat(stat+"_robust_z", points, alphas, true))
	}
	page.AddCharts(plotComposite(points))
	f, err := os.Create(filepath.Join(outDir, "GoBSAseq_best_of_both.html"))
	if err != nil {
		return err
	}
	defer f.Close()
	_ = blocks
	return page.Render(f)
}

func plotStat(title string, points []smoothPoint, alphas []float64, z bool) *charts.Line {
	stat := strings.TrimSuffix(title, "_robust_z")
	line := charts.NewLine()
	line.SetGlobalOptions(charts.WithTitleOpts(opts.Title{Title: title}), charts.WithXAxisOpts(opts.XAxis{Name: "genome window", AxisLabel: &opts.AxisLabel{Show: opts.Bool(false)}}))
	var x []string
	var y []opts.LineData
	for _, p := range points {
		x = append(x, fmt.Sprintf("%s:%d", p.Chrom, p.Center))
		v := p.Stats[stat]
		if z {
			v = p.Z[stat]
		}
		y = append(y, opts.LineData{Value: clean(v)})
	}
	line.SetXAxis(x).AddSeries(stat, y)
	if z {
		addFlat(line, x, "z=2", 2)
		addFlat(line, x, "z=-2", -2)
		addFlat(line, x, "z=3", 3)
		addFlat(line, x, "z=-3", -3)
	} else if len(points) > 0 {
		for _, a := range alphas {
			addFlat(line, x, "upper p"+ftoa(a), points[0].Thresholds[stat].Upper[a])
			if signedStat(stat) {
				addFlat(line, x, "lower p"+ftoa(a), points[0].Thresholds[stat].Lower[a])
			}
		}
	}
	return line
}

func plotComposite(points []smoothPoint) *charts.Line {
	line := charts.NewLine()
	line.SetGlobalOptions(charts.WithTitleOpts(opts.Title{Title: "max absolute robust z"}), charts.WithXAxisOpts(opts.XAxis{AxisLabel: &opts.AxisLabel{Show: opts.Bool(false)}}))
	var x []string
	var y []opts.LineData
	for _, p := range points {
		maxZ := 0.0
		for _, stat := range statNames {
			maxZ = math.Max(maxZ, math.Abs(p.Z[stat]))
		}
		x = append(x, fmt.Sprintf("%s:%d", p.Chrom, p.Center))
		y = append(y, opts.LineData{Value: maxZ})
	}
	line.SetXAxis(x).AddSeries("max_abs_z", y)
	addFlat(line, x, "z=3", 3)
	return line
}

func addFlat(line *charts.Line, x []string, name string, value float64) {
	data := make([]opts.LineData, len(x))
	for i := range data {
		data[i] = opts.LineData{Value: clean(value)}
	}
	line.AddSeries(name, data)
}

func clean(v float64) interface{} {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return nil
	}
	return v
}

func percentile(vals []float64, p float64) float64 {
	if len(vals) == 0 {
		return math.NaN()
	}
	if p <= 0 {
		return vals[0]
	}
	if p >= 1 {
		return vals[len(vals)-1]
	}
	pos := p * float64(len(vals)-1)
	lo, hi := int(math.Floor(pos)), int(math.Ceil(pos))
	if lo == hi {
		return vals[lo]
	}
	return vals[lo] + (vals[hi]-vals[lo])*(pos-float64(lo))
}

func signedStat(stat string) bool {
	return stat == "delta_si" || stat == "brm"
}

func modeName(m mode) string {
	switch m {
	case modeTwoParentsTwoBulks:
		return "two parents + two bulks"
	case modeBulksOnly:
		return "bulks only"
	case modeTwoParentsOneBulk:
		return "two parents + one bulk"
	default:
		return "one parent + one bulk"
	}
}

func logStage(msg string) {
	color.Cyan("==> %s", msg)
}

func formatAlphas(alphas []float64) string {
	parts := make([]string, 0, len(alphas))
	for _, a := range alphas {
		parts = append(parts, ftoa(a))
	}
	return strings.Join(parts, ",")
}

func ftoa(v float64) string {
	if math.IsNaN(v) {
		return "NA"
	}
	return strconv.FormatFloat(v, 'g', 8, 64)
}

func i64(v int64) string { return strconv.FormatInt(v, 10) }
func itoa(v int) string  { return strconv.Itoa(v) }

func maxInt64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}
