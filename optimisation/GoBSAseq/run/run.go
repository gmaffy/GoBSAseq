package run

import (
	"bufio"
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
	"time"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/opts"
	"gonum.org/v1/gonum/stat/distuv"
)

const eps = 1e-12

var statNames = []string{"snp_index_high", "snp_index_low", "delta_snp_index", "g_statistic", "euclidean_distance", "lod", "beta_binomial_bf", "brm"}

type sampleDepth struct {
	Ref int
	Alt int
	DP  int
}

type marker struct {
	Chrom string
	Pos   int64
	Ref   string
	Alt   string

	High sampleDepth
	Low  sampleDepth
	One  sampleDepth

	HighIndex float64
	LowIndex  float64
	Delta     float64
	G         float64
	ED        float64
	LOD       float64
	BF        float64
	BRM       float64
}

type smoothPoint struct {
	Chrom      string
	Start      int64
	End        int64
	Center     int64
	MarkerN    int
	Stats      map[string]float64
	Thresholds map[string]map[float64]float64
	Z          map[string]float64
}

type qtl struct {
	Chrom     string
	Start     int64
	End       int64
	Stat      string
	Threshold string
	PeakPos   int64
	PeakValue float64
	Markers   int
}

type analysisMode int

const (
	modeTwoParentsTwoBulks analysisMode = iota
	modeBulksOnly
	modeTwoParentsOneBulk
	modeOneParentOneBulk
)

func Run(cfg utils.AnalysisConfig, hf utils.HardFilterConfig) error {
	start := time.Now()
	if cfg.VCF == "" {
		return errors.New("variant file is required")
	}
	if cfg.WindowSize <= 0 || cfg.StepSize <= 0 {
		return errors.New("window-size and step-size must be positive")
	}
	if cfg.Rep < 1 {
		cfg.Rep = 1
	}
	if len(cfg.Alphas) == 0 {
		cfg.Alphas = []float64{0.05, 0.01}
	}
	if cfg.OutputDir == "" {
		cfg.OutputDir = "."
	}
	if err := os.MkdirAll(cfg.OutputDir, 0o755); err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Join(cfg.OutputDir, "plots"), 0o755); err != nil {
		return err
	}

	mode, err := inferMode(cfg)
	if err != nil {
		return err
	}
	logStage("Starting GoBSAseq")
	color.White("  VCF: %s", cfg.VCF)
	color.White("  Output directory: %s", cfg.OutputDir)
	color.White("  Mode: %s", modeName(mode))
	color.White("  Window/step: %d/%d bp", cfg.WindowSize, cfg.StepSize)
	color.White("  Simulations: %d reps; alpha levels: %s", cfg.Rep, formatAlphas(cfg.Alphas))

	logStage("Reading and filtering variants")
	markers, err := readMarkers(cfg, hf, mode)
	if err != nil {
		return err
	}
	if len(markers) == 0 {
		return errors.New("no informative markers passed filters")
	}
	color.White("  Informative markers retained: %d across %d chromosomes", len(markers), countChromosomes(markers))

	logStage("Calculating marker-level statistics")
	sortMarkers(markers)
	assignBRM(markers, cfg.WindowSize)

	logStage("Smoothing statistics with Gaussian kernel weighting")
	points := smoothMarkers(markers, cfg.WindowSize, cfg.StepSize)
	color.White("  Smoothed windows generated: %d", len(points))

	logStage("Simulating null thresholds")
	thresholds := simulateThresholds(markers, cfg, mode)
	attachThresholds(points, thresholds)

	logStage("Normalising smoothed statistics with robust z-scores")
	attachRobustZ(points)

	logStage("Detecting QTL intervals")
	statQTLs := detectThresholdQTLs(points, cfg.MinQTLWidth, cfg.MergeDistance)
	zQTLs := detectZQTLs(points, cfg.MinQTLWidth, cfg.MergeDistance)
	color.White("  Simulated-threshold QTLs: %d", len(statQTLs))
	color.White("  Robust-z final QTLs: %d", len(zQTLs))

	logStage("Writing result tables")
	if err := writeMarkers(filepath.Join(cfg.OutputDir, "variant_statistics.tsv"), markers); err != nil {
		return err
	}
	color.White("  Wrote %s", filepath.Join(cfg.OutputDir, "variant_statistics.tsv"))
	if err := writeSmooth(filepath.Join(cfg.OutputDir, "smoothed_statistics.tsv"), points); err != nil {
		return err
	}
	color.White("  Wrote %s", filepath.Join(cfg.OutputDir, "smoothed_statistics.tsv"))
	if err := writeThresholds(filepath.Join(cfg.OutputDir, "thresholds.tsv"), thresholds); err != nil {
		return err
	}
	color.White("  Wrote %s", filepath.Join(cfg.OutputDir, "thresholds.tsv"))
	if err := writeQTLs(filepath.Join(cfg.OutputDir, "qtls_by_simulated_thresholds.tsv"), statQTLs); err != nil {
		return err
	}
	color.White("  Wrote %s", filepath.Join(cfg.OutputDir, "qtls_by_simulated_thresholds.tsv"))
	if err := writeQTLs(filepath.Join(cfg.OutputDir, "final_qtls_by_robust_z.tsv"), zQTLs); err != nil {
		return err
	}
	color.White("  Wrote %s", filepath.Join(cfg.OutputDir, "final_qtls_by_robust_z.tsv"))

	logStage("Rendering go-echarts plots")
	if err := plotAll(cfg.OutputDir, points); err != nil {
		return err
	}
	color.White("  Wrote HTML plots to %s", filepath.Join(cfg.OutputDir, "plots"))

	color.Green("GoBSAseq complete in %s: %d informative markers, %d smoothed windows, %d simulated-threshold QTLs, %d robust-z QTLs", time.Since(start).Round(time.Millisecond), len(markers), len(points), len(statQTLs), len(zQTLs))
	return nil
}

func logStage(msg string) {
	color.Cyan("==> %s", msg)
}

func modeName(mode analysisMode) string {
	switch mode {
	case modeTwoParentsTwoBulks:
		return "two parents + two bulks"
	case modeBulksOnly:
		return "bulks only"
	case modeTwoParentsOneBulk:
		return "two parents + one bulk"
	case modeOneParentOneBulk:
		return "one parent + one bulk"
	default:
		return "unknown"
	}
}

func formatAlphas(alphas []float64) string {
	parts := make([]string, 0, len(alphas))
	for _, alpha := range alphas {
		parts = append(parts, ftoa(alpha))
	}
	return strings.Join(parts, ",")
}

func countChromosomes(markers []marker) int {
	seen := map[string]struct{}{}
	for _, m := range markers {
		seen[m.Chrom] = struct{}{}
	}
	return len(seen)
}

func inferMode(cfg utils.AnalysisConfig) (analysisMode, error) {
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
		return modeBulksOnly, fmt.Errorf("unsupported sample setup: provide two bulks, two parents plus two bulks, or parent(s) plus one bulk")
	}
}

func readMarkers(cfg utils.AnalysisConfig, hf utils.HardFilterConfig, mode analysisMode) ([]marker, error) {
	r, closeFn, err := openMaybeGzip(cfg.VCF)
	if err != nil {
		return nil, err
	}
	defer closeFn()

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 1024), 64*1024*1024)
	sampleIndex := map[string]int{}
	var markers []marker
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "##") {
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			parts := strings.Split(line, "\t")
			for i := 9; i < len(parts); i++ {
				sampleIndex[parts[i]] = i
			}
			continue
		}
		if strings.HasPrefix(line, "#") || strings.TrimSpace(line) == "" {
			continue
		}
		parts := strings.Split(line, "\t")
		if len(parts) < 10 {
			continue
		}
		if strings.Contains(parts[4], ",") {
			continue
		}
		isSNP := len(parts[3]) == 1 && len(parts[4]) == 1
		if !passesHardFilters(parts, hf, isSNP) {
			continue
		}
		pos, err := strconv.ParseInt(parts[1], 10, 64)
		if err != nil {
			continue
		}
		format := strings.Split(parts[8], ":")
		get := func(name string) (sampleDepth, string, bool) {
			idx, ok := sampleIndex[name]
			if !ok || idx >= len(parts) {
				return sampleDepth{}, "", false
			}
			return parseSample(format, parts[idx])
		}

		trackAllele := 1
		if mode == modeTwoParentsTwoBulks || mode == modeTwoParentsOneBulk {
			hi, hiGT, okH := get(cfg.HighParentName)
			lo, loGT, okL := get(cfg.LowParentName)
			if !okH || !okL || hi.DP < cfg.HighParentDepth || lo.DP < cfg.LowParentDepth {
				continue
			}
			hiAllele, okH := homozygousAllele(hiGT)
			loAllele, okL := homozygousAllele(loGT)
			if !okH || !okL || hiAllele == loAllele {
				continue
			}
			trackAllele = hiAllele
		}

		m := marker{Chrom: parts[0], Pos: pos, Ref: parts[3], Alt: parts[4]}
		if mode == modeTwoParentsTwoBulks || mode == modeBulksOnly {
			hi, _, okH := get(cfg.HighBulkName)
			lo, _, okL := get(cfg.LowBulkName)
			if !okH || !okL || hi.DP < cfg.HighBulkDepth || lo.DP < cfg.LowBulkDepth {
				continue
			}
			m.High = orientDepth(hi, trackAllele)
			m.Low = orientDepth(lo, trackAllele)
			calcTwoBulkStats(&m)
		} else {
			bulkName, minDepth := oneBulkName(cfg)
			one, _, ok := get(bulkName)
			if !ok || one.DP < minDepth {
				continue
			}
			m.One = orientDepth(one, trackAllele)
			calcOneBulkStats(&m)
		}
		markers = append(markers, m)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return markers, nil
}

func openMaybeGzip(path string) (io.Reader, func(), error) {
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

func parseSample(format []string, value string) (sampleDepth, string, bool) {
	fields := strings.Split(value, ":")
	data := map[string]string{}
	for i, k := range format {
		if i < len(fields) {
			data[k] = fields[i]
		}
	}
	gt := data["GT"]
	ad := strings.Split(data["AD"], ",")
	if len(ad) < 2 {
		return sampleDepth{}, gt, false
	}
	ref, errR := strconv.Atoi(ad[0])
	alt, errA := strconv.Atoi(ad[1])
	if errR != nil || errA != nil {
		return sampleDepth{}, gt, false
	}
	dp := ref + alt
	if rawDP, ok := data["DP"]; ok {
		if v, err := strconv.Atoi(rawDP); err == nil && v > 0 {
			dp = v
		}
	}
	return sampleDepth{Ref: ref, Alt: alt, DP: dp}, gt, true
}

func homozygousAllele(gt string) (int, bool) {
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) < 2 || parts[0] == "." || parts[1] == "." || parts[0] != parts[1] {
		return 0, false
	}
	allele, err := strconv.Atoi(parts[0])
	if err != nil || allele < 0 || allele > 1 {
		return 0, false
	}
	return allele, true
}

func orientDepth(d sampleDepth, trackAllele int) sampleDepth {
	if trackAllele == 0 {
		return sampleDepth{Ref: d.Alt, Alt: d.Ref, DP: d.DP}
	}
	return d
}

func passesHardFilters(parts []string, hf utils.HardFilterConfig, isSNP bool) bool {
	qual, _ := strconv.ParseFloat(parts[5], 64)
	info := parseInfo(parts[7])
	if isSNP {
		return qual >= hf.SNP_QUAL_Min &&
			infoMin(info, "QD", hf.SNP_QD_Min) &&
			infoMax(info, "SOR", hf.SNP_SOR_Max) &&
			infoMax(info, "FS", hf.SNP_FS_Max) &&
			infoMin(info, "MQ", hf.SNP_MQ_Min) &&
			infoMin(info, "MQRankSum", hf.SNP_MQRankSum_Min) &&
			infoMin(info, "ReadPosRankSum", hf.SNP_ReadPosRankSum_Min)
	}
	return qual >= hf.INDEL_QUAL_Min &&
		infoMin(info, "QD", hf.INDEL_QD_Min) &&
		infoMax(info, "FS", hf.INDEL_FS_Max) &&
		infoMax(info, "SOR", hf.INDEL_SOR_Max) &&
		infoMin(info, "ReadPosRankSum", hf.INDEL_ReadPosRankSum_Min)
}

func parseInfo(raw string) map[string]float64 {
	out := map[string]float64{}
	for _, item := range strings.Split(raw, ";") {
		kv := strings.SplitN(item, "=", 2)
		if len(kv) != 2 {
			continue
		}
		if v, err := strconv.ParseFloat(strings.Split(kv[1], ",")[0], 64); err == nil {
			out[kv[0]] = v
		}
	}
	return out
}

func infoMin(info map[string]float64, key string, min float64) bool {
	v, ok := info[key]
	return !ok || v >= min
}

func infoMax(info map[string]float64, key string, max float64) bool {
	v, ok := info[key]
	return !ok || v <= max
}

func oneBulkName(cfg utils.AnalysisConfig) (string, int) {
	switch {
	case cfg.OneBulkName != "":
		return cfg.OneBulkName, cfg.OneBulkDepth
	case cfg.HighBulkName != "":
		return cfg.HighBulkName, cfg.HighBulkDepth
	default:
		return cfg.LowBulkName, cfg.LowBulkDepth
	}
}

func calcTwoBulkStats(m *marker) {
	m.HighIndex = ratio(m.High.Alt, m.High.Ref+m.High.Alt)
	m.LowIndex = ratio(m.Low.Alt, m.Low.Ref+m.Low.Alt)
	m.Delta = m.HighIndex - m.LowIndex
	m.G = gStatistic(m.High.Alt, m.High.Ref, m.Low.Alt, m.Low.Ref)
	m.ED = math.Sqrt(math.Pow(m.HighIndex-m.LowIndex, 2) + math.Pow((1-m.HighIndex)-(1-m.LowIndex), 2))
	m.LOD = lodFromG(m.G, 1)
	m.BF = betaBinomialBF(m.High.Alt, m.High.Ref, m.Low.Alt, m.Low.Ref)
	m.BRM = m.Delta
}

func calcOneBulkStats(m *marker) {
	m.High = m.One
	m.HighIndex = ratio(m.One.Alt, m.One.Ref+m.One.Alt)
	m.LowIndex = math.NaN()
	m.Delta = m.HighIndex - 0.5
	m.G = gOneBulk(m.One.Alt, m.One.Ref, 0.5)
	m.ED = math.Sqrt(math.Pow(m.HighIndex-0.5, 2) + math.Pow((1-m.HighIndex)-0.5, 2))
	m.LOD = lodFromG(m.G, 1)
	m.BF = betaBinomialOneBF(m.One.Alt, m.One.Ref, 0.5)
	m.BRM = m.Delta
}

func ratio(a, b int) float64 {
	if b <= 0 {
		return math.NaN()
	}
	return float64(a) / float64(b)
}

func gStatistic(a1, r1, a2, r2 int) float64 {
	row1 := float64(a1 + r1)
	row2 := float64(a2 + r2)
	colA := float64(a1 + a2)
	colR := float64(r1 + r2)
	total := row1 + row2
	if total <= 0 {
		return 0
	}
	obs := []float64{float64(a1), float64(r1), float64(a2), float64(r2)}
	exp := []float64{row1 * colA / total, row1 * colR / total, row2 * colA / total, row2 * colR / total}
	g := 0.0
	for i := range obs {
		if obs[i] > 0 && exp[i] > 0 {
			g += 2 * obs[i] * math.Log(obs[i]/exp[i])
		}
	}
	return g
}

func gOneBulk(alt, ref int, p float64) float64 {
	n := float64(alt + ref)
	obs := []float64{float64(alt), float64(ref)}
	exp := []float64{n * p, n * (1 - p)}
	g := 0.0
	for i := range obs {
		if obs[i] > 0 && exp[i] > 0 {
			g += 2 * obs[i] * math.Log(obs[i]/exp[i])
		}
	}
	return g
}

func lodFromG(g float64, df float64) float64 {
	if g <= 0 {
		return 0
	}
	p := 1 - distuv.ChiSquared{K: df}.CDF(g)
	if p < eps {
		p = eps
	}
	return -math.Log10(p)
}

func betaBinomialBF(a1, r1, a2, r2 int) float64 {
	sep := logBetaBinomial(a1, r1, 1, 1) + logBetaBinomial(a2, r2, 1, 1)
	shared := logBetaBinomial(a1+a2, r1+r2, 1, 1)
	return (sep - shared) / math.Ln10
}

func betaBinomialOneBF(alt, ref int, p float64) float64 {
	flexible := logBetaBinomial(alt, ref, 1, 1)
	fixed := logBinomialPMF(alt, alt+ref, p)
	return (flexible - fixed) / math.Ln10
}

func logBetaBinomial(alt, ref int, alpha, beta float64) float64 {
	a := float64(alt)
	r := float64(ref)
	l1, _ := math.Lgamma(a + alpha)
	l2, _ := math.Lgamma(r + beta)
	l3, _ := math.Lgamma(a + r + alpha + beta)
	l4, _ := math.Lgamma(alpha + beta)
	l5, _ := math.Lgamma(alpha)
	l6, _ := math.Lgamma(beta)
	return l1 + l2 - l3 + l4 - l5 - l6
}

func logBinomialPMF(k, n int, p float64) float64 {
	lN, _ := math.Lgamma(float64(n + 1))
	lK, _ := math.Lgamma(float64(k + 1))
	lNK, _ := math.Lgamma(float64(n - k + 1))
	return lN - lK - lNK + float64(k)*math.Log(p) + float64(n-k)*math.Log(1-p)
}

func sortMarkers(markers []marker) {
	sort.Slice(markers, func(i, j int) bool {
		if markers[i].Chrom == markers[j].Chrom {
			return markers[i].Pos < markers[j].Pos
		}
		return markers[i].Chrom < markers[j].Chrom
	})
}

func assignBRM(markers []marker, window int) {
	byChrom := splitByChrom(markers)
	for _, idxs := range byChrom {
		for _, idx := range idxs {
			markers[idx].BRM = localLinear(markers, idxs, markers[idx].Pos, float64(window)/3, func(m marker) float64 { return m.Delta })
		}
	}
}

func smoothMarkers(markers []marker, window, step int) []smoothPoint {
	byChrom := splitByChrom(markers)
	var points []smoothPoint
	bandwidth := math.Max(1, float64(window)/3)
	for chrom, idxs := range byChrom {
		minPos := markers[idxs[0]].Pos
		maxPos := markers[idxs[len(idxs)-1]].Pos
		for center := minPos; center <= maxPos; center += int64(step) {
			start := center - int64(window)/2
			if start < 1 {
				start = 1
			}
			end := center + int64(window)/2
			stats := map[string]float64{}
			n := 0
			for _, name := range statNames {
				if name == "brm" {
					stats[name] = localLinear(markers, idxs, center, bandwidth, statGetter(name))
				} else {
					v, count := kernelMean(markers, idxs, center, bandwidth, statGetter(name))
					stats[name] = v
					if count > n {
						n = count
					}
				}
			}
			if n > 0 {
				points = append(points, smoothPoint{Chrom: chrom, Start: start, End: end, Center: center, MarkerN: n, Stats: stats})
			}
		}
	}
	sort.Slice(points, func(i, j int) bool {
		if points[i].Chrom == points[j].Chrom {
			return points[i].Center < points[j].Center
		}
		return points[i].Chrom < points[j].Chrom
	})
	return points
}

func splitByChrom(markers []marker) map[string][]int {
	out := map[string][]int{}
	for i, m := range markers {
		out[m.Chrom] = append(out[m.Chrom], i)
	}
	return out
}

func statGetter(name string) func(marker) float64 {
	return func(m marker) float64 {
		switch name {
		case "snp_index_high":
			return m.HighIndex
		case "snp_index_low":
			return m.LowIndex
		case "delta_snp_index":
			return m.Delta
		case "g_statistic":
			return m.G
		case "euclidean_distance":
			return m.ED
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

func kernelMean(markers []marker, idxs []int, center int64, bandwidth float64, get func(marker) float64) (float64, int) {
	sumW, sumY := 0.0, 0.0
	count := 0
	for _, idx := range idxs {
		dist := float64(markers[idx].Pos - center)
		w := math.Exp(-0.5 * math.Pow(dist/bandwidth, 2))
		if w < 0.001 {
			continue
		}
		y := get(markers[idx])
		if math.IsNaN(y) || math.IsInf(y, 0) {
			continue
		}
		sumW += w
		sumY += w * y
		count++
	}
	if sumW == 0 {
		return math.NaN(), count
	}
	return sumY / sumW, count
}

func localLinear(markers []marker, idxs []int, center int64, bandwidth float64, get func(marker) float64) float64 {
	var sw, sx, sy, sxx, sxy float64
	for _, idx := range idxs {
		x := float64(markers[idx].Pos - center)
		w := math.Exp(-0.5 * math.Pow(x/bandwidth, 2))
		if w < 0.001 {
			continue
		}
		y := get(markers[idx])
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

func simulateThresholds(markers []marker, cfg utils.AnalysisConfig, mode analysisMode) map[string]map[float64]float64 {
	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	values := map[string][]float64{}
	for _, name := range statNames {
		values[name] = []float64{}
	}
	progressEvery := cfg.Rep / 10
	if progressEvery < 1 {
		progressEvery = 1
	}
	for rep := 0; rep < cfg.Rep; rep++ {
		sim := make([]marker, len(markers))
		for i, m := range markers {
			sim[i] = simulateMarker(m, rng, mode)
		}
		assignBRM(sim, cfg.WindowSize)
		points := smoothMarkers(sim, cfg.WindowSize, cfg.StepSize)
		for _, p := range points {
			for _, name := range statNames {
				v := p.Stats[name]
				if math.IsNaN(v) || math.IsInf(v, 0) {
					continue
				}
				if signedStat(name) {
					v = math.Abs(v)
				}
				values[name] = append(values[name], v)
			}
		}
		if (rep+1)%progressEvery == 0 || rep+1 == cfg.Rep {
			color.White("  Simulations completed: %d/%d", rep+1, cfg.Rep)
		}
	}
	out := map[string]map[float64]float64{}
	for _, name := range statNames {
		out[name] = map[float64]float64{}
		sort.Float64s(values[name])
		for _, alpha := range cfg.Alphas {
			out[name][alpha] = percentile(values[name], 1-alpha)
			color.White("  Threshold %s at %.2f percentile: %s", name, (1-alpha)*100, ftoa(out[name][alpha]))
		}
	}
	return out
}

func simulateMarker(m marker, rng *rand.Rand, mode analysisMode) marker {
	out := m
	if mode == modeTwoParentsTwoBulks || mode == modeBulksOnly {
		p := ratio(m.High.Alt+m.Low.Alt, m.High.Ref+m.High.Alt+m.Low.Ref+m.Low.Alt)
		hAlt := binomial(rng, m.High.Ref+m.High.Alt, p)
		lAlt := binomial(rng, m.Low.Ref+m.Low.Alt, p)
		out.High = sampleDepth{Alt: hAlt, Ref: m.High.Ref + m.High.Alt - hAlt, DP: m.High.DP}
		out.Low = sampleDepth{Alt: lAlt, Ref: m.Low.Ref + m.Low.Alt - lAlt, DP: m.Low.DP}
		calcTwoBulkStats(&out)
		return out
	}
	n := m.One.Ref + m.One.Alt
	alt := binomial(rng, n, 0.5)
	out.One = sampleDepth{Alt: alt, Ref: n - alt, DP: m.One.DP}
	calcOneBulkStats(&out)
	return out
}

func binomial(rng *rand.Rand, n int, p float64) int {
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	x := 0
	for i := 0; i < n; i++ {
		if rng.Float64() < p {
			x++
		}
	}
	return x
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
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	if lo == hi {
		return vals[lo]
	}
	return vals[lo] + (vals[hi]-vals[lo])*(pos-float64(lo))
}

func attachThresholds(points []smoothPoint, thresholds map[string]map[float64]float64) {
	for i := range points {
		points[i].Thresholds = thresholds
	}
}

func attachRobustZ(points []smoothPoint) {
	for _, name := range statNames {
		var vals []float64
		for _, p := range points {
			v := p.Stats[name]
			if !math.IsNaN(v) && !math.IsInf(v, 0) {
				vals = append(vals, v)
			}
		}
		med := median(vals)
		var dev []float64
		for _, v := range vals {
			dev = append(dev, math.Abs(v-med))
		}
		mad := median(dev)
		if mad < eps {
			mad = 1
		}
		for i := range points {
			if points[i].Z == nil {
				points[i].Z = map[string]float64{}
			}
			points[i].Z[name] = (points[i].Stats[name] - med) / (1.4826 * mad)
		}
	}
}

func median(vals []float64) float64 {
	if len(vals) == 0 {
		return math.NaN()
	}
	cp := append([]float64(nil), vals...)
	sort.Float64s(cp)
	mid := len(cp) / 2
	if len(cp)%2 == 0 {
		return (cp[mid-1] + cp[mid]) / 2
	}
	return cp[mid]
}

func signedStat(name string) bool {
	return name == "delta_snp_index" || name == "brm"
}

func detectThresholdQTLs(points []smoothPoint, minWidth, mergeDist int64) []qtl {
	var all []qtl
	for _, stat := range statNames {
		if stat == "snp_index_high" || stat == "snp_index_low" {
			continue
		}
		alphas := sortedAlphas(points, stat)
		for _, alpha := range alphas {
			label := fmt.Sprintf("p%.2g", alpha)
			all = append(all, detectIntervals(points, stat, label, minWidth, mergeDist, func(p smoothPoint) (float64, bool) {
				v := p.Stats[stat]
				t := p.Thresholds[stat][alpha]
				if signedStat(stat) {
					v = math.Abs(v)
				}
				return v, !math.IsNaN(t) && v >= t
			})...)
		}
	}
	return all
}

func sortedAlphas(points []smoothPoint, stat string) []float64 {
	if len(points) == 0 || points[0].Thresholds == nil {
		return nil
	}
	var alphas []float64
	for a := range points[0].Thresholds[stat] {
		alphas = append(alphas, a)
	}
	sort.Float64s(alphas)
	return alphas
}

func detectZQTLs(points []smoothPoint, minWidth, mergeDist int64) []qtl {
	var all []qtl
	for _, stat := range statNames {
		if stat == "snp_index_high" || stat == "snp_index_low" {
			continue
		}
		for _, cutoff := range []float64{2, 3} {
			label := fmt.Sprintf("z%.0f", cutoff)
			all = append(all, detectIntervals(points, stat, label, minWidth, mergeDist, func(p smoothPoint) (float64, bool) {
				v := p.Z[stat]
				return math.Abs(v), math.Abs(v) >= cutoff
			})...)
		}
	}
	return all
}

func detectIntervals(points []smoothPoint, stat, threshold string, minWidth, mergeDist int64, pass func(smoothPoint) (float64, bool)) []qtl {
	var out []qtl
	var cur *qtl
	for _, p := range points {
		v, ok := pass(p)
		if !ok {
			if cur != nil {
				out = appendInterval(out, *cur, minWidth, mergeDist)
				cur = nil
			}
			continue
		}
		if cur == nil || cur.Chrom != p.Chrom || p.Start-cur.End > mergeDist {
			if cur != nil {
				out = appendInterval(out, *cur, minWidth, mergeDist)
			}
			cur = &qtl{Chrom: p.Chrom, Start: p.Start, End: p.End, Stat: stat, Threshold: threshold, PeakPos: p.Center, PeakValue: v, Markers: p.MarkerN}
			continue
		}
		if p.End > cur.End {
			cur.End = p.End
		}
		cur.Markers += p.MarkerN
		if v > cur.PeakValue {
			cur.PeakValue = v
			cur.PeakPos = p.Center
		}
	}
	if cur != nil {
		out = appendInterval(out, *cur, minWidth, mergeDist)
	}
	return out
}

func appendInterval(out []qtl, item qtl, minWidth, _ int64) []qtl {
	if item.End-item.Start+1 >= minWidth {
		return append(out, item)
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
	_ = w.Write([]string{"chrom", "pos", "ref", "alt", "high_ref_depth", "high_alt_depth", "low_ref_depth", "low_alt_depth", "snp_index_high", "snp_index_low", "delta_snp_index", "g_statistic", "euclidean_distance", "lod", "beta_binomial_bf", "brm"})
	for _, m := range markers {
		_ = w.Write([]string{m.Chrom, i64(m.Pos), m.Ref, m.Alt, itoa(m.High.Ref), itoa(m.High.Alt), itoa(m.Low.Ref), itoa(m.Low.Alt), ftoa(m.HighIndex), ftoa(m.LowIndex), ftoa(m.Delta), ftoa(m.G), ftoa(m.ED), ftoa(m.LOD), ftoa(m.BF), ftoa(m.BRM)})
	}
	return w.Error()
}

func writeSmooth(path string, points []smoothPoint) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()
	header := []string{"chrom", "start", "end", "center", "marker_count"}
	for _, name := range statNames {
		header = append(header, name, name+"_robust_z")
	}
	_ = w.Write(header)
	for _, p := range points {
		row := []string{p.Chrom, i64(p.Start), i64(p.End), i64(p.Center), itoa(p.MarkerN)}
		for _, name := range statNames {
			row = append(row, ftoa(p.Stats[name]), ftoa(p.Z[name]))
		}
		_ = w.Write(row)
	}
	return w.Error()
}

func writeThresholds(path string, thresholds map[string]map[float64]float64) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()
	_ = w.Write([]string{"stat", "alpha", "percentile", "threshold"})
	for _, stat := range statNames {
		var alphas []float64
		for a := range thresholds[stat] {
			alphas = append(alphas, a)
		}
		sort.Float64s(alphas)
		for _, a := range alphas {
			_ = w.Write([]string{stat, ftoa(a), ftoa(1 - a), ftoa(thresholds[stat][a])})
		}
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
	_ = w.Write([]string{"chrom", "start", "end", "stat", "threshold", "peak_pos", "peak_value", "supporting_window_markers"})
	for _, q := range qtls {
		_ = w.Write([]string{q.Chrom, i64(q.Start), i64(q.End), q.Stat, q.Threshold, i64(q.PeakPos), ftoa(q.PeakValue), itoa(q.Markers)})
	}
	return w.Error()
}

func plotAll(outDir string, points []smoothPoint) error {
	for _, stat := range statNames {
		if err := plotStat(filepath.Join(outDir, "plots", stat+".html"), points, stat, false); err != nil {
			return err
		}
		if err := plotStat(filepath.Join(outDir, "plots", stat+"_robust_z.html"), points, stat, true); err != nil {
			return err
		}
	}
	return nil
}

func plotStat(path string, points []smoothPoint, stat string, z bool) error {
	line := charts.NewLine()
	title := stat
	if z {
		title += " robust z-score"
	}
	line.SetGlobalOptions(
		charts.WithTitleOpts(opts.Title{Title: title}),
		charts.WithXAxisOpts(opts.XAxis{Name: "genome windows", AxisLabel: &opts.AxisLabel{Show: opts.Bool(false)}}),
		charts.WithYAxisOpts(opts.YAxis{Name: stat}),
		charts.WithLegendOpts(opts.Legend{Show: opts.Bool(true)}),
	)
	var x []string
	var chroms []string
	for _, p := range points {
		if len(chroms) == 0 || chroms[len(chroms)-1] != p.Chrom {
			chroms = append(chroms, p.Chrom)
		}
		x = append(x, fmt.Sprintf("%s:%d", p.Chrom, p.Center))
	}
	series := map[string][]opts.LineData{}
	for _, chrom := range chroms {
		series[chrom] = make([]opts.LineData, len(points))
		for i := range series[chrom] {
			series[chrom][i] = opts.LineData{Value: nil}
		}
	}
	for i, p := range points {
		v := p.Stats[stat]
		if z {
			v = p.Z[stat]
		}
		series[p.Chrom][i] = opts.LineData{Value: cleanFloat(v)}
	}
	line.SetXAxis(x)
	for _, chrom := range chroms {
		line.AddSeries(chrom, series[chrom])
	}
	if z {
		addFlatSeries(line, x, "z=2", 2)
		addFlatSeries(line, x, "z=-2", -2)
		addFlatSeries(line, x, "z=3", 3)
		addFlatSeries(line, x, "z=-3", -3)
	} else if len(points) > 0 && points[0].Thresholds != nil {
		var alphas []float64
		for a := range points[0].Thresholds[stat] {
			alphas = append(alphas, a)
		}
		sort.Float64s(alphas)
		for _, a := range alphas {
			addFlatSeries(line, x, fmt.Sprintf("p%.2g", a), points[0].Thresholds[stat][a])
			if signedStat(stat) {
				addFlatSeries(line, x, fmt.Sprintf("-p%.2g", a), -points[0].Thresholds[stat][a])
			}
		}
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return line.Render(f)
}

func addFlatSeries(line *charts.Line, x []string, name string, value float64) {
	data := make([]opts.LineData, len(x))
	for i := range data {
		data[i] = opts.LineData{Value: cleanFloat(value)}
	}
	line.AddSeries(name, data)
}

func cleanFloat(v float64) interface{} {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return nil
	}
	return v
}

func ftoa(v float64) string {
	if math.IsNaN(v) {
		return "NA"
	}
	if math.IsInf(v, 1) {
		return "Inf"
	}
	if math.IsInf(v, -1) {
		return "-Inf"
	}
	return strconv.FormatFloat(v, 'g', 8, 64)
}

func i64(v int64) string {
	return strconv.FormatInt(v, 10)
}

func itoa(v int) string {
	return strconv.Itoa(v)
}
