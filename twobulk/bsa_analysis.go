package twobulk

//
//import (
//	"bufio"
//	"fmt"
//	"math"
//	"os"
//	"os/exec"
//	"runtime"
//	"sort"
//	"strconv"
//	"strings"
//	"sync"
//	"sync/atomic"
//	"time"
//)
//
//// SNP represents a single SNP record with allele depth information
//type SNP struct {
//	Chromosome string
//	Position   int
//	Ref        string
//	Alt        string
//	Qual       float64
//	Sample1    SampleData
//	Sample2    SampleData
//}
//
//// SampleData holds allele count information for one sample
//type SampleData struct {
//	RefCount int
//	AltCount int
//	DP       int
//	GT       string
//}
//
//// Result holds all calculated statistics for a SNP
//type Result struct {
//	SNP
//	SNPIndex1   float64
//	SNPIndex2   float64
//	DeltaSNP    float64
//	ED          float64
//	ED4         float64  // ED^4 (used for enhanced peak detection)
//	GStat       float64
//	DeltaSmooth float64
//	EDSmooth    float64
//	GStatSmooth float64
//}
//
//// Config holds program configuration
//type Config struct {
//	VCFPath      string
//	Sample1Name  string
//	Sample2Name  string
//	WindowSize   int      // Physical window size in bp
//	NumWorkers   int
//	MinDepth     int
//	MinQuality   float64
//	OutputPrefix string
//	MinSNPIndex  float64  // Minimum SNP index to consider
//}
//
//func main() {
//	if len(os.Args) < 4 {
//		fmt.Fprintf(os.Stderr, `BSA Analysis Tool - High Performance Go Implementation
//
//Usage: bsa <vcf.gz> <sample1> <sample2> [output_prefix]
//
//Arguments:
//  vcf.gz        Path to bgzipped VCF file (must be indexed)
//  sample1       Name of first bulk sample
//  sample2       Name of second bulk sample
//  output_prefix Output file prefix (default: bsa_output)
//
//Example:
//  bsa variants.vcf.gz bulk_high bulk_low my_analysis
//
//Output files:
//  <prefix>_results.tsv    - All SNPs with raw and smoothed statistics
//  <prefix>_summary.tsv    - Chromosome-wide summary with peak positions
//  <prefix>_peaks.tsv      - Detected QTL peaks
//
//Statistics calculated:
//  - SNP Index (allele frequency) for each bulk
//  - Delta SNP Index (difference between bulks)
//  - Euclidean Distance (ED)
//  - G-statistic (likelihood ratio test)
//
//Smoothing methods:
//  - Delta SNP: Tricube kernel (Nadaraya-Watson)
//  - ED: Gaussian kernel-weighted smoothing
//  - G-statistic: Tricube kernel (likelihood-based)
//`)
//		os.Exit(1)
//	}
//
//	config := Config{
//		VCFPath:      os.Args[1],
//		Sample1Name:  os.Args[2],
//		Sample2Name:  os.Args[3],
//		WindowSize:   2000000, // 2Mb default (adjustable)
//		NumWorkers:   runtime.NumCPU(),
//		MinDepth:     10,
//		MinQuality:   20,
//		OutputPrefix: "bsa_output",
//		MinSNPIndex:  0.0,
//	}
//	if len(os.Args) > 4 {
//		config.OutputPrefix = os.Args[4]
//	}
//
//	fmt.Printf("BSA Analysis - High Performance Implementation\n")
//	fmt.Printf("==============================================\n")
//	fmt.Printf("VCF File:      %s\n", config.VCFPath)
//	fmt.Printf("Sample 1:      %s\n", config.Sample1Name)
//	fmt.Printf("Sample 2:      %s\n", config.Sample2Name)
//	fmt.Printf("Workers:       %d (CPU cores)\n", config.NumWorkers)
//	fmt.Printf("Window Size:   %d bp\n", config.WindowSize)
//	fmt.Printf("Min Depth:     %d\n", config.MinDepth)
//	fmt.Printf("Min Quality:   %.1f\n", config.MinQuality)
//	fmt.Println()
//
//	start := time.Now()
//
//	// Step 1: Parse VCF using bcftools
//	fmt.Println("[1/4] Parsing VCF with bcftools...")
//	snps, err := parseVCF(config)
//	if err != nil {
//		fmt.Fprintf(os.Stderr, "Error parsing VCF: %v\n", err)
//		os.Exit(1)
//	}
//	fmt.Printf("      Parsed %d SNPs (%.2fs)\n", len(snps), time.Since(start).Seconds())
//
//	if len(snps) == 0 {
//		fmt.Println("No valid SNPs found. Check sample names and VCF format.")
//		os.Exit(1)
//	}
//
//	// Step 2: Calculate statistics in parallel
//	fmt.Println("[2/4] Calculating statistics...")
//	step2 := time.Now()
//	results := calculateStatsParallel(snps, config)
//	fmt.Printf("      Calculated for %d SNPs (%.2fs)\n", len(results), time.Since(step2).Seconds())
//
//	// Step 3: Apply smoothing
//	fmt.Println("[3/4] Applying smoothing...")
//	step3 := time.Now()
//	results = smoothAll(results, config)
//	fmt.Printf("      Smoothed %d results (%.2fs)\n", len(results), time.Since(step3).Seconds())
//
//	// Step 4: Write outputs
//	fmt.Println("[4/4] Writing output files...")
//	step4 := time.Now()
//	if err := writeResults(results, config); err != nil {
//		fmt.Fprintf(os.Stderr, "Error writing results: %v\n", err)
//		os.Exit(1)
//	}
//	fmt.Printf("      Files written (%.2fs)\n", time.Since(step4).Seconds())
//
//	fmt.Printf("\nTotal time: %.2fs\n", time.Since(start).Seconds())
//	fmt.Println("Done!")
//}
//
//// parseVCF reads VCF using bcftools query for maximum performance
//func parseVCF(config Config) ([]SNP, error) {
//	// Use bcftools query to extract only needed fields
//	// Format: CHROM POS REF ALT QUAL [GT AD DP]
//	cmd := exec.Command("bcftools", "query", "-s", config.Sample1Name+","+config.Sample2Name,
//		"-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t%AD\t%DP\n]\n",
//		config.VCFPath)
//
//	stdout, err := cmd.StdoutPipe()
//	if err != nil {
//		return nil, fmt.Errorf("bcftools stdout pipe: %w", err)
//	}
//
//	stderr, err := cmd.StderrPipe()
//	if err != nil {
//		return nil, fmt.Errorf("bcftools stderr pipe: %w", err)
//	}
//
//	if err := cmd.Start(); err != nil {
//		return nil, fmt.Errorf("bcftools start: %w", err)
//	}
//
//	// Capture stderr for debugging
//	var stderrBuf strings.Builder
//	go func() {
//		scanner := bufio.NewScanner(stderr)
//		for scanner.Scan() {
//			stderrBuf.WriteString(scanner.Text() + "\n")
//		}
//	}()
//
//	var snps []SNP
//	scanner := bufio.NewScanner(stdout)
//	// Increase buffer for large lines
//	const maxCapacity = 1024 * 1024 // 1MB
//	buf := make([]byte, maxCapacity)
//	scanner.Buffer(buf, maxCapacity)
//
//	lineNum := 0
//	for scanner.Scan() {
//		line := scanner.Text()
//		lineNum++
//
//		if line == "" || strings.HasPrefix(line, "#") {
//			continue
//		}
//
//		fields := strings.Split(line, "\t")
//		if len(fields) < 8 {
//			continue
//		}
//
//		pos, err := strconv.Atoi(fields[1])
//		if err != nil {
//			continue
//		}
//
//		qual, err := strconv.ParseFloat(fields[4], 64)
//		if err != nil {
//			qual = 0
//		}
//
//		if qual < config.MinQuality {
//			continue
//		}
//
//		// Parse sample 1 data (fields 5,6,7)
//		s1 := parseSampleData(fields[5:8])
//		// Parse sample 2 data (fields 8,9,10)
//		s2 := parseSampleData(fields[8:11])
//
//		if s1.DP < config.MinDepth || s2.DP < config.MinDepth {
//			continue
//		}
//
//		// Skip if no alt reads in either sample
//		if s1.AltCount == 0 && s2.AltCount == 0 {
//			continue
//		}
//
//		snps = append(snps, SNP{
//			Chromosome: fields[0],
//			Position:   pos,
//			Ref:        fields[2],
//			Alt:        fields[3],
//			Qual:       qual,
//			Sample1:    s1,
//			Sample2:    s2,
//		})
//	}
//
//	if err := scanner.Err(); err != nil {
//		return nil, fmt.Errorf("scanner error: %w", err)
//	}
//
//	if err := cmd.Wait(); err != nil {
//		stderrStr := stderrBuf.String()
//		if stderrStr != "" {
//			return nil, fmt.Errorf("bcftools error: %v\nstderr: %s", err, stderrStr)
//		}
//		return nil, fmt.Errorf("bcftools error: %w", err)
//	}
//
//	return snps, nil
//}
//
//func parseSampleData(fields []string) SampleData {
//	if len(fields) < 3 {
//		return SampleData{}
//	}
//
//	gt := fields[0]
//	adStr := strings.Split(fields[1], ",")
//	dp, _ := strconv.Atoi(fields[2])
//
//	var refCount, altCount int
//	if len(adStr) >= 2 {
//		refCount, _ = strconv.Atoi(adStr[0])
//		altCount, _ = strconv.Atoi(adStr[1])
//	}
//
//	return SampleData{
//		GT:       gt,
//		RefCount: refCount,
//		AltCount: altCount,
//		DP:       dp,
//	}
//}
//
//// calculateStatsParallel calculates all statistics concurrently using worker pool
//func calculateStatsParallel(snps []SNP, config Config) []Result {
//	n := len(snps)
//	results := make([]Result, n)
//
//	// Use atomic counter for progress tracking
//	var processed int64
//	var mu sync.Mutex
//	lastProgress := 0
//
//	// Determine optimal chunk size
//	chunkSize := 1000
//	if n < config.NumWorkers*chunkSize {
//		chunkSize = (n + config.NumWorkers - 1) / config.NumWorkers
//		if chunkSize < 100 {
//			chunkSize = 100
//		}
//	}
//
//	numChunks := (n + chunkSize - 1) / chunkSize
//
//	var wg sync.WaitGroup
//	for i := 0; i < numChunks; i++ {
//		start := i * chunkSize
//		end := start + chunkSize
//		if end > n {
//			end = n
//		}
//
//		wg.Add(1)
//		go func(s, e int) {
//			defer wg.Done()
//			for j := s; j < e; j++ {
//				results[j] = calculateSingleStats(snps[j])
//			}
//
//			// Progress tracking
//			current := atomic.AddInt64(&processed, int64(e-s))
//			progress := int(float64(current) / float64(n) * 100)
//			mu.Lock()
//			if progress > lastProgress && progress%10 == 0 {
//				lastProgress = progress
//				fmt.Printf("      Progress: %d%%\n", progress)
//			}
//			mu.Unlock()
//		}(start, end)
//	}
//
//	wg.Wait()
//	return results
//}
//
//// calculateSingleStats computes all per-SNP statistics
//func calculateSingleStats(snp SNP) Result {
//	s1 := snp.Sample1
//	s2 := snp.Sample2
//
//	// SNP Index = alt allele frequency (0 to 1)
//	snpIndex1 := float64(s1.AltCount) / float64(s1.DP)
//	snpIndex2 := float64(s2.AltCount) / float64(s2.DP)
//
//	// Delta SNP Index
//	delta := snpIndex1 - snpIndex2
//
//	// Euclidean Distance (ED)
//	// ED = sqrt((p1 - p2)^2 + (q1 - q2)^2)
//	// where p = ref frequency, q = alt frequency
//	p1 := float64(s1.RefCount) / float64(s1.DP)
//	q1 := float64(s1.AltCount) / float64(s1.DP)
//	p2 := float64(s2.RefCount) / float64(s2.DP)
//	q2 := float64(s2.AltCount) / float64(s2.DP)
//
//	dp := p1 - p2
//	dq := q1 - q2
//	ed := math.Sqrt(dp*dp + dq*dq)
//	ed4 := math.Pow(ed, 4) // ED^4 for enhanced peak detection
//
//	// G-statistic (likelihood ratio test)
//	// G = 2 * sum(O_i * ln(O_i / E_i))
//	total1 := float64(s1.DP)
//	total2 := float64(s2.DP)
//	totalRef := float64(s1.RefCount + s2.RefCount)
//	totalAlt := float64(s1.AltCount + s2.AltCount)
//	grandTotal := total1 + total2
//
//	// Expected counts under null hypothesis
//	eRef1 := total1 * totalRef / grandTotal
//	eAlt1 := total1 * totalAlt / grandTotal
//	eRef2 := total2 * totalRef / grandTotal
//	eAlt2 := total2 * totalAlt / grandTotal
//
//	gStat := 0.0
//	obs := []float64{float64(s1.RefCount), float64(s1.AltCount), float64(s2.RefCount), float64(s2.AltCount)}
//	exp := []float64{eRef1, eAlt1, eRef2, eAlt2}
//
//	for i := 0; i < 4; i++ {
//		if obs[i] > 0 && exp[i] > 0 {
//			gStat += 2 * obs[i] * math.Log(obs[i]/exp[i])
//		}
//	}
//
//	return Result{
//		SNP:       snp,
//		SNPIndex1: snpIndex1,
//		SNPIndex2: snpIndex2,
//		DeltaSNP:  delta,
//		ED:        ed,
//		ED4:       ed4,
//		GStat:     gStat,
//	}
//}
//
//// smoothAll applies all smoothing methods per chromosome in parallel
//func smoothAll(results []Result, config Config) []Result {
//	if len(results) == 0 {
//		return results
//	}
//
//	// Group by chromosome
//	chromGroups := groupByChromosome(results)
//
//	var wg sync.WaitGroup
//	var mu sync.Mutex
//	smoothed := make([]Result, 0, len(results))
//
//	for chrom, snps := range chromGroups {
//		wg.Add(1)
//		go func(c string, s []Result) {
//			defer wg.Done()
//
//			// Sort by position within chromosome
//			sort.Slice(s, func(i, j int) bool {
//				return s[i].Position < s[j].Position
//			})
//
//			// Apply smoothing methods
//			smoothDeltaTricube(s, config.WindowSize)
//			smoothEDGaussian(s, config.WindowSize)
//			smoothGStatTricube(s, config.WindowSize)
//
//			mu.Lock()
//			smoothed = append(smoothed, s...)
//			mu.Unlock()
//		}(chrom, snps)
//	}
//
//	wg.Wait()
//	return smoothed
//}
//
//func groupByChromosome(results []Result) map[string][]Result {
//	groups := make(map[string][]Result)
//	for _, r := range results {
//		groups[r.Chromosome] = append(groups[r.Chromosome], r)
//	}
//	return groups
//}
//
//// smoothDeltaTricube applies tricube kernel smoothing to Delta SNP index
//// Uses Nadaraya-Watson estimator with tricube weighting
//func smoothDeltaTricube(snps []Result, windowSize int) {
//	n := len(snps)
//	if n == 0 {
//		return
//	}
//
//	halfWindow := windowSize / 2
//
//	for i := 0; i < n; i++ {
//		var weightedSum, weightSum float64
//		centerPos := snps[i].Position
//
//		for j := 0; j < n; j++ {
//			dist := float64(snps[j].Position - centerPos)
//			if math.Abs(dist) > float64(halfWindow) {
//				continue
//			}
//
//			// Tricube kernel: w(u) = (1 - |u|^3)^3 where u = distance / bandwidth
//			u := math.Abs(dist) / float64(halfWindow)
//			if u >= 1 {
//				continue
//			}
//
//			weight := math.Pow(1-math.Pow(u, 3), 3)
//
//			weightedSum += snps[j].DeltaSNP * weight
//			weightSum += weight
//		}
//
//		if weightSum > 0 {
//			snps[i].DeltaSmooth = weightedSum / weightSum
//		}
//	}
//}
//
//// smoothEDGaussian applies Gaussian kernel-weighted smoothing to ED
//func smoothEDGaussian(snps []Result, windowSize int) {
//	n := len(snps)
//	if n == 0 {
//		return
//	}
//
//	halfWindow := windowSize / 2
//	// Sigma = bandwidth / 3 for Gaussian kernel
//	sigma := float64(halfWindow) / 3.0
//	// Precompute denominator for efficiency
//	twoSigmaSq := 2 * sigma * sigma
//
//	for i := 0; i < n; i++ {
//		var weightedSum, weightSum float64
//		centerPos := snps[i].Position
//
//		for j := 0; j < n; j++ {
//			dist := float64(snps[j].Position - centerPos)
//			if math.Abs(dist) > float64(halfWindow) {
//				continue
//			}
//
//			// Gaussian kernel: w(d) = exp(-d^2 / (2*sigma^2))
//			weight := math.Exp(-(dist * dist) / twoSigmaSq)
//
//			weightedSum += snps[j].ED * weight
//			weightSum += weight
//		}
//
//		if weightSum > 0 {
//			snps[i].EDSmooth = weightedSum / weightSum
//		}
//	}
//}
//
//// smoothGStatTricube applies tricube smoothing to G-statistic
//// This is the standard approach for likelihood-based statistics in BSA
//func smoothGStatTricube(snps []Result, windowSize int) {
//	n := len(snps)
//	if n == 0 {
//		return
//	}
//
//	halfWindow := windowSize / 2
//
//	for i := 0; i < n; i++ {
//		var weightedSum, weightSum float64
//		centerPos := snps[i].Position
//
//		for j := 0; j < n; j++ {
//			dist := float64(snps[j].Position - centerPos)
//			if math.Abs(dist) > float64(halfWindow) {
//				continue
//			}
//
//			// Tricube kernel
//			u := math.Abs(dist) / float64(halfWindow)
//			if u >= 1 {
//				continue
//			}
//
//			weight := math.Pow(1-math.Pow(u, 3), 3)
//
//			weightedSum += snps[j].GStat * weight
//			weightSum += weight
//		}
//
//		if weightSum > 0 {
//			snps[i].GStatSmooth = weightedSum / weightSum
//		}
//	}
//}
//
//// writeResults writes all output files
//func writeResults(results []Result, config Config) error {
//	// Sort results by chromosome and position
//	sort.Slice(results, func(i, j int) bool {
//		if results[i].Chromosome != results[j].Chromosome {
//			return results[i].Chromosome < results[j].Chromosome
//		}
//		return results[i].Position < results[j].Position
//	})
//
//	// Write main results file
//	mainFile := config.OutputPrefix + "_results.tsv"
//	f, err := os.Create(mainFile)
//	if err != nil {
//		return fmt.Errorf("create results file: %w", err)
//	}
//	defer f.Close()
//
//	w := bufio.NewWriterSize(f, 1024*1024) // 1MB buffer
//	fmt.Fprintf(w, "# BSA Analysis Results\n")
//	fmt.Fprintf(w, "# Sample1: %s, Sample2: %s\n", config.Sample1Name, config.Sample2Name)
//	fmt.Fprintf(w, "# Window: %d bp\n", config.WindowSize)
//	fmt.Fprintf(w, "Chromosome\tPosition\tRef\tAlt\tQual\t")
//	fmt.Fprintf(w, "SNPIndex1\tSNPIndex2\tDeltaSNP\t")
//	fmt.Fprintf(w, "ED\tED4\tGStat\t")
//	fmt.Fprintf(w, "DeltaSmooth\tEDSmooth\tGStatSmooth\n")
//
//	for _, r := range results {
//		fmt.Fprintf(w, "%s\t%d\t%s\t%s\t%.1f\t",
//			r.Chromosome, r.Position, r.Ref, r.Alt, r.Qual)
//		fmt.Fprintf(w, "%.4f\t%.4f\t%.4f\t",
//			r.SNPIndex1, r.SNPIndex2, r.DeltaSNP)
//		fmt.Fprintf(w, "%.4f\t%.4f\t%.4f\t",
//			r.ED, r.ED4, r.GStat)
//		fmt.Fprintf(w, "%.4f\t%.4f\t%.4f\n",
//			r.DeltaSmooth, r.EDSmooth, r.GStatSmooth)
//	}
//	w.Flush()
//	fmt.Printf("      %s\n", mainFile)
//
//	// Write chromosome summary
//	if err := writeSummary(results, config); err != nil {
//		return err
//	}
//
//	// Write peaks
//	if err := writePeaks(results, config); err != nil {
//		return err
//	}
//
//	return nil
//}
//
//func writeSummary(results []Result, config Config) error {
//	chromStats := make(map[string]struct {
//		maxDelta     float64
//		maxDeltaPos  int
//		maxED        float64
//		maxEDPos     int
//		maxGStat     float64
//		maxGStatPos  int
//		snpCount     int
//		minPos       int
//		maxPos       int
//	})
//
//	for _, r := range results {
//		stats := chromStats[r.Chromosome]
//		stats.snpCount++
//
//		if stats.minPos == 0 || r.Position < stats.minPos {
//			stats.minPos = r.Position
//		}
//		if r.Position > stats.maxPos {
//			stats.maxPos = r.Position
//		}
//
//		absDelta := math.Abs(r.DeltaSmooth)
//		if absDelta > math.Abs(stats.maxDelta) {
//			stats.maxDelta = r.DeltaSmooth
//			stats.maxDeltaPos = r.Position
//		}
//		if r.EDSmooth > stats.maxED {
//			stats.maxED = r.EDSmooth
//			stats.maxEDPos = r.Position
//		}
//		if r.GStatSmooth > stats.maxGStat {
//			stats.maxGStat = r.GStatSmooth
//			stats.maxGStatPos = r.Position
//		}
//
//		chromStats[r.Chromosome] = stats
//	}
//
//	summaryFile := config.OutputPrefix + "_summary.tsv"
//	f, err := os.Create(summaryFile)
//	if err != nil {
//		return err
//	}
//	defer f.Close()
//
//	w := bufio.NewWriter(f)
//	fmt.Fprintf(w, "# Chromosome Summary Statistics\n")
//	fmt.Fprintf(w, "Chromosome\tSNP_Count\tStart\tEnd\t")
//	fmt.Fprintf(w, "Max_Delta\tMax_Delta_Pos\t")
//	fmt.Fprintf(w, "Max_ED\tMax_ED_Pos\t")
//	fmt.Fprintf(w, "Max_GStat\tMax_GStat_Pos\n")
//
//	var chroms []string
//	for c := range chromStats {
//		chroms = append(chroms, c)
//	}
//	sort.Strings(chroms)
//
//	for _, c := range chroms {
//		s := chromStats[c]
//		fmt.Fprintf(w, "%s\t%d\t%d\t%d\t", c, s.snpCount, s.minPos, s.maxPos)
//		fmt.Fprintf(w, "%.4f\t%d\t", s.maxDelta, s.maxDeltaPos)
//		fmt.Fprintf(w, "%.4f\t%d\t", s.maxED, s.maxEDPos)
//		fmt.Fprintf(w, "%.4f\t%d\n", s.maxGStat, s.maxGStatPos)
//	}
//	w.Flush()
//	fmt.Printf("      %s\n", summaryFile)
//	return nil
//}
//
//// writePeaks identifies and writes significant peaks
//func writePeaks(results []Result, config Config) error {
//	// Simple peak detection: local maxima in smoothed statistics
//	peaksFile := config.OutputPrefix + "_peaks.tsv"
//	f, err := os.Create(peaksFile)
//	if err != nil {
//		return err
//	}
//	defer f.Close()
//
//	w := bufio.NewWriter(f)
//	fmt.Fprintf(w, "# QTL Peaks (local maxima in smoothed statistics)\n")
//	fmt.Fprintf(w, "Chromosome\tPosition\tType\tValue\n")
//
//	// Group by chromosome
//	chromGroups := groupByChromosome(results)
//
//	for _, snps := range chromGroups {
//		if len(snps) < 3 {
//			continue
//		}
//
//		// Find peaks for each statistic
//		for i := 1; i < len(snps)-1; i++ {
//			// Delta peak
//			if math.Abs(snps[i].DeltaSmooth) > math.Abs(snps[i-1].DeltaSmooth) &&
//				math.Abs(snps[i].DeltaSmooth) > math.Abs(snps[i+1].DeltaSmooth) &&
//				math.Abs(snps[i].DeltaSmooth) > 0.3 {
//				fmt.Fprintf(w, "%s\t%d\tDeltaSNP\t%.4f\n",
//					snps[i].Chromosome, snps[i].Position, snps[i].DeltaSmooth)
//			}
//
//			// ED peak
//			if snps[i].EDSmooth > snps[i-1].EDSmooth &&
//				snps[i].EDSmooth > snps[i+1].EDSmooth &&
//				snps[i].EDSmooth > 0.3 {
//				fmt.Fprintf(w, "%s\t%d\tED\t%.4f\n",
//					snps[i].Chromosome, snps[i].Position, snps[i].EDSmooth)
//			}
//
//			// GStat peak
//			if snps[i].GStatSmooth > snps[i-1].GStatSmooth &&
//				snps[i].GStatSmooth > snps[i+1].GStatSmooth &&
//				snps[i].GStatSmooth > 10 {
//				fmt.Fprintf(w, "%s\t%d\tGStat\t%.4f\n",
//					snps[i].Chromosome, snps[i].Position, snps[i].GStatSmooth)
//			}
//		}
//	}
//
//	w.Flush()
//	fmt.Printf("      %s\n", peaksFile)
//	return nil
//}
