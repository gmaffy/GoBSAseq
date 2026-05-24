package cmd

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/best_of_both/run"
	"github.com/gmaffy/GoBSAseq/best_of_both/utils"
	"github.com/spf13/cobra"
)

var rootCmd = &cobra.Command{
	Use:   "GoBSAseq-best [variant.vcf]",
	Short: "Best-of-both BSA-seq analysis pipeline",
	Long:  "Best-of-both BSA-seq analysis pipeline combining robust VCF parsing, adaptive thresholds, consensus QTLs, BRM validation, clean logging, and interactive sample selection.",
	Run: func(cmd *cobra.Command, args []string) {
		variant, _ := cmd.Flags().GetString("variant")
		if strings.TrimSpace(variant) == "" {
			if len(args) == 1 {
				variant = args[0]
			} else {
				_ = cmd.Help()
				return
			}
		}

		cfg, hf, err := configFromFlags(cmd, variant)
		if err != nil {
			color.Red("%v", err)
			return
		}
		if err := run.Run(cfg, hf); err != nil {
			color.Red("%v", err)
		}
	},
}

func Execute() {
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

func configFromFlags(cmd *cobra.Command, variant string) (utils.AnalysisConfig, utils.HardFilterConfig, error) {
	parents, _ := cmd.Flags().GetString("parents")
	bulks, _ := cmd.Flags().GetString("bulks")
	parentsDepth, _ := cmd.Flags().GetString("parents-depth")
	bulksDepth, _ := cmd.Flags().GetString("bulks-depth")
	bulkSizes, _ := cmd.Flags().GetString("bulk-sizes")
	windowSize, _ := cmd.Flags().GetInt64("window-size")
	stepSize, _ := cmd.Flags().GetInt64("step-size")
	population, _ := cmd.Flags().GetString("population")
	rep, _ := cmd.Flags().GetInt("rep")
	alphas, _ := cmd.Flags().GetFloat64Slice("alpha")
	minQTL, _ := cmd.Flags().GetInt64("min-qtl-length")
	mergeDist, _ := cmd.Flags().GetInt64("merge-distance")
	outputDir, _ := cmd.Flags().GetString("out")

	selections := sampleSelections{}
	if !cmd.Flags().Changed("parents") && !cmd.Flags().Changed("bulks") && strings.TrimSpace(parents) == "" && strings.TrimSpace(bulks) == "" {
		prompted, err := promptSamplesFromVCF(variant)
		if err != nil {
			return utils.AnalysisConfig{}, utils.HardFilterConfig{}, err
		}
		selections = prompted
	} else {
		parentNames := splitCSV(parents)
		bulkNames := splitCSV(bulks)
		switch len(parentNames) {
		case 1:
			selections.OneParent = parentNames[0]
		case 2:
			selections.HighParent, selections.LowParent = parentNames[0], parentNames[1]
		default:
			if len(parentNames) > 2 {
				return utils.AnalysisConfig{}, utils.HardFilterConfig{}, fmt.Errorf("parents must be high,low or one sample")
			}
		}
		switch len(bulkNames) {
		case 1:
			selections.OneBulk = bulkNames[0]
		case 2:
			selections.HighBulk, selections.LowBulk = bulkNames[0], bulkNames[1]
		default:
			if len(bulkNames) > 2 {
				return utils.AnalysisConfig{}, utils.HardFilterConfig{}, fmt.Errorf("bulks must be high,low or one sample")
			}
		}
	}

	parentDepths, err := parseIntList(parentsDepth, 2)
	if err != nil {
		return utils.AnalysisConfig{}, utils.HardFilterConfig{}, fmt.Errorf("parents-depth: %w", err)
	}
	bulkDepths, err := parseIntList(bulksDepth, 2)
	if err != nil {
		return utils.AnalysisConfig{}, utils.HardFilterConfig{}, fmt.Errorf("bulks-depth: %w", err)
	}
	sizes, err := parseIntList(bulkSizes, 2)
	if err != nil {
		return utils.AnalysisConfig{}, utils.HardFilterConfig{}, fmt.Errorf("bulk-sizes: %w", err)
	}

	cfg := utils.AnalysisConfig{
		VCF: variant, Population: population, WindowSize: int(windowSize), StepSize: int(stepSize),
		Rep: rep, Alphas: alphas, MinQTLWidth: minQTL, MergeDistance: mergeDist, OutputDir: outputDir,
		HighParentName: selections.HighParent, LowParentName: selections.LowParent, OneParentName: selections.OneParent,
		HighBulkName: selections.HighBulk, LowBulkName: selections.LowBulk, OneBulkName: selections.OneBulk,
		HighParentDepth: parentDepths[0], LowParentDepth: parentDepths[1], OneParentDepth: parentDepths[0],
		HighBulkDepth: bulkDepths[0], LowBulkDepth: bulkDepths[1], OneBulkDepth: bulkDepths[0],
		HighBulkSize: sizes[0], LowBulkSize: sizes[1], OneBulkSize: sizes[0],
	}

	cfg.SnpEffDB, _ = cmd.Flags().GetString("snpEffDB")
	cfg.Gff, _ = cmd.Flags().GetString("gff")
	cfg.Cds, _ = cmd.Flags().GetString("cds")
	cfg.Ref, _ = cmd.Flags().GetString("reference")
	cfg.Protein, _ = cmd.Flags().GetString("protein")
	cfg.GeneDesc, _ = cmd.Flags().GetString("gene-descriptions")
	cfg.Prg, _ = cmd.Flags().GetString("prg")

	hf := utils.HardFilterConfig{}
	hf.SNP_QD_Min, _ = cmd.Flags().GetFloat64("min-QD-SNP")
	hf.SNP_QUAL_Min, _ = cmd.Flags().GetFloat64("min-QUAL-SNP")
	hf.SNP_SOR_Max, _ = cmd.Flags().GetFloat64("min-SOR-SNP")
	hf.SNP_FS_Max, _ = cmd.Flags().GetFloat64("min-FS-SNP")
	hf.SNP_MQ_Min, _ = cmd.Flags().GetFloat64("min-MQ-SNP")
	hf.SNP_MQRankSum_Min, _ = cmd.Flags().GetFloat64("min-MQRankSum-SNP")
	hf.SNP_ReadPosRankSum_Min, _ = cmd.Flags().GetFloat64("min-ReadPosRankSum-SNP")
	hf.INDEL_QD_Min, _ = cmd.Flags().GetFloat64("min-QD-INDEL")
	hf.INDEL_QUAL_Min, _ = cmd.Flags().GetFloat64("min-QUAL-INDEL")
	hf.INDEL_FS_Max, _ = cmd.Flags().GetFloat64("max-FS-INDEL")
	hf.INDEL_ReadPosRankSum_Min, _ = cmd.Flags().GetFloat64("min-ReadPosRankSum-INDEL")
	hf.INDEL_SOR_Max, _ = cmd.Flags().GetFloat64("max-SOR-INDEL")
	return cfg, hf, nil
}

type sampleSelections struct {
	HighParent string
	LowParent  string
	OneParent  string
	HighBulk   string
	LowBulk    string
	OneBulk    string
}

func splitCSV(raw string) []string {
	var out []string
	for _, part := range strings.Split(raw, ",") {
		part = strings.TrimSpace(part)
		if part != "" && !strings.EqualFold(part, "none") {
			out = append(out, part)
		}
	}
	return out
}

func parseIntList(raw string, want int) ([]int, error) {
	parts := splitCSV(raw)
	if len(parts) == 0 {
		return nil, fmt.Errorf("no values")
	}
	out := make([]int, want)
	for i := range out {
		src := parts[0]
		if i < len(parts) {
			src = parts[i]
		}
		v, err := strconv.Atoi(src)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

func promptSamplesFromVCF(vcf string) (sampleSelections, error) {
	samples, err := readVCFSamples(vcf)
	if err != nil {
		return sampleSelections{}, err
	}
	if len(samples) == 0 {
		return sampleSelections{}, fmt.Errorf("no sample columns found in VCF header")
	}
	color.Cyan("No parents or bulks were passed. Select samples from the VCF header.")
	reader := bufio.NewReader(os.Stdin)
	used := map[string]string{}
	var s sampleSelections
	if s.HighParent, err = promptSample(reader, samples, "resistant/high parent", used); err != nil {
		return s, err
	}
	if s.LowParent, err = promptSample(reader, samples, "susceptible/low parent", used); err != nil {
		return s, err
	}
	if s.HighBulk, err = promptSample(reader, samples, "high bulk", used); err != nil {
		return s, err
	}
	if s.LowBulk, err = promptSample(reader, samples, "low bulk", used); err != nil {
		return s, err
	}
	if s.HighParent != "" && s.LowParent == "" {
		s.OneParent, s.HighParent = s.HighParent, ""
	}
	if s.LowParent != "" && s.HighParent == "" {
		s.OneParent, s.LowParent = s.LowParent, ""
	}
	if s.HighBulk != "" && s.LowBulk == "" {
		s.OneBulk, s.HighBulk = s.HighBulk, ""
	}
	if s.LowBulk != "" && s.HighBulk == "" {
		s.OneBulk, s.LowBulk = s.LowBulk, ""
	}
	if s.HighBulk == "" && s.LowBulk == "" && s.OneBulk == "" {
		return s, fmt.Errorf("at least one bulk sample is required")
	}
	return s, nil
}

func promptSample(reader *bufio.Reader, samples []string, role string, used map[string]string) (string, error) {
	for {
		color.White("")
		color.White("Select %s:", role)
		color.White("  0) None")
		for i, sample := range samples {
			if prev, ok := used[sample]; ok {
				color.White("  %d) %s (already selected as %s)", i+1, sample, prev)
			} else {
				color.White("  %d) %s", i+1, sample)
			}
		}
		color.White("Enter number or sample name: ")
		raw, err := reader.ReadString('\n')
		if err != nil && err != io.EOF {
			return "", err
		}
		answer := strings.TrimSpace(raw)
		if answer == "" || answer == "0" || strings.EqualFold(answer, "none") {
			return "", nil
		}
		if idx, err := strconv.Atoi(answer); err == nil {
			if idx < 0 || idx > len(samples) {
				color.Red("Please choose a number between 0 and %d", len(samples))
				continue
			}
			if idx == 0 {
				return "", nil
			}
			answer = samples[idx-1]
		}
		if !contains(samples, answer) {
			color.Red("Sample %q is not present in the VCF header", answer)
			continue
		}
		if prev, ok := used[answer]; ok {
			color.Red("Sample %q is already selected as %s", answer, prev)
			continue
		}
		used[answer] = role
		return answer, nil
	}
}

func contains(values []string, value string) bool {
	for _, item := range values {
		if item == value {
			return true
		}
	}
	return false
}

func readVCFSamples(path string) ([]string, error) {
	r, closeFn, err := openMaybeGzip(path)
	if err != nil {
		return nil, err
	}
	defer closeFn()
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 1024), 1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#CHROM") {
			fields := strings.Split(line, "\t")
			if len(fields) <= 9 {
				return nil, nil
			}
			return fields[9:], nil
		}
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return nil, fmt.Errorf("VCF #CHROM header was not found")
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

func init() {
	rootCmd.Flags().StringP("variant", "V", "", "Variant File")
	rootCmd.Flags().StringP("parents", "P", "", "parent names: high,low")
	rootCmd.Flags().StringP("bulks", "B", "", "bulk names: high,low")
	rootCmd.Flags().StringP("parents-depth", "p", "5,5", "Parents minimum depth")
	rootCmd.Flags().StringP("bulks-depth", "b", "40,40", "Bulk minimum depth")
	rootCmd.Flags().StringP("bulk-sizes", "S", "20,20", "Bulk sizes")
	rootCmd.Flags().Int64P("window-size", "w", 2000000, "Window size")
	rootCmd.Flags().Int64P("step-size", "s", 100000, "Step size")
	rootCmd.Flags().StringP("population", "m", "F2", "Population type: F2, F3, BC, RIL")
	rootCmd.Flags().Int("rep", 1000, "Number of simulations")
	rootCmd.Flags().Float64Slice("alpha", []float64{0.05, 0.01}, "Significance levels")
	rootCmd.Flags().Int64("min-qtl-length", 100000, "Minimum QTL length")
	rootCmd.Flags().Int64("merge-distance", 500000, "Merge distance for QTLs")
	rootCmd.Flags().StringP("out", "o", ".", "Output directory")
	rootCmd.Flags().Float64("min-QD-SNP", 2.0, "SNP QD minimum")
	rootCmd.Flags().Float64("min-QUAL-SNP", 30.0, "SNP QUAL minimum")
	rootCmd.Flags().Float64("min-SOR-SNP", 3.0, "SNP SOR maximum")
	rootCmd.Flags().Float64("min-FS-SNP", 60.0, "SNP FS maximum")
	rootCmd.Flags().Float64("min-MQ-SNP", 40.0, "SNP MQ minimum")
	rootCmd.Flags().Float64("min-MQRankSum-SNP", -12.5, "SNP MQRankSum minimum")
	rootCmd.Flags().Float64("min-ReadPosRankSum-SNP", -8.0, "SNP ReadPosRankSum minimum")
	rootCmd.Flags().Float64("min-QD-INDEL", 2.0, "INDEL QD minimum")
	rootCmd.Flags().Float64("min-QUAL-INDEL", 30.0, "INDEL QUAL minimum")
	rootCmd.Flags().Float64("max-FS-INDEL", 200.0, "INDEL FS maximum")
	rootCmd.Flags().Float64("min-ReadPosRankSum-INDEL", -20.0, "INDEL ReadPosRankSum minimum")
	rootCmd.Flags().Float64("max-SOR-INDEL", 10.0, "INDEL SOR maximum")
	rootCmd.Flags().String("snpEffDB", "", "snpEff database")
	rootCmd.Flags().String("gff", "", "GFF3 file path")
	rootCmd.Flags().String("cds", "", "CDS FASTA path")
	rootCmd.Flags().StringP("reference", "r", "", "reference FASTA path")
	rootCmd.Flags().String("protein", "", "protein FASTA path")
	rootCmd.Flags().String("gene-descriptions", "", "gene descriptions path")
	rootCmd.Flags().String("prg", "", "PRG BLAST path")
}
