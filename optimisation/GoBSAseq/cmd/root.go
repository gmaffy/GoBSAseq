/*
Copyright © 2026 NAME HERE mafireyi@gmail.com
*/
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
	"github.com/gmaffy/GoBSAseq/run"
	"github.com/gmaffy/GoBSAseq/utils"

	"github.com/spf13/cobra"
)

var rootCmd = &cobra.Command{
	Use:   "GoBSAseq [variant.vcf]",
	Short: "Pipeline for BSAseq analysis implemented in Go",
	Long:  `Pipeline for BSAseq analysis implemented in Go.`,

	Run: func(cmd *cobra.Command, args []string) {
		variant, _ := cmd.Flags().GetString("variant")
		if strings.TrimSpace(variant) == "" {
			if len(args) == 1 {
				variant = args[0]
			} else {
				cmd.Help()
				return
			}
		}

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

		minQD_SNP, _ := cmd.Flags().GetFloat64("min-QD-SNP")
		minQUAL_SNP, _ := cmd.Flags().GetFloat64("min-QUAL-SNP")
		minSOR_SNP, _ := cmd.Flags().GetFloat64("min-SOR-SNP")
		minFS_SNP, _ := cmd.Flags().GetFloat64("min-FS-SNP")
		minMQ_SNP, _ := cmd.Flags().GetFloat64("min-MQ-SNP")
		minMQRank, _ := cmd.Flags().GetFloat64("min-MQRankSum-SNP")
		minReadPosRank, _ := cmd.Flags().GetFloat64("min-ReadPosRankSum-SNP")

		minQD_INDEL, _ := cmd.Flags().GetFloat64("min-QD-INDEL")
		minQUAL_INDEL, _ := cmd.Flags().GetFloat64("min-QUAL-INDEL")
		maxFS_INDEL, _ := cmd.Flags().GetFloat64("max-FS-INDEL")
		minReadPosRank_INDEL, _ := cmd.Flags().GetFloat64("min-ReadPosRankSum-INDEL")
		maxSOR_INDEL, _ := cmd.Flags().GetFloat64("max-SOR-INDEL")

		snpEffDB, _ := cmd.Flags().GetString("snpEffDB")
		gff, _ := cmd.Flags().GetString("gff")
		cds, _ := cmd.Flags().GetString("cds")
		ref, _ := cmd.Flags().GetString("reference")
		protein, _ := cmd.Flags().GetString("protein")
		geneDescriptions, _ := cmd.Flags().GetString("gene-descriptions")
		prg, _ := cmd.Flags().GetString("prg")

		var promptedSelections *sampleSelections
		promptSamples := !cmd.Flags().Changed("parents") && !cmd.Flags().Changed("bulks") && strings.TrimSpace(parents) == "" && strings.TrimSpace(bulks) == ""
		if promptSamples {
			selections, err := promptSamplesFromVCF(variant)
			if err != nil {
				color.Red("%v", err)
				return
			}
			promptedSelections = &selections
		}

		parentNamesLst := splitCSV(parents)
		bulkNamesLst := splitCSV(bulks)
		bulksDepthLst := strings.Split(bulksDepth, ",")
		bulkSizesLst := strings.Split(bulkSizes, ",")
		parentsDepthLst := strings.Split(parentsDepth, ",")

		highParentName := ""
		lowParentName := ""
		oneParentName := ""

		highBulkName := ""
		lowBulkName := ""
		oneBulkName := ""

		highBulkDepth := 0
		lowBulkDepth := 0
		oneBulkDepth := 0
		highParentDepth := 0
		lowParentDepth := 0
		oneParentDepth := 0
		highBulkSize := 0
		lowBulkSize := 0
		oneBulkSize := 0
		winSize := 0

		var err error

		// ========================================== Get parents =================================================== //
		if promptedSelections != nil {
			highParentName = promptedSelections.HighParent
			lowParentName = promptedSelections.LowParent
			if highParentName != "" && lowParentName == "" {
				oneParentName = highParentName
				highParentName = ""
			} else if lowParentName != "" && highParentName == "" {
				oneParentName = lowParentName
				lowParentName = ""
			}
		} else if len(parentNamesLst) > 0 {
			if len(parentNamesLst) > 2 {
				color.Red("parents is supposed to be in the form a,b (where a and b are integers)")
				return
			} else if len(parentNamesLst) == 1 {
				oneParentName = parentNamesLst[0]

			} else {
				highParentName = parentNamesLst[0]
				lowParentName = parentNamesLst[1]
			}

		}

		// ========================================== Get Bulk Names =============================================== //
		if promptedSelections != nil {
			highBulkName = promptedSelections.HighBulk
			lowBulkName = promptedSelections.LowBulk
			if highBulkName != "" && lowBulkName == "" {
				oneBulkName = highBulkName
				highBulkName = ""
			} else if lowBulkName != "" && highBulkName == "" {
				oneBulkName = lowBulkName
				lowBulkName = ""
			}
		} else if len(bulkNamesLst) > 0 {
			if len(bulkNamesLst) > 2 {
				color.Red("bulks is supposed to be in the form a,b (where a and b are integers)")
				return
			} else if len(bulkNamesLst) == 1 {
				oneBulkName = bulkNamesLst[0]

			} else {
				highBulkName = bulkNamesLst[0]
				lowBulkName = bulkNamesLst[1]
			}
		}

		// ========================================== Get bulk depths =============================================== //
		if len(bulksDepth) > 0 {
			if len(bulksDepthLst) > 2 {
				color.Red("bulksDepth is supposed to be in the form a,b (where a and b are integers)")
				return
			} else if len(bulksDepthLst) == 1 {
				oneBulkDepth, err = strconv.Atoi(bulksDepthLst[0])
				if err != nil {
					color.Red("bulksDepth is supposed to be in the form a,b (where a and b are integers)")
					return
				}

			} else {
				highBulkDepth, err = strconv.Atoi(bulksDepthLst[0])
				if err != nil {
					color.Red("bulksDepth is supposed to be in the form a,b (where a and b are integers)")
					return
				}
				lowBulkDepth, err = strconv.Atoi(bulksDepthLst[1])
				if err != nil {
					color.Red("bulksDepth is supposed to be in the form a,b (where a and b are integers)")
					return
				}
			}
		}

		// =========================================== Get Parent Depths============================================= //

		if len(parentsDepthLst) > 0 {
			if len(parentsDepthLst) > 2 {
				color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
				return
			} else if len(parentsDepthLst) == 1 {
				oneParentDepth, err = strconv.Atoi(parentsDepthLst[0])
				if err != nil {
					color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
					return
				}
			} else {
				highParentDepth, err = strconv.Atoi(parentsDepthLst[0])
				if err != nil {
					color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
					return
				}
				lowParentDepth, err = strconv.Atoi(parentsDepthLst[1])
				if err != nil {
					color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
				}

			}

		}

		// ========================================== Get Bulk Sizes ================================================ //
		if len(bulkSizesLst) > 0 {
			if len(bulkSizesLst) > 2 {
				color.Red("bulkSizes is supposed to be in the form a,b (where a and b are integers)")
				return
			} else if len(bulkSizesLst) == 1 {
				oneBulkSize, err = strconv.Atoi(bulkSizesLst[0])

			} else {
				highBulkSize, err = strconv.Atoi(bulkSizesLst[0])
				if err != nil {
					color.Red("bulkSizes is supposed to be in the form a,b (where a and b are integers)")
					return
				}
				lowBulkSize, err = strconv.Atoi(bulkSizesLst[1])
				if err != nil {
					color.Red("bulkSizes is supposed to be in the form a,b (where a and b are integers)")
					return
				}
			}
		}

		// ============================================== window size ================================================//
		winSize, err = strconv.Atoi(fmt.Sprintf("%d", windowSize))
		if err != nil {
			color.Red("windowSize is supposed to be an integer")
			return
		}

		// ============================================= Step size ================================================ //
		step, err := strconv.Atoi(fmt.Sprintf("%d", stepSize))
		if err != nil {
			color.Red("stepSize is supposed to be an integer")
			return
		}

		// ==================================== Hard Filter config ================================================== //

		qd_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minQD_SNP), 64)
		if err != nil {
			color.Red("min-QD-SNP is supposed to be a float")
			return
		}
		qual_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minQUAL_SNP), 64)
		if err != nil {
			color.Red("min-QUAL-SNP is supposed to be a float")
			return
		}
		sor_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minSOR_SNP), 64)
		if err != nil {
			color.Red("min-SOR-SNP is supposed to be a float")
			return
		}
		fs_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minFS_SNP), 64)
		if err != nil {
			color.Red("min-FS-SNP is supposed to be a float")
			return
		}
		mq_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minMQ_SNP), 64)
		if err != nil {
			color.Red("min-MQ-SNP is supposed to be a float")
			return
		}
		mqRank_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minMQRank), 64)
		if err != nil {
			color.Red("min-MQRankSum-SNP is supposed to be a float")
			return
		}
		readPosRank_snp, err := strconv.ParseFloat(fmt.Sprintf("%f", minReadPosRank), 64)
		if err != nil {
			color.Red("min-ReadPosRankSum-SNP is supposed to be a float")
			return
		}

		qd_indel, err := strconv.ParseFloat(fmt.Sprintf("%f", minQD_INDEL), 64)
		if err != nil {
			color.Red("min-QD-INDEL is supposed to be a float")
			return
		}
		qual_indel, err := strconv.ParseFloat(fmt.Sprintf("%f", minQUAL_INDEL), 64)
		if err != nil {
			color.Red("min-QUAL-INDEL is supposed to be a float")
			return
		}
		fs_indel, err := strconv.ParseFloat(fmt.Sprintf("%f", maxFS_INDEL), 64)
		if err != nil {
			color.Red("max-FS-INDEL is supposed to be a float")
			return
		}
		readPosRank_indel, err := strconv.ParseFloat(fmt.Sprintf("%f", minReadPosRank_INDEL), 64)
		if err != nil {
			color.Red("min-ReadPosRankSum-INDEL is supposed to be a float")
			return
		}

		sor_indel, err := strconv.ParseFloat(fmt.Sprintf("%f", maxSOR_INDEL), 64)
		if err != nil {
			color.Red("max-SOR-INDEL is supposed to be a float")
			return
		}

		hfConfig := utils.HardFilterConfig{
			SNP_QD_Min:               qd_snp,
			SNP_QUAL_Min:             qual_snp,
			SNP_SOR_Max:              sor_snp,
			SNP_FS_Max:               fs_snp,
			SNP_MQ_Min:               mq_snp,
			SNP_MQRankSum_Min:        mqRank_snp,
			SNP_ReadPosRankSum_Min:   readPosRank_snp,
			INDEL_QD_Min:             qd_indel,
			INDEL_QUAL_Min:           qual_indel,
			INDEL_FS_Max:             fs_indel,
			INDEL_ReadPosRankSum_Min: readPosRank_indel,
			INDEL_SOR_Max:            sor_indel,
		}

		a_config := utils.AnalysisConfig{
			VCF:           variant,
			Population:    population,
			WindowSize:    winSize,
			StepSize:      step,
			Rep:           rep,
			Alphas:        alphas,
			MinQTLWidth:   minQTL,
			MergeDistance: mergeDist,
			OutputDir:     outputDir,
			//HighParentIdx:    -1,
			HighParentName:  highParentName,
			HighParentDepth: highParentDepth,
			//LowParentIdx:    -1,
			LowParentName:  lowParentName,
			LowParentDepth: lowParentDepth,
			//OneParentIdx:    -1,
			OneParentName:  oneParentName,
			OneParentDepth: oneParentDepth,
			//HighBulkIdx:    -1,
			HighBulkName:  highBulkName,
			HighBulkDepth: highBulkDepth,
			HighBulkSize:  highBulkSize,
			//LowBulkIdx:    -1,
			LowBulkName:  lowBulkName,
			LowBulkDepth: lowBulkDepth,
			LowBulkSize:  lowBulkSize,
			//OneBulkIdx:    -1,
			OneBulkName:  oneBulkName,
			OneBulkDepth: oneBulkDepth,
			OneBulkSize:  oneBulkSize,

			SnpEffDB: snpEffDB,
			Ref:      ref,
			Protein:  protein,
			Gff:      gff,
			Cds:      cds,
			GeneDesc: geneDescriptions,
			Prg:      prg,
		}

		err = run.Run(a_config, hfConfig)
		if err != nil {
			return
		}
		//fmt.Printf("variant: %s\n", variant)
	},
}

func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}

type sampleSelections struct {
	HighParent string
	LowParent  string
	HighBulk   string
	LowBulk    string
}

func splitCSV(raw string) []string {
	var out []string
	for _, part := range strings.Split(raw, ",") {
		part = strings.TrimSpace(part)
		if part == "" || strings.EqualFold(part, "none") {
			continue
		}
		out = append(out, part)
	}
	return out
}

func promptSamplesFromVCF(vcf string) (sampleSelections, error) {
	samples, err := readVCFSamples(vcf)
	if err != nil {
		return sampleSelections{}, err
	}
	if len(samples) == 0 {
		return sampleSelections{}, fmt.Errorf("no VCF samples found in header")
	}

	color.Cyan("No parents or bulks were passed. Select samples from the VCF header.")
	reader := bufio.NewReader(os.Stdin)
	var s sampleSelections
	used := map[string]string{}
	if s.HighParent, err = promptSample(reader, samples, "resistant/high parent", used); err != nil {
		return sampleSelections{}, err
	}
	if s.LowParent, err = promptSample(reader, samples, "susceptible/low parent", used); err != nil {
		return sampleSelections{}, err
	}
	if s.HighBulk, err = promptSample(reader, samples, "high bulk", used); err != nil {
		return sampleSelections{}, err
	}
	if s.LowBulk, err = promptSample(reader, samples, "low bulk", used); err != nil {
		return sampleSelections{}, err
	}
	if s.HighBulk == "" && s.LowBulk == "" {
		return sampleSelections{}, fmt.Errorf("at least one bulk sample is required")
	}
	return s, nil
}

func promptSample(reader *bufio.Reader, samples []string, role string, used map[string]string) (string, error) {
	for {
		color.White("")
		color.White("Select %s:", role)
		color.White("  0) None")
		for i, sample := range samples {
			if prevRole, ok := used[sample]; ok {
				color.White("  %d) %s (already selected as %s)", i+1, sample, prevRole)
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
		if answer == "" || strings.EqualFold(answer, "none") || answer == "0" {
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
		if !containsSample(samples, answer) {
			color.Red("Sample %q is not present in the VCF header", answer)
			continue
		}
		if prevRole, ok := used[answer]; ok {
			color.Red("Sample %q is already selected as %s", answer, prevRole)
			continue
		}
		used[answer] = role
		return answer, nil
	}
}

func containsSample(samples []string, sample string) bool {
	for _, item := range samples {
		if item == sample {
			return true
		}
	}
	return false
}

func readVCFSamples(path string) ([]string, error) {
	r, closeFn, err := openVCFHeaderReader(path)
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
	return nil, fmt.Errorf("VCF header line beginning with #CHROM was not found")
}

func openVCFHeaderReader(path string) (io.Reader, func(), error) {
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
		return gz, func() {
			_ = gz.Close()
			_ = f.Close()
		}, nil
	}
	return f, func() { _ = f.Close() }, nil
}

func init() {
	rootCmd.Flags().StringP("variant", "V", "", "Variant File")
	rootCmd.Flags().StringP("parents", "P", "", "parent names")
	rootCmd.Flags().StringP("bulks", "B", "", "bulk names")
	rootCmd.Flags().StringP("parents-depth", "p", "5,5", "Parents Min Depth")
	rootCmd.Flags().StringP("bulks-depth", "b", "40,40", "Low Parent Min Depth")
	rootCmd.Flags().StringP("bulk-sizes", "S", "20,20", "High Bulk Min Depth")
	rootCmd.Flags().Int64P("window-size", "w", 2000000, "Window Size")
	rootCmd.Flags().Int64P("step-size", "s", 100000, "Step Size")
	rootCmd.Flags().StringP("population", "m", "F2", "Population type (F2, F3, BC, RIL)")
	rootCmd.Flags().Bool("recurrent", false, "BCAltIsRecurrent: if true, alt allele is recurrent in BC")
	rootCmd.Flags().Int("rep", 1000, "Number of simulations")
	rootCmd.Flags().Float64Slice("alpha", []float64{0.05, 0.01}, "Significance levels (comma-separated)")
	rootCmd.Flags().Int64("min-qtl-length", 100000, "Minimum QTL length")
	rootCmd.Flags().Int64("merge-distance", 500000, "Merge distance for QTLs")
	rootCmd.Flags().StringP("out", "o", ".", "Output directory")

	rootCmd.Flags().Float64("min-QD-SNP", 2.0, "QualByDepth SNPs") // SNP_QD_Min             float64 // default 2.0   – QualByDepth
	rootCmd.Flags().Float64("min-QUAL-SNP", 30.0, "Variant quality SNPs")
	rootCmd.Flags().Float64("min-SOR-SNP", 3.0, "StrandOddsRatio SNPs")
	rootCmd.Flags().Float64("min-FS-SNP", 60.0, "FisherStrand SNPs")
	rootCmd.Flags().Float64("min-MQ-SNP", 40.0, "RMSMappingQuality SNPs")
	rootCmd.Flags().Float64("min-MQRankSum-SNP", -12.5, "MappingQualityRank SNPs")
	rootCmd.Flags().Float64("min-ReadPosRankSum-SNP", -8.0, "ReadPosRank SNPs")

	rootCmd.Flags().Float64("min-QD-INDEL", 2.0, "QualByDepth INDELs")               // INDEL_QD_Min             float64 // default 2.0   – QualByDepth
	rootCmd.Flags().Float64("min-QUAL-INDEL", 30.0, "Variant quality INDELs")        //INDEL_QUAL_Min           float64 // default 30.0  – variant quality score
	rootCmd.Flags().Float64("max-FS-INDEL", 200.0, "FisherStrand INDELs")            //INDEL_FS_Max             float64 // default 200.0 – FisherStrand
	rootCmd.Flags().Float64("min-ReadPosRankSum-INDEL", -20.0, "ReadPosRank INDELs") //INDEL_ReadPosRankSum_Min float64 // default -20.0 – ReadPosRankSumTest
	rootCmd.Flags().Float64("max-SOR-INDEL", 10.0, "StrandOddsRatio INDELs")

	rootCmd.Flags().String("snpEffDB", "", "snpEff database")
	rootCmd.Flags().String("gff", "", "gff3 file path")
	rootCmd.Flags().String("cds", "", "cds file path")
	rootCmd.Flags().StringP("reference", "r", "", "cds fasta path")
	rootCmd.Flags().String("protein", "", "protein fasta path")
	rootCmd.Flags().String("gene-descriptions", "", "gene descriptions file path")
	rootCmd.Flags().String("prg", "", "prg blast file path ")

}
