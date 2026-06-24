/*
Copyright © 2026 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/run"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/spf13/cobra"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "GoBSAseq",
	Short: "Pipeline for BSAseq analysis implemented in Go",
	Long:  `Pipeline for BSAseq analysis implemented in Go.`,
	// Uncomment the following line if your bare application
	// has an action associated with it:
	Run: func(cmd *cobra.Command, args []string) {

		variant, _ := cmd.Flags().GetString("variant")
		parents, _ := cmd.Flags().GetString("parents")
		bulks, _ := cmd.Flags().GetString("bulks")
		parentsBams, _ := cmd.Flags().GetString("parents-bams")
		bulksBams, _ := cmd.Flags().GetString("bulks-bams")
		hpReads, _ := cmd.Flags().GetString("hp-reads")
		lpReads, _ := cmd.Flags().GetString("lp-reads")
		hbReads, _ := cmd.Flags().GetString("hb-reads")
		lbReads, _ := cmd.Flags().GetString("lb-reads")

		parentsDepth, _ := cmd.Flags().GetString("parents-depth")
		bulksDepth, _ := cmd.Flags().GetString("bulks-depth")
		bulkSizes, _ := cmd.Flags().GetString("bulk-sizes")
		windowSize, _ := cmd.Flags().GetInt64("window-size")
		stepSize, _ := cmd.Flags().GetInt64("step-size")
		population, _ := cmd.Flags().GetString("population")
		rep, _ := cmd.Flags().GetInt("rep")
		brmAlpha, _ := cmd.Flags().GetFloat64("brm-alpha")
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

		merger, _ := cmd.Flags().GetString("merger")
		caller, _ := cmd.Flags().GetString("caller")
		noMerging, _ := cmd.Flags().GetBool("no-merging")
		deepVariantVersion, _ := cmd.Flags().GetString("deepvariant-version")
		modelType, _ := cmd.Flags().GetString("model-type")
		verbose, _ := cmd.Flags().GetBool("verbose")
		threads, _ := cmd.Flags().GetInt("threads")

		parentNamesLst := strings.Split(parents, ",")
		bulkNamesLst := strings.Split(bulks, ",")           //splitArg(bulks)
		bulksDepthLst := strings.Split(bulksDepth, ",")     //splitArg(bulksDepth)
		bulkSizesLst := strings.Split(bulkSizes, ",")       //splitArg(bulkSizes)
		parentsDepthLst := strings.Split(parentsDepth, ",") //splitArg(parentsDepth)

		parentBamsLst := strings.Split(parentsBams, ",")
		bulkBamsLst := strings.Split(bulksBams, ",")
		hpReadsLst := strings.Split(hpReads, ",")
		lpReadsLst := strings.Split(lpReads, ",")
		hbReadsLst := strings.Split(hbReads, ",")
		lbReadsLst := strings.Split(lbReads, ",")

		highParentName := ""
		lowParentName := ""

		hpBam := ""
		lpBam := ""
		hbBam := ""
		lbBam := ""

		hpFwdReads := ""
		hpRevReads := ""
		lpFwdReads := ""
		lpRevReads := ""
		hbFwdReads := ""
		hbRevReads := ""
		lbFwdReads := ""
		lbRevReads := ""

		highBulkName := ""
		lowBulkName := ""

		highBulkDepth := 0
		lowBulkDepth := 0

		highParentDepth := 0
		lowParentDepth := 0

		highBulkSize := 0
		lowBulkSize := 0

		winSize := 0

		var err error

		// ========================================== Get parents =================================================== //

		//------------------------------------------- Parent Names -------------------------------------------------- //
		if len(parentNamesLst) > 0 {
			if len(parentNamesLst) != 2 {
				color.Red("parents is supposed to be in the form a,b (where a and b are parent names)")
				return
			}

			highParentName = parentNamesLst[0]
			lowParentName = parentNamesLst[1]

		}

		// ---------------------------------------- Parent Bams ----------------------------------------------------- //
		if len(parentBamsLst) > 0 {
			if len(parentBamsLst) != 2 {
				color.Red("parentsBams is supposed to be in the form a,b (where a and b are parent bam files)")
				return
			}
			hpBam = parentBamsLst[0]
			lpBam = parentBamsLst[1]

		}

		// --------------------------------------- Parent Reads ----------------------------------------------------- //
		if len(hpReadsLst) > 0 {
			if len(hpReadsLst) != 2 {
				color.Red("hpReads is supposed to be in the form a,b (where a and b are parent reads)")
				return
			}
			hpFwdReads = hpReadsLst[0]
			hpRevReads = hpReadsLst[1]
		}

		if len(lpReadsLst) > 0 {
			if len(lpReadsLst) != 2 {
				color.Red("lpReads is supposed to be in the form a,b (where a and b are parent reads)")
				return
			}
			lpFwdReads = lpReadsLst[0]
			lpRevReads = lpReadsLst[1]
		}

		// ========================================== Get Bulk Names =============================================== //
		if len(bulkNamesLst) > 0 {
			if len(bulkNamesLst) != 2 {
				color.Red("bulks is supposed to be in the form a,b (where a and b are bulk names)")
				return
			}
			highBulkName = bulkNamesLst[0]
			lowBulkName = bulkNamesLst[1]

		}

		// ---------------------------------------- Bulk Bams ----------------------------------------------------- //
		if len(bulkBamsLst) > 0 {
			if len(bulkBamsLst) != 2 {
				color.Red("bulkBams is supposed to be in the form a,b (where a and b are bam files)")
				return
			}
			hbBam = bulkBamsLst[0]
			lbBam = bulkBamsLst[1]
		}

		// ---------------------------------------- Bulk Reads ----------------------------------------------------- //
		if len(hbReadsLst) > 0 {
			if len(hbReadsLst) != 2 {
				color.Red("hbReads is supposed to be in the form fwdReads,revReads")
				return
			}
			hbFwdReads = hbReadsLst[0]
			hbRevReads = hbReadsLst[1]
		}

		if len(lbReadsLst) > 0 {
			if len(lbReadsLst) != 2 {
				color.Red("lbReads is supposed to be in the form fwdReads,revReads")
				return
			}
			lbFwdReads = lbReadsLst[0]
			lbRevReads = lbReadsLst[1]
		}

		// ========================================== Get bulk depths =============================================== //
		if len(bulksDepth) > 0 {
			if len(bulksDepthLst) != 2 {
				color.Red("bulksDepth is supposed to be in the form a,b (where a and b are integers)")
				return
			}

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

		// =========================================== Get Parent Depths============================================= //

		if len(parentsDepthLst) > 0 {
			if len(parentsDepthLst) != 2 {
				color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
				return
			}

			highParentDepth, err = strconv.Atoi(parentsDepthLst[0])
			if err != nil {
				color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
				return
			}
			lowParentDepth, err = strconv.Atoi(parentsDepthLst[1])
			if err != nil {
				color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
				return
			}

		}

		// ========================================== Get Bulk Sizes ================================================ //
		if len(bulkSizesLst) > 0 {
			if len(bulkSizesLst) != 2 {
				color.Red("bulkSizes is supposed to be in the form a,b (where a and b are integers)")
				return
			}
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

		// ============================================== window size ================================================//
		winSize = int(windowSize)
		step := int(stepSize)

		qd_snp := minQD_SNP
		qual_snp := minQUAL_SNP
		sor_snp := minSOR_SNP
		fs_snp := minFS_SNP
		mq_snp := minMQ_SNP
		mqRank_snp := minMQRank
		readPosRank_snp := minReadPosRank

		qd_indel := minQD_INDEL
		qual_indel := minQUAL_INDEL
		fs_indel := maxFS_INDEL
		readPosRank_indel := minReadPosRank_INDEL
		sor_indel := maxSOR_INDEL

		// --------------------------------- Output dir --------------------------------------------------//

		resultsDir, err := utils.CreateResultsDir(outputDir)
		if err != nil {
			fmt.Println("Error creating results directory:", err)
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
			VCF: variant,

			HighParFwdReads: hpFwdReads,
			HighParRevReads: hpRevReads,
			LowParFwdReads:  lpFwdReads,
			LowParRevReads:  lpRevReads,

			HighBulkFwdReads: hbFwdReads,
			HighBulkRevReads: hbRevReads,
			LowBulkFwdReads:  lbFwdReads,
			LowBulkRevReads:  lbRevReads,

			HighParBam:  hpBam,
			LowParBam:   lpBam,
			HighBulkBam: hbBam,
			LowBulkBam:  lbBam,

			Population:    population,
			WindowSize:    winSize,
			StepSize:      step,
			Rep:           rep,
			BrmAlpha:      brmAlpha,
			MinQTLWidth:   minQTL,
			MergeDistance: mergeDist,
			OutputDir:     resultsDir,

			HighParentName:  highParentName,
			HighParentDepth: highParentDepth,

			LowParentName:  lowParentName,
			LowParentDepth: lowParentDepth,

			HighBulkName:  highBulkName,
			HighBulkDepth: highBulkDepth,
			HighBulkSize:  highBulkSize,

			LowBulkName:  lowBulkName,
			LowBulkDepth: lowBulkDepth,
			LowBulkSize:  lowBulkSize,

			SnpEffDB: snpEffDB,
			Ref:      ref,
			Protein:  protein,
			Gff:      gff,
			Cds:      cds,
			GeneDesc: geneDescriptions,
			Prg:      prg,

			Caller:             caller,
			Merger:             merger,
			NoMerging:          noMerging,
			DeepVariantVersion: deepVariantVersion,
			ModelType:          modelType,
			Verbose:            verbose,
			Threads:            threads,
		}

		err = run.Run(&a_config, hfConfig)
		if err != nil {
			return
		}

	},
}

func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}

func init() {
	// ---------------------------------------- Input files ----------------------------------------------------- //
	rootCmd.Flags().StringP("variant", "V", "", "Variant File")
	rootCmd.Flags().StringP("parents", "P", ",", "parent names (hp,lp)")
	rootCmd.Flags().StringP("bulks", "B", ",", "bulk names (hb,lb)")
	rootCmd.Flags().String("parents-bams", ",", "parent bam files (hp,lp)")
	rootCmd.Flags().String("bulks-bams", ",", "bulk bam files (hb,lb)")
	rootCmd.Flags().String("hp-reads", ",", "high parent reads (fwd,rev)")
	rootCmd.Flags().String("lp-reads", ",", "low parent reads (fwd,rev)")
	rootCmd.Flags().String("hb-reads", ",", "high bulk reads (fwd,rev)")
	rootCmd.Flags().String("lb-reads", ",", "low bulk reads (fwd,rev)")

	// ================================================ Parameters ================================================== //
	// -------------------------------------------------- Inputs ---------------------------------------------------- //
	rootCmd.Flags().StringP("parents-depth", "p", "5,5", "Parents Min Depth")
	rootCmd.Flags().StringP("bulks-depth", "b", "40,40", "Low Parent Min Depth")
	rootCmd.Flags().StringP("bulk-sizes", "S", "20,20", "High Bulk Min Depth")

	// ------------------------------------------------- Smoothing -------------------------------------------------- //
	rootCmd.Flags().Int64P("window-size", "w", 2000000, "Window Size")
	rootCmd.Flags().Int64P("step-size", "s", 100000, "Step Size")

	// ------------------------------------------------- Threshold -------------------------------------------------- //
	rootCmd.Flags().Float64("brm-alpha", 0.05, "Significance level")
	rootCmd.Flags().StringP("population", "m", "F2", "Population type (F2, F3, RIL, BC1H, BC1L, BC2H, BC2L, ...)")
	rootCmd.Flags().Int("rep", 1000, "Number of simulations")

	// ----------------------------------------------- Filtering ---------------------------------------------------- //
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

	// ------------------------------------- Gene space analysis --------------------------------------------------- //

	rootCmd.Flags().String("snpEffDB", "", "snpEff database")
	rootCmd.Flags().String("gff", "", "gff3 file path")
	rootCmd.Flags().String("cds", "", "cds file path")
	rootCmd.Flags().StringP("reference", "r", "", "reference file path")
	rootCmd.Flags().String("protein", "", "protein fasta path")
	rootCmd.Flags().String("gene-descriptions", "", "gene descriptions file path")
	rootCmd.Flags().String("prg", "", "prg blast file path ")

	// ------------------------------------------- OutDir ----------------------------------------------------------- //
	rootCmd.Flags().StringP("out", "o", ".", "Output directory")

	// ------------------------------------------ Variant Calling --------------------------------------------------- //
	rootCmd.Flags().String("caller", "deepvariant", "Variant caller to use. Options: gatk or deepvariant")
	rootCmd.Flags().String("merger", "glnexus", "GVCF merger to use. Options: gatk or glnexus")
	rootCmd.Flags().Bool("no-merging", false, "do not merge gvcfs.")
	rootCmd.Flags().String("deepvariant-version", "1.10.0", "DeepVariant version")
	rootCmd.Flags().String("model-type", "WGS", "DeepVariant Model Type: WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA")
	rootCmd.Flags().BoolP("verbose", "v", false, "Verbose")
	
	rootCmd.Flags().IntP("threads", "t", 4, "Number of threads")

}
