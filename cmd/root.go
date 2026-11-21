/*
Copyright © 2026 Godwin Mafireyi <mafireyi@gmail.com>
*/
package cmd

import (
	"bufio"
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

	Run: func(cmd *cobra.Command, args []string) {

		ref := ""
		gff := ""
		cds := ""
		protein := ""
		geneDescriptions := ""
		prg := ""
		vcfFile := ""

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

		var err error
		configFilePath, _ := cmd.Flags().GetString("config")
		if configFilePath != "" {
			configFile, err := os.Open(configFilePath)
			if err != nil {
				color.Red("Error opening config file: %v", err)
				return
			}
			defer configFile.Close()

			scanner := bufio.NewScanner(configFile)
			for scanner.Scan() {
				line := strings.TrimSpace(scanner.Text())

				if line == "" {
					continue
				}

				parts := strings.SplitN(line, ":", 2)
				if len(parts) != 2 {
					continue
				}

				key := strings.TrimSpace(parts[0])
				value := strings.TrimSpace(parts[1])

				switch key {
				case "Reference":
					ref = value
				case "gff":
					gff = value
				case "gene-descriptions":
					geneDescriptions = value
				case "prg":
					prg = value
				case "HighParentReads":
					hpReads := strings.Split(value, " ")
					if len(hpReads) == 3 {
						hpFwdReads = hpReads[0]
						hpRevReads = hpReads[1]
						highParentName = hpReads[2]
					}
				case "LowParentReads":
					lpReads := strings.Split(value, " ")
					if len(lpReads) == 3 {
						lpFwdReads = lpReads[0]
						lpRevReads = lpReads[1]
						lowParentName = lpReads[2]
					}
				case "HighBulkReads":
					hbReads := strings.Split(value, " ")
					if len(hbReads) == 3 {
						hbFwdReads = hbReads[0]
						hbRevReads = hbReads[1]
						highBulkName = hbReads[2]
					}
				case "LowBulkReads":
					lbReads := strings.Split(value, " ")
					if len(lbReads) == 3 {
						lbFwdReads = lbReads[0]
						lbRevReads = lbReads[1]
						lowBulkName = lbReads[2]
					}
				case "HighParentBam":
					hpBam = value
				case "LowParentBam":
					lpBam = value
				case "HighBulkBam":
					hbBam = value
				case "LowBulkBam":
					lbBam = value
				case "VCF":
					vcfFile = value
				}
			}

		} else {
			vcfFile, _ = cmd.Flags().GetString("variant")

			// ============================================  Bams ================================================ //
			parentsBams, _ := cmd.Flags().GetString("parents-bams")
			parentBamsLst := strings.Split(parentsBams, ",")
			if len(parentBamsLst) == 2 {
				hpBam = parentBamsLst[0]
				lpBam = parentBamsLst[1]
			}

			bulksBams, _ := cmd.Flags().GetString("bulks-bams")
			bulkBamsLst := strings.Split(bulksBams, ",")
			if len(bulkBamsLst) == 2 {
				hbBam = bulkBamsLst[0]
				lbBam = bulkBamsLst[1]
			}

			// =========================================== Reads ==================================================== //

			hpReads, _ := cmd.Flags().GetString("hp-reads")
			hpReadsLst := strings.Split(hpReads, ",")
			if len(hpReadsLst) == 2 {
				hpFwdReads = hpReadsLst[0]
				hpRevReads = hpReadsLst[1]
			}

			lpReads, _ := cmd.Flags().GetString("lp-reads")
			lpReadsLst := strings.Split(lpReads, ",")
			if len(lpReadsLst) == 2 {
				lpFwdReads = lpReadsLst[0]
				lpRevReads = lpReadsLst[1]
			}

			hbReads, _ := cmd.Flags().GetString("hb-reads")
			hbReadsLst := strings.Split(hbReads, ",")
			if len(hbReadsLst) == 2 {
				hbFwdReads = hbReadsLst[0]
				hbRevReads = hbReadsLst[1]

			}

			lbReads, _ := cmd.Flags().GetString("lb-reads")
			lbReadsLst := strings.Split(lbReads, ",")
			if len(lbReadsLst) == 2 {
				lbFwdReads = lbReadsLst[0]
				lbRevReads = lbReadsLst[1]
			}

			// =============================================== Names ================================================ //

			parents, _ := cmd.Flags().GetString("parents")
			parentNamesLst := strings.Split(parents, ",")
			if len(parentNamesLst) == 2 {
				highParentName = parentNamesLst[0]
				lowParentName = parentNamesLst[1]
			}

			bulks, _ := cmd.Flags().GetString("bulks")
			bulkNamesLst := strings.Split(bulks, ",")
			if len(bulkNamesLst) == 2 {
				highBulkName = bulkNamesLst[0]
				lowBulkName = bulkNamesLst[1]
			}

			// ======================================= Gene Space Analysis ========================================== //

			gff, _ = cmd.Flags().GetString("gff")
			cds, _ = cmd.Flags().GetString("cds")
			ref, _ = cmd.Flags().GetString("reference")
			protein, _ = cmd.Flags().GetString("protein")
			geneDescriptions, _ = cmd.Flags().GetString("gene-descriptions")
			prg, _ = cmd.Flags().GetString("prg")

		}

		parentsDepth, _ := cmd.Flags().GetString("parents-depth")
		bulksDepth, _ := cmd.Flags().GetString("bulks-depth")
		bulksMaxDepth, _ := cmd.Flags().GetString("bulks-max-depth")
		parentsMaxDepth, _ := cmd.Flags().GetString("parents-max-depth")
		minGQ, _ := cmd.Flags().GetInt("min-GQ")
		noSplitMultiallelic, _ := cmd.Flags().GetBool("no-split-multiallelic")
		minPooledMAF, _ := cmd.Flags().GetFloat64("min-pooled-maf")
		requireNonParent, _ := cmd.Flags().GetBool("require-nonparent-allele")
		oneBulkP0Mode, _ := cmd.Flags().GetString("one-bulk-p0")
		oneBulkP0Mode = strings.ToLower(strings.TrimSpace(oneBulkP0Mode))
		if oneBulkP0Mode != "empirical" && oneBulkP0Mode != "theoretical" {
			color.Red("Error: --one-bulk-p0 must be empirical or theoretical, got %q", oneBulkP0Mode)
			os.Exit(1)
		}
		cmPerMb, _ := cmd.Flags().GetFloat64("cm-per-mb")
		nullReps, _ := cmd.Flags().GetInt("null-reps")
		if cmPerMb < 0 || nullReps < 20 {
			color.Red("Error: --cm-per-mb must be ≥ 0 and --null-reps ≥ 20")
			os.Exit(1)
		}
		refParent, _ := cmd.Flags().GetString("ref-parent")
		refParent = strings.ToLower(strings.TrimSpace(refParent))
		if refParent != "high" && refParent != "low" && refParent != "none" {
			color.Red("Error: --ref-parent must be high, low or none, got %q", refParent)
			os.Exit(1)
		}
		bulkSizes, _ := cmd.Flags().GetString("bulk-sizes")
		windowSize, _ := cmd.Flags().GetInt64("window-size")
		stepSize, _ := cmd.Flags().GetInt64("step-size")
		population, _ := cmd.Flags().GetString("population")
		if _, err := utils.ParsePopulation(population); err != nil {
			color.Red("Error: --population: %v", err)
			os.Exit(1)
		}
		brmAlpha, _ := cmd.Flags().GetFloat64("brm-alpha")
		mergeDistance, _ := cmd.Flags().GetInt64("merge-distance")
		minQTLMarkers, _ := cmd.Flags().GetInt("min-qtl-markers")
		geneSpaceSpan, _ := cmd.Flags().GetBool("gene-space-span")
		outputDir, _ := cmd.Flags().GetString("out")

		lightFiltering, _ := cmd.Flags().GetBool("light-filtering")
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

		merger, _ := cmd.Flags().GetString("merger")
		caller, _ := cmd.Flags().GetString("caller")
		noMerging, _ := cmd.Flags().GetBool("no-merging")
		deepVariantVersion, _ := cmd.Flags().GetString("deepvariant-version")
		modelType, _ := cmd.Flags().GetString("model-type")
		verbose, _ := cmd.Flags().GetBool("verbose")
		threads, _ := cmd.Flags().GetInt("threads")

		bulksDepthLst := strings.Split(bulksDepth, ",")     //splitArg(bulksDepth)
		bulkSizesLst := strings.Split(bulkSizes, ",")       //splitArg(bulkSizes)
		parentsDepthLst := strings.Split(parentsDepth, ",") //splitArg(parentsDepth)

		highBulkDepth := 0
		lowBulkDepth := 0

		highParentDepth := 0
		lowParentDepth := 0

		highBulkSize := 0
		lowBulkSize := 0

		winSize := 0

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

		// ===================================== Get Max Depth caps (optional) ====================================== //
		highBulkMaxDepth := 0
		lowBulkMaxDepth := 0
		highParentMaxDepth := 0
		lowParentMaxDepth := 0

		parsePair := func(name, raw string) (int, int, bool) {
			if strings.TrimSpace(raw) == "" {
				return 0, 0, true // empty = disabled
			}
			parts := strings.Split(raw, ",")
			if len(parts) != 2 {
				color.Red("%s is supposed to be in the form a,b (where a and b are integers)", name)
				return 0, 0, false
			}
			a, errA := strconv.Atoi(strings.TrimSpace(parts[0]))
			b, errB := strconv.Atoi(strings.TrimSpace(parts[1]))
			if errA != nil || errB != nil {
				color.Red("%s is supposed to be in the form a,b (where a and b are integers)", name)
				return 0, 0, false
			}
			return a, b, true
		}

		var okPair bool
		if highBulkMaxDepth, lowBulkMaxDepth, okPair = parsePair("bulks-max-depth", bulksMaxDepth); !okPair {
			return
		}
		if highParentMaxDepth, lowParentMaxDepth, okPair = parsePair("parents-max-depth", parentsMaxDepth); !okPair {
			return
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
			LightFilter:              lightFiltering,
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
			VCF: vcfFile,

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

			Population:             population,
			WindowSize:             winSize,
			StepSize:               step,
			BrmAlpha:               brmAlpha,
			MergeDistance:          mergeDistance,
			MinQTLMarkers:          minQTLMarkers,
			GeneSpaceSpan:          geneSpaceSpan,
			MinGQ:                  minGQ,
			SplitMultiallelic:      !noSplitMultiallelic,
			MinPooledMAF:           minPooledMAF,
			RequireNonParentAllele: requireNonParent,
			RefParent:              refParent,
			OneBulkP0Mode:          oneBulkP0Mode,
			CMPerMb:                cmPerMb,
			NullReps:               nullReps,
			OutputDir:              resultsDir,

			HighParentName:     highParentName,
			HighParentDepth:    highParentDepth,
			HighParentMaxDepth: highParentMaxDepth,

			LowParentName:     lowParentName,
			LowParentDepth:    lowParentDepth,
			LowParentMaxDepth: lowParentMaxDepth,

			HighBulkName:     highBulkName,
			HighBulkDepth:    highBulkDepth,
			HighBulkMaxDepth: highBulkMaxDepth,
			HighBulkSize:     highBulkSize,

			LowBulkName:     lowBulkName,
			LowBulkDepth:    lowBulkDepth,
			LowBulkMaxDepth: lowBulkMaxDepth,
			LowBulkSize:     lowBulkSize,

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
			color.Red("Error: %v", err)
			os.Exit(1)
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
	rootCmd.Flags().StringP("config", "c", "", "Config file path")
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
	rootCmd.Flags().StringP("parents-depth", "p", "5,5", "Parents min depth (hp,lp)")
	rootCmd.Flags().StringP("bulks-depth", "b", "40,40", "Bulks min depth (hb,lb)")
	rootCmd.Flags().StringP("bulk-sizes", "S", "20,20", "Bulk sizes (hb,lb)")
	rootCmd.Flags().String("bulks-max-depth", "", "Bulks max-depth cap hb,lb to drop coverage outliers (empty = disabled)")
	rootCmd.Flags().String("parents-max-depth", "", "Parents max-depth cap hp,lp to drop coverage outliers (empty = disabled)")
	rootCmd.Flags().Int("min-GQ", 0, "Minimum genotype quality (GQ) for PARENT samples; 0 = disabled, applied only when GQ present (not applied to bulks, to avoid AF-correlated bias)")
	rootCmd.Flags().Bool("no-split-multiallelic", false, "Do not decompose multi-allelic records into biallelic records before filtering")
	rootCmd.Flags().Bool("require-nonparent-allele", false, "One parent + one bulk modes: drop sites where the bulk has <2 reads of the allele the parent does not carry (removes non-segregating sites, but also QTL markers where the bulk is fixed for the parent's allele)")
	rootCmd.Flags().Float64("cm-per-mb", 0, "Recombination rate for the HMM LOD, in cM per Mb. 0 = assume 1 Morgan (one crossover per meiosis) over each chromosome")
	rootCmd.Flags().Int("null-reps", 200, "Number of simulated no-QTL genomes for the genome-wide HMM LOD and CompositeZ thresholds")
	rootCmd.Flags().String("one-bulk-p0", "empirical", "One-bulk modes: null allele frequency of the high-parent allele. empirical = genome-wide median per polarity (absorbs reference mapping bias); theoretical = the population's expected value")
	rootCmd.Flags().String("ref-parent", "none", "Bulks-only (2b) mode: which parent the reference genome is (high, low or none). high/low fixes allele polarity so ΔSI is signed; none (unknown, or the reference is not a parent) makes 2b use |ΔSI|. Use high/low only if the reference is that parent or the same inbred line. Ignored when a parent sample is given")
	rootCmd.Flags().Float64("min-pooled-maf", -1, "Segregation floor for 2b/hp2b/lp2b: minimum pooled minor allele frequency of the two bulks (-1 = auto: 0.4 x expected minor AF, 0.2 for F2; 0 = off)")

	// ------------------------------------------------- Smoothing -------------------------------------------------- //
	rootCmd.Flags().Int64P("window-size", "w", 2000000, "Window Size")
	rootCmd.Flags().Int64P("step-size", "s", 100000, "Smoothing grid spacing (bp) on dense chromosomes: kernels are evaluated on this grid and interpolated to markers; capped at window/40 (σ/20)")
	rootCmd.Flags().Int64("min-qtl-width", 0, "Ignored: width follows from smoothing; use --min-qtl-markers")
	_ = rootCmd.Flags().MarkDeprecated("min-qtl-width", "region width follows from smoothing; use --min-qtl-markers")
	rootCmd.Flags().Int64("merge-distance", 0, "Maximum gap in bp to merge same-direction QTL support (0 uses smoothing sigma)")
	rootCmd.Flags().Int("min-qtl-markers", 5, "Minimum significant markers supporting a QTL")

	// ------------------------------------------------- Threshold -------------------------------------------------- //
	rootCmd.Flags().Float64("brm-alpha", 0.05, "Significance level")
	rootCmd.Flags().StringP("population", "m", "F2", "Population type (F2, F3, RIL, BC1H, BC1L, BC2H, BC2L, ...)")
	rootCmd.Flags().Int("rep", 1000, "Ignored: thresholds come from the genome-wide null (--null-reps)")
	_ = rootCmd.Flags().MarkDeprecated("rep", "thresholds now come from the genome-wide null simulation; use --null-reps")

	// ----------------------------------------------- Filtering ---------------------------------------------------- //
	rootCmd.Flags().Bool("light-filtering", false, "Enable light filtering")
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
	rootCmd.Flags().Bool("gene-space-span", false, "Gene space analysis on the final QTLs' span (union of HMM LOD, BRM and CompositeZ regions) instead of their consensus core")
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
