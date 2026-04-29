/*
Copyright © 2026 NAME HERE mafireyi@gmail.com
*/
package cmd

import (
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/run"

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
		if !cmd.Flags().Changed("variant") {
			cmd.Help()
			return
		}

		variant, _ := cmd.Flags().GetString("variant")
		//parents, _ := cmd.Flags().GetString("parents")
		//bulks, _ := cmd.Flags().GetString("bulks")
		parentsDepth, _ := cmd.Flags().GetString("parents-depth")
		bulksDepth, _ := cmd.Flags().GetString("bulks-depth")
		bulkSizes, _ := cmd.Flags().GetString("bulk-sizes")
		windowSize, _ := cmd.Flags().GetInt64("window-size")
		population, _ := cmd.Flags().GetString("population")
		recurrent, _ := cmd.Flags().GetBool("recurrent")
		rep, _ := cmd.Flags().GetInt("rep")
		alpha, _ := cmd.Flags().GetFloat64("alpha")
		minQTL, _ := cmd.Flags().GetInt64("min-qtl-length")
		mergeDist, _ := cmd.Flags().GetInt64("merge-distance")
		outputDir, _ := cmd.Flags().GetString("out")

		//bulksDepthLst := strings.Split(bulksDepth, ",")
		bulksDepthLst := strings.Split(bulksDepth, ",")
		bulkSizesLst := strings.Split(bulkSizes, ",")
		parentsDepthLst := strings.Split(parentsDepth, ",")

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

		// ========================================== Get bulk depths =============================================== //
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

		// =========================================== Get Parent Depths============================================= //
		if len(parentsDepthLst) > 2 {
			color.Red("parentsDepth is supposed to be in the form a,b (where a and b are integers)")
			return
		} else if len(parentsDepthLst) == 1 {
			oneParentDepth, err = strconv.Atoi(parentsDepthLst[0])

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

		// ========================================== Get Bulk Sizes ================================================ //

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

		// ============================================== window size ================================================//
		winSize, err = strconv.Atoi(fmt.Sprintf("%d", windowSize))
		if err != nil {
			color.Red("windowSize is supposed to be an integer")
			return
		}

		err = run.Run(variant, highParentDepth, lowParentDepth, oneParentDepth, highBulkDepth, lowBulkDepth, oneBulkDepth, highBulkSize, lowBulkSize, oneBulkSize, winSize, population, recurrent, rep, alpha, minQTL, mergeDist, outputDir)
		if err != nil {
			return
		}
		fmt.Printf("variant: %s\n", variant)
	},
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	err := rootCmd.Execute()
	if err != nil {
		os.Exit(1)
	}
}

func init() {
	rootCmd.Flags().StringP("variant", "V", "", "Variant File")
	rootCmd.Flags().StringP("parents", "P", "", "parent names")
	rootCmd.Flags().StringP("bulks", "B", "", "bulk names")
	rootCmd.Flags().StringP("parents-depth", "p", "5,5", "Parents Min Depth")
	rootCmd.Flags().StringP("bulks-depth", "b", "40,40", "Low Parent Min Depth")
	rootCmd.Flags().StringP("bulk-sizes", "S", "20,20", "High Bulk Min Depth")
	rootCmd.Flags().Int64P("window-size", "w", 2000000, "Window Size")
	rootCmd.Flags().StringP("population", "m", "F2", "Population type (F2, F3, BC, RIL)")
	rootCmd.Flags().Bool("recurrent", false, "BCAltIsRecurrent: if true, alt allele is recurrent in BC")
	rootCmd.Flags().Int("rep", 1000, "Number of simulations")
	rootCmd.Flags().Float64("alpha", 0.05, "Significance level")
	rootCmd.Flags().Int64("min-qtl-length", 100000, "Minimum QTL length")
	rootCmd.Flags().Int64("merge-distance", 500000, "Merge distance for QTLs")
	rootCmd.Flags().StringP("out", "o", "bsaseq_results", "Output directory/prefix")
}
