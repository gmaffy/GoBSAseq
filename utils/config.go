package utils

import (
	"fmt"
	"math/rand/v2"

	"github.com/brentp/vcfgo"
)

type AnalysisConfig struct {
	VCF string
	Rdr *vcfgo.Reader
	// Parameters for the analysis
	Population    string
	WindowSize    int
	StepSize      int
	Rep           int
	Alphas        []float64
	MinQTLWidth   int64
	MergeDistance int64
	OutputDir     string
	//SampleNames      []string // for TSV header
	HighParentIdx   int
	HighParentName  string
	HighParentDepth int

	LowParentIdx   int
	LowParentName  string
	LowParentDepth int

	OneParentIdx   int
	OneParentName  string
	OneParentDepth int

	HighBulkIdx   int
	HighBulkName  string
	HighBulkDepth int
	HighBulkSize  int

	LowBulkIdx   int
	LowBulkName  string
	LowBulkDepth int
	LowBulkSize  int

	OneBulkIdx   int
	OneBulkName  string
	OneBulkDepth int
	OneBulkSize  int
}

func SimulateAF(popStruc string, bulkSize float64, rep int) float64 {

	var prob []float64

	switch popStruc {
	case "F2":
		prob = []float64{0.25, 0.5, 0.25}
	case "RIL":
		prob = []float64{0.5, 0.0, 0.5}
	case "BC":
		prob = []float64{0.5, 0.5, 0.0}
	default:
		fmt.Println("Invalid population structure")
		return 0.0
	}

	var totalFreq float64
	for i := 0; i < rep; i++ {
		var sumFreq float64
		for j := 0; j < int(bulkSize); j++ {
			// Use a simple weighted random choice
			r := rand.Float64()
			var allele float64
			if r < prob[0] {
				allele = 0.0
			} else if r < prob[0]+prob[1] {
				allele = 0.5
			} else {
				allele = 1.0
			}
			sumFreq += allele
		}
		totalFreq += sumFreq / bulkSize
	}
	return totalFreq / float64(rep)

}
