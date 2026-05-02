package utils

import (
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
	OutputFile    string
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
