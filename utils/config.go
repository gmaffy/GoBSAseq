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

	SnpEffDB string
	Ref      string
	Protein  string
	Gff      string
	Cds      string
	GeneDesc string
	Prg      string
}

type HardFilterConfig struct {
	// SNP thresholds (GATK Best Practices)
	SNP_QD_Min             float64 // default 2.0   – QualByDepth
	SNP_QUAL_Min           float64 // default 30.0  – variant quality score
	SNP_SOR_Max            float64 // default 3.0   – StrandOddsRatio
	SNP_FS_Max             float64 // default 60.0  – FisherStrand
	SNP_MQ_Min             float64 // default 40.0  – RMSMappingQuality
	SNP_MQRankSum_Min      float64 // default -12.5 – MappingQualityRankSumTest
	SNP_ReadPosRankSum_Min float64 // default -8.0  – ReadPosRankSumTest

	INDEL_QD_Min             float64 // default 2.0   – QualByDepth
	INDEL_QUAL_Min           float64 // default 30.0  – variant quality score
	INDEL_FS_Max             float64 // default 200.0 – FisherStrand
	INDEL_ReadPosRankSum_Min float64 // default -20.0 – ReadPosRankSumTest

	// SaveFilteredVCF writes the hard-filtered records (PASS only) to a
	// bgzf-compressed VCF at FilteredVCFPath when true.

}
