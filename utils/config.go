package utils

import "github.com/brentp/vcfgo"

type HardFilterConfig struct {
	SNP_QD_Min               float64
	SNP_QUAL_Min             float64
	SNP_SOR_Max              float64
	SNP_FS_Max               float64
	SNP_MQ_Min               float64
	SNP_MQRankSum_Min        float64
	SNP_ReadPosRankSum_Min   float64
	INDEL_QD_Min             float64
	INDEL_QUAL_Min           float64
	INDEL_FS_Max             float64
	INDEL_ReadPosRankSum_Min float64
	INDEL_SOR_Max            float64
}

type AnalysisConfig struct {
	VCF string
	Rdr *vcfgo.Reader
	// Parameters for the analysis
	Population    string
	WindowSize    int
	StepSize      int
	Rep           int
	BrmAlpha      float64
	MinQTLWidth   int64
	MergeDistance int64
	OutputDir     string

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
