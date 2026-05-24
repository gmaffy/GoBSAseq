package utils

type HardFilterConfig struct {
	SNP_QD_Min             float64
	SNP_QUAL_Min           float64
	SNP_SOR_Max            float64
	SNP_FS_Max             float64
	SNP_MQ_Min             float64
	SNP_MQRankSum_Min      float64
	SNP_ReadPosRankSum_Min float64

	INDEL_QD_Min             float64
	INDEL_QUAL_Min           float64
	INDEL_FS_Max             float64
	INDEL_ReadPosRankSum_Min float64
	INDEL_SOR_Max            float64
}

type AnalysisConfig struct {
	VCF           string
	Population    string
	WindowSize    int
	StepSize      int
	Rep           int
	Alphas        []float64
	MinQTLWidth   int64
	MergeDistance int64
	OutputDir     string

	HighParentName  string
	HighParentDepth int
	LowParentName   string
	LowParentDepth  int
	OneParentName   string
	OneParentDepth  int

	HighBulkName  string
	HighBulkDepth int
	HighBulkSize  int
	LowBulkName   string
	LowBulkDepth  int
	LowBulkSize   int
	OneBulkName   string
	OneBulkDepth  int
	OneBulkSize   int

	SnpEffDB string
	Ref      string
	Protein  string
	Gff      string
	Cds      string
	GeneDesc string
	Prg      string
}
