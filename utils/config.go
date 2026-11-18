package utils

import "github.com/brentp/vcfgo"

type HardFilterConfig struct {
	LightFilter              bool
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

	HighParFwdReads  string
	HighParRevReads  string
	LowParFwdReads   string
	LowParRevReads   string
	HighBulkFwdReads string
	HighBulkRevReads string
	LowBulkFwdReads  string
	LowBulkRevReads  string

	HighParBam  string
	LowParBam   string
	HighBulkBam string
	LowBulkBam  string

	// Parameters for the analysis
	Population string
	WindowSize int
	StepSize   int
	Rep        int
	BrmAlpha   float64
	// MinGQ is the minimum genotype quality (FORMAT/GQ) required of the parent
	// samples at a site. Zero disables the check and is applied only when GQ is
	// present (e.g. DeepVariant VCFs). It is deliberately NOT applied to bulks:
	// in a pool, GQ reflects confidence in a diploid genotype call that is a
	// fiction for a bulk and is lowest at intermediate allele fractions, so a
	// bulk GQ floor would be an allele-frequency-correlated filter that biases
	// the SNP-index. Parents are real genotypes, so GQ is meaningful there.
	MinGQ int
	// SplitMultiallelic decomposes multi-allelic records into biallelic ones
	// before filtering (bcftools norm -m- equivalent). Enabled by default.
	SplitMultiallelic bool
	// Region calling controls. A zero value selects biologically conservative
	// defaults derived from the Gaussian smoothing bandwidth.
	MinQTLWidth   int64
	MergeDistance int64
	MinQTLMarkers int
	OutputDir     string

	HighParentIdx      int
	HighParentName     string
	HighParentDepth    int
	HighParentMaxDepth int

	LowParentIdx      int
	LowParentName     string
	LowParentDepth    int
	LowParentMaxDepth int

	OneParentIdx   int
	OneParentName  string
	OneParentDepth int

	HighBulkIdx      int
	HighBulkName     string
	HighBulkDepth    int
	HighBulkMaxDepth int
	HighBulkSize     int

	LowBulkIdx      int
	LowBulkName     string
	LowBulkDepth    int
	LowBulkMaxDepth int
	LowBulkSize     int

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

	Merger             string
	Caller             string
	NoMerging          bool
	DeepVariantVersion string
	ModelType          string
	Verbose            bool
	Threads            int
}

// EffectiveMergeDistance returns the maximum bp gap bridged when merging adjacent
// QTL intervals. Zero MergeDistance defaults to the Gaussian smoothing sigma.
func (cfg AnalysisConfig) EffectiveMergeDistance() int64 {
	if cfg.MergeDistance > 0 {
		return cfg.MergeDistance
	}
	if cfg.WindowSize > 0 {
		return int64(cfg.WindowSize / 2)
	}
	return 0
}
