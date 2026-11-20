package utils

import (
	"fmt"
	"strings"
)

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
	// Rdr *vcfgo.Reader // LEGACY: the vcfgo reader used by the old filter.
	// VCFInfoIDs lists the INFO keys declared in the input VCF header.
	VCFInfoIDs map[string]bool

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
	// MinPooledMAF is the segregation floor for two-bulk modes without two
	// informative parents (2b, hp2b, lp2b): the pooled minor allele frequency of
	// both bulks must be at least this. < 0 = automatic (0.4 × the population's
	// expected minor allele frequency, 0.2 for F2); 0 = disabled.
	MinPooledMAF float64
	// RequireNonParentAllele (one-parent one-bulk modes) drops sites where the
	// bulk shows fewer than 2 reads of the allele the parent does not carry.
	// Off by default: such sites may be non-segregating, but they may also be
	// QTL markers where the bulk is fixed for the parent's allele.
	RequireNonParentAllele bool
	// RefParent says which parent the reference genome is, for bulks-only (2b)
	// mode: "high" (REF is the high-parent allele), "low" (ALT is) or "none"
	// (unknown, or the reference is not a parent). With "none" the sign of ΔSI
	// is arbitrary at each SNP, so 2b uses |ΔSI| instead (see SignedDeltaSI).
	// Ignored in modes with a parent sample.
	RefParent string
	// OneBulkP0Mode sets the one-bulk null allele frequency: "empirical" (default;
	// also when empty) estimates it per polarity from the data to absorb
	// reference mapping bias, "theoretical" uses the population's expected value.
	OneBulkP0Mode string
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
	// OneBulkSize is not set from the CLI; use BulkSizes(bsaType) instead.
	OneBulkSize int

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

// SignedDeltaSI reports whether ΔSI has a meaningful sign for bsaType: true
// whenever polarity is known (a parent sample, or --ref-parent high|low in 2b).
// When false (2b with --ref-parent none) the pipeline stores, smooths and
// thresholds |ΔSI|, has no lower-tail ΔSI track and skips two-bulk BRM.
func (cfg AnalysisConfig) SignedDeltaSI(bsaType string) bool {
	return bsaType != "2b" || cfg.RefParent == "high" || cfg.RefParent == "low"
}

// BulkSizes returns the number of individuals in each bulk used by bsaType.
// high/low are the two-bulk sizes; one is the size of the single bulk in
// one-bulk modes (LowBulkSize for *lb modes, HighBulkSize otherwise). Bulk
// size sets the pool-sampling variance of every null model, so a missing size
// is an error rather than a silent default.
func (cfg AnalysisConfig) BulkSizes(bsaType string) (high, low, one int, err error) {
	hasHigh := strings.Contains(bsaType, "hb") || strings.Contains(bsaType, "2b")
	hasLow := strings.Contains(bsaType, "lb") || strings.Contains(bsaType, "2b")
	if hasHigh && cfg.HighBulkSize <= 0 {
		return 0, 0, 0, fmt.Errorf("high bulk size must be > 0 for %s (set -S/--bulk-sizes), got %d", bsaType, cfg.HighBulkSize)
	}
	if hasLow && cfg.LowBulkSize <= 0 {
		return 0, 0, 0, fmt.Errorf("low bulk size must be > 0 for %s (set -S/--bulk-sizes), got %d", bsaType, cfg.LowBulkSize)
	}
	high, low = cfg.HighBulkSize, cfg.LowBulkSize
	switch {
	case hasHigh && !hasLow:
		one = high
	case hasLow && !hasHigh:
		one = low
	}
	return high, low, one, nil
}
