package utils

import "github.com/brentp/vcfgo"

// ---------------------------------------------------------------------------
// 1. Two bulks + two parents  (the standard BSA-seq design)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// 2. Two parents + one bulk
// ---------------------------------------------------------------------------

// passesTwoParentOneBulkFilter validates sites for a design with both parental
// lines and a single bulk (placeholder — logic to be implemented).
func passesTwoParentOneBulkFilter(v *vcfgo.Variant, cfg AnalysisConfig) bool {
	// TODO: implement
	// Expected requirements (to confirm with team):
	//   - Exactly one non-symbolic ALT allele.
	//   - Both parents and the single bulk present, genotyped, REF/ALT only.
	//   - Both parents homozygous and fixed for different alleles.
	//   - All three samples meet minimum depth thresholds.
	_ = v
	_ = cfg
	return false
}

// ---------------------------------------------------------------------------
// 3. Two bulks only (no parental lines)
// ---------------------------------------------------------------------------

// passesBulksOnlyFilter validates sites for a design with two bulks and no
// parental lines (placeholder — logic to be implemented).
func passesBulksOnlyFilter(v *vcfgo.Variant, cfg AnalysisConfig) bool {
	// TODO: implement
	// Expected requirements (to confirm with team):
	//   - Exactly one non-symbolic ALT allele.
	//   - Both bulk samples present, genotyped, and carrying only REF or ALT.
	//   - Both bulks meet minimum depth thresholds.
	//   - (No parental homozygosity constraint — absent from this design.)
	_ = v
	_ = cfg
	return false
}

// ---------------------------------------------------------------------------
// 4. One parent + two bulks
// ---------------------------------------------------------------------------

// passesOneParentBulksOnlyFilter validates sites for a design with one
// parental line and two bulks (placeholder — logic to be implemented).
func passesOneParentBulksOnlyFilter(v *vcfgo.Variant, cfg AnalysisConfig) bool {
	// TODO: implement
	// Expected requirements (to confirm with team):
	//   - Exactly one non-symbolic ALT allele.
	//   - One parent and both bulks present, genotyped, REF/ALT only.
	//   - The single parent homozygous.
	//   - All three samples meet minimum depth thresholds.
	_ = v
	_ = cfg
	return false
}

// ---------------------------------------------------------------------------
// Shared sample validation
// ---------------------------------------------------------------------------
