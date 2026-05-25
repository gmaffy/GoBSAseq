package utils

import (
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/tabix"
	"github.com/brentp/vcfgo"
)

func getFloat(v *vcfgo.Variant, key string) (float64, bool) {
	raw, err := v.Info().Get(key)
	if err != nil || raw == nil {
		return 0, false
	}
	switch val := raw.(type) {
	case float32:
		return float64(val), true
	case float64:
		return val, true
	case int:
		return float64(val), true
	default:
		return 0, false
	}
}

// ---------------------------------------------------------------------------
// Variant-class helpers
// ---------------------------------------------------------------------------

// isHomozygous returns true when all alleles in gt are identical and
// non-missing.
func isHomozygous(gt []int) bool {
	if len(gt) == 0 {
		return false
	}
	first := gt[0]
	if first < 0 {
		return false
	}
	for _, a := range gt[1:] {
		if a < 0 || a != first {
			return false
		}
	}
	return true
}

// ---------------------------------------------------------------------------
// Tabix helpers
// ---------------------------------------------------------------------------

type vcfRecord struct {
	chrom string
	start int // 0-based
	end   int // 0-based half-open
}

func (r vcfRecord) RefName() string { return r.chrom }
func (r vcfRecord) Start() int      { return r.start }
func (r vcfRecord) End() int        { return r.end }

// newTabixIndex returns a tabix index configured for VCF.
func newTabixIndex() *tabix.Index {
	idx := tabix.New()
	idx.Format = 2 // VCF
	idx.NameColumn = 1
	idx.BeginColumn = 2
	idx.EndColumn = 0
	idx.MetaChar = '#'
	return idx
}

// addTabixRecord records the bgzf chunk for a variant in idx.
func addTabixRecord(idx *tabix.Index, v *vcfgo.Variant, chunk bgzf.Chunk) error {
	rec := vcfRecord{
		chrom: v.Chromosome,
		start: int(v.Pos) - 1,
		end:   int(v.Pos) - 1 + len(v.Ref()),
	}
	return idx.Add(rec, chunk, true, true)
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

func TricubeWeight(d, D float64) float64 {
	if D <= 0 || d >= D {
		return 0
	}
	x := 1 - math.Pow(d/D, 3)
	return x * x * x
}
