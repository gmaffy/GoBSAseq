package stats

import (
	"io"
	"math"
	"math/rand"
	"os"
	"sort"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
	"gonum.org/v1/gonum/stat/distuv"
)

// plantSim is a plant-level F2 generator, independent of the HMM's pool model.
type plantSim struct {
	nChrom  int
	lenBp   float64
	markers []int64
}

// gamete returns, for each marker, whether the gamete carries the high allele.
func (p plantSim) gamete(rng *rand.Rand) []bool {
	g := make([]bool, len(p.markers))
	state := rng.Intn(2) == 1
	nx := distuv.Poisson{Lambda: 1, Src: rng}.Rand() // 1 Morgan
	xs := make([]float64, int(nx))
	for i := range xs {
		xs[i] = rng.Float64() * p.lenBp
	}
	sort.Float64s(xs)
	k := 0
	for i, pos := range p.markers {
		for k < len(xs) && xs[k] < float64(pos) {
			state = !state
			k++
		}
		g[i] = state
	}
	return g
}

func (p plantSim) dataset(rng *rand.Rand, nPlants int, rho float64, qtl bool) []BSAstats {
	var raw []BSAstats
	calc := rawCalc{roles: siteRoles{hp: -1, lp: -1, hb: 0, lb: 1}, signedDSI: true}
	for c := 0; c < p.nChrom; c++ {
		copies := [2][]int{make([]int, len(p.markers)), make([]int, len(p.markers))}
		for b := 0; b < 2; b++ {
			for range nPlants {
				var g1, g2 []bool
				for {
					g1, g2 = p.gamete(rng), p.gamete(rng)
					if !qtl || c != 0 {
						break
					}
					// QTL at the chromosome middle: high bulk needs ≥1 high allele, low bulk none (recessive-ish selection)
					mid := len(p.markers) / 2
					n := 0
					if g1[mid] {
						n++
					}
					if g2[mid] {
						n++
					}
					if (b == 0 && n == 2) || (b == 1 && n == 0) {
						break
					}
				}
				for i := range p.markers {
					if g1[i] {
						copies[b][i]++
					}
					if g2[i] {
						copies[b][i]++
					}
				}
			}
		}
		for i, pos := range p.markers {
			site := filter.Site{Chrom: string(rune('A' + c)), Pos: pos, Ref: "A", Alt: "T", NSamples: 2}
			for b := 0; b < 2; b++ {
				d := int(distuv.Poisson{Lambda: distuv.Gamma{Alpha: 8, Beta: 8.0 / 50, Src: rng}.Rand(), Src: rng}.Rand())
				f := float64(copies[b][i]) / float64(2*nPlants)
				pr := f*0.99 + (1-f)*0.01
				if rho > 0 { // beta-binomial overdispersion
					a := pr * (1 - rho) / rho
					pr = distuv.Beta{Alpha: a, Beta: a * (1 - pr) / pr, Src: rng}.Rand()
				}
				h := binomSample(rng, d, pr)
				site.Samples[b] = filter.SampleAD{Ref: int32(h), Alt: int32(d - h), HasAD: true}
			}
			raw = append(raw, calc.site(&site))
		}
	}
	return raw
}

// TestHMMLODCalibration checks the genome-wide null against data from an
// independent plant-level generator (F2 plants with Poisson crossovers,
// gamma-Poisson depths, optional beta-binomial overdispersion): about 5% of
// no-QTL genomes should exceed the p95 HMM LOD threshold, 1% the p99, and 1%
// of sites the p99 CompositeZ threshold. Slow (~90 s): set
// GOBSASEQ_CALIBRATION=1 to run.
func TestHMMLODCalibration(t *testing.T) {
	if os.Getenv("GOBSASEQ_CALIBRATION") == "" {
		t.Skip("set GOBSASEQ_CALIBRATION=1 to run the calibration check")
	}
	color.Output = io.Discard
	p := plantSim{nChrom: 3, lenBp: 30e6}
	for x := 1500.0; x < p.lenBp; x += 3000 {
		p.markers = append(p.markers, int64(x))
	}
	cfg := utils.AnalysisConfig{WindowSize: 2_000_000, Population: "F2", HighBulkSize: 20, LowBulkSize: 20, RefParent: "high", NullReps: 200}
	sc, _ := newScorer(cfg, "2p2b")
	for _, rho := range []float64{0, 0.01} {
		rng := rand.New(rand.NewSource(11))
		tmpl := p.dataset(rng, 20, rho, false)
		zt, ht, err := genomeNullThresholds(cfg, "2p2b", tmpl)
		if err != nil {
			t.Fatal(err)
		}
		var over95, over99, czOver99, czN int
		const nNull = 100
		for range nNull {
			raw := p.dataset(rng, 20, rho, false)
			sm, _ := sc.score(raw, false, 1)
			mx := 0.0
			for i := range raw {
				mx = math.Max(mx, raw[i].HMMLOD)
			}
			if mx > ht.LODP95 {
				over95++
			}
			if mx > ht.LODP99 {
				over99++
			}
			for i := range sm {
				czN++
				if sm[i].CompositeZ > zt.CompositeZP99 {
					czOver99++
				}
			}
		}
		// power: one QTL dataset
		qraw := p.dataset(rng, 20, rho, true)
		sc.score(qraw, false, 1)
		qmax := 0.0
		for i := range qraw {
			if qraw[i].CHROM == "A" {
				qmax = math.Max(qmax, qraw[i].HMMLOD)
			}
		}
		if over99 > 5 || over95 > 12 {
			t.Errorf("rho=%.3f: null genomes over p95 %d/%d, over p99 %d/%d: thresholds too low", rho, over95, nNull, over99, nNull)
		}
		t.Logf("rho=%.3f: HMM LOD thresholds p95 %.2f p99 %.2f | null genomes over p95: %d/%d, over p99: %d/%d | CompositeZ sites over p99: %.2f%% | QTL chromosome max LOD %.1f",
			rho, ht.LODP95, ht.LODP99, over95, nNull, over99, nNull, 100*float64(czOver99)/float64(czN), qmax)
	}
}
