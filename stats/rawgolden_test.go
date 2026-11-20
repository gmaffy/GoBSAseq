package stats

import (
	"bytes"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
)

var updateGolden = flag.Bool("update-golden", false, "rewrite stats/testdata/raw_golden from the current code")

var allModes = []struct {
	mode string
	n    int
}{
	{"2p2b", 4}, {"2phb", 3}, {"2plb", 3}, {"hp2b", 3}, {"lp2b", 3},
	{"hphb", 2}, {"hplb", 2}, {"lphb", 2}, {"lplb", 2}, {"2b", 2},
}

// goldenSites builds deterministic sites covering polarity flips, zero-depth
// and one-allele bulks, other-ALT reads and odd genotypes. 600 sites so the
// one-bulk empirical p0 path (≥200 per polarity class) is exercised.
func goldenSites(nSamples int, seed int64) []filter.Site {
	rng := rand.New(rand.NewSource(seed))
	sites := make([]filter.Site, 600)
	for i := range sites {
		s := filter.Site{Chrom: fmt.Sprintf("chr%d", 1+i/300), Pos: int64(1000 + 37*i), Ref: "A", Alt: "T",
			HighIsAlt: rng.Intn(2) == 1, NSamples: uint8(nSamples)}
		if i%50 == 7 {
			s.Ref, s.Alt = "AT", "A"
		}
		for k := 0; k < nSamples; k++ {
			d := int32(5 + rng.Intn(120))
			ref := int32(rng.Intn(int(d) + 1))
			sa := filter.SampleAD{Ref: ref, Alt: d - ref, DP: d, GQ: -1, HasAD: true}
			switch rng.Intn(12) {
			case 0:
				sa.Ref, sa.Alt = 0, 0 // zero depth
			case 1:
				sa.Alt = 0 // one allele only
			case 2:
				sa.Ref = 0
			case 3:
				sa.Other = int32(1 + rng.Intn(4))
			}
			switch rng.Intn(6) {
			case 0:
				sa.GTRef = 2
			case 1:
				sa.GTAlt = 2
			case 2:
				sa.GTMissing = 2
			case 3:
				sa.GTRef, sa.GTOther = 1, 1
			default:
				sa.GTRef, sa.GTAlt = 1, 1
			}
			s.Samples[k] = sa
		}
		sites[i] = s
	}
	return sites
}

// numericDump writes every numeric BSAstats field at full precision, so the
// golden test also pins values that the TSV rounds.
func numericDump(raw []BSAstats) []byte {
	var b bytes.Buffer
	g := func(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }
	for _, s := range raw {
		fmt.Fprintln(&b, s.CHROM, s.POS, s.REF, s.ALT, s.HighBulkL, s.HighBulkH, s.LowBulkL, s.LowBulkH,
			g(s.HighSI), g(s.LowSI), g(s.DeltaSI), g(s.Gstat), g(s.ED),
			g(s.OneBulkP0), g(s.OneBulkAFDev), g(s.OneBulkGstat), s.HighIsAlt,
			s.HighBulkOther, s.LowBulkOther, s.HighDepth, s.LowDepth, s.Depth)
	}
	return b.Bytes()
}

// TestRawStatsGolden pins RawStats output (TSV bytes and full-precision
// values) for all 10 modes. Refactors must leave it unchanged; regenerate with
// go test ./stats -run Golden -update-golden only for intended result changes.
func TestRawStatsGolden(t *testing.T) {
	for mi, m := range allModes {
		for _, rp := range []string{"none", "low"} {
			if rp == "low" && m.mode != "2b" {
				continue
			}
			name := m.mode
			if m.mode == "2b" {
				name += "_" + rp
			}
			t.Run(name, func(t *testing.T) {
				idxs := make([]int, m.n)
				for i := range idxs {
					idxs[i] = i
				}
				cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), Population: "F2", RefParent: rp}
				raw, err := RawStats(cfg, m.mode, idxs, goldenSites(m.n, int64(mi+1)))
				if err != nil {
					t.Fatal(err)
				}
				tsv, err := os.ReadFile(filepath.Join(cfg.OutputDir, "stats", "GoBSAseq."+m.mode+".raw.tsv"))
				if err != nil {
					t.Fatal(err)
				}
				got := map[string][]byte{name + ".tsv": tsv, name + ".values": numericDump(raw)}
				for file, data := range got {
					path := filepath.Join("testdata", "raw_golden", file)
					if *updateGolden {
						if err := os.WriteFile(path, data, 0o644); err != nil {
							t.Fatal(err)
						}
						continue
					}
					want, err := os.ReadFile(path)
					if err != nil {
						t.Fatalf("%v (run with -update-golden to create)", err)
					}
					if !bytes.Equal(data, want) {
						t.Errorf("%s differs from golden", file)
					}
				}
			})
		}
	}
}
