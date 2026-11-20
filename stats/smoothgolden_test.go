package stats

import (
	"bytes"
	"fmt"
	"io"
	"math/rand"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/utils"
)

// smoothInput builds deterministic raw stats over two 20 Mb chromosomes with
// irregular spacing, a sparse stretch, duplicate positions (split
// multi-allelic) and a 3-marker contig that the sparse-kernel rule drops.
func smoothInput(seed int64, nPerChrom int) []BSAstats {
	rng := rand.New(rand.NewSource(seed))
	var raw []BSAstats
	add := func(chrom string, pos int64) {
		hd, ld := 5+rng.Intn(150), 5+rng.Intn(150)
		hsi, lsi := rng.Float64(), rng.Float64()
		// Draw in the original order (the retired LOD/BBLogBF draws are kept
		// and discarded) so the input stays comparable with earlier goldens.
		g, ed := rng.ExpFloat64()*3, rng.Float64()*0.2
		_, _ = rng.ExpFloat64(), rng.NormFloat64()
		g1 := rng.ExpFloat64() * 2
		_, _ = rng.ExpFloat64(), rng.NormFloat64()
		raw = append(raw, BSAstats{
			CHROM: chrom, POS: pos, REF: "A", ALT: "T",
			HighSI: hsi, LowSI: lsi, DeltaSI: hsi - lsi,
			Gstat: g, ED: ed,
			OneBulkP0: 0.5, OneBulkAFDev: hsi - 0.5, OneBulkGstat: g1,
			HighBulkH: int(hsi * float64(hd)), HighBulkL: hd - int(hsi*float64(hd)),
			LowBulkH: int(lsi * float64(ld)), LowBulkL: ld - int(lsi*float64(ld)),
			HighDepth: hd, LowDepth: ld, Depth: effectiveDepth(hd, ld),
		})
	}
	for _, chrom := range []string{"chr2", "chr10"} { // lexical order differs from VCF order
		var pos []int64
		for len(pos) < nPerChrom {
			p := int64(1 + rng.Intn(20_000_000))
			if p > 8_000_000 && p < 14_000_000 && rng.Intn(50) != 0 {
				continue // sparse stretch
			}
			pos = append(pos, p)
		}
		sort.Slice(pos, func(a, b int) bool { return pos[a] < pos[b] })
		for i, p := range pos {
			add(chrom, p)
			if i%97 == 0 {
				add(chrom, p) // duplicate position
			}
		}
	}
	for _, p := range []int64{100, 5000, 9000} {
		add("contig_tiny", p)
	}
	return raw
}

func dumpSmoothed(sm []SmoothedStats) []byte {
	var b bytes.Buffer
	g := func(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }
	for _, s := range sm {
		fmt.Fprintln(&b, s.CHROM, s.POS, s.Depth, s.HighDepth, s.LowDepth, g(s.NEff),
			g(s.SmHighSI), g(s.SmLowSI), g(s.SmDeltaSI), g(s.SmGstat), g(s.SmED),
			g(s.OneBulkP0), g(s.SmAFDev), g(s.SmOneBulkG), g(s.HMMLOD),
			g(s.ZHighSI), g(s.ZLowSI), g(s.ZDeltaSI), g(s.ZGstat), g(s.ZED),
			g(s.ZAFDev), g(s.ZOneBulkG), g(s.ZHMMLOD), g(s.CompositeZ), g(s.MaxAbsZ))
	}
	return b.Bytes()
}

// TestSmoothGolden pins SmoothAndNormalise output (TSV bytes and
// full-precision values) for two-bulk signed, two-bulk |ΔSI|, high-bulk-only
// and low-bulk-only layouts. Regenerate with -update-golden only for intended
// result changes.
func TestSmoothGolden(t *testing.T) {
	color.Output = io.Discard
	for _, tc := range []struct{ name, mode, refParent string }{
		{"2p2b", "2p2b", ""}, {"2b_abs", "2b", "none"}, {"hphb", "hphb", ""}, {"2plb", "2plb", ""},
	} {
		t.Run(tc.name, func(t *testing.T) {
			cfg := utils.AnalysisConfig{OutputDir: t.TempDir(), WindowSize: 2_000_000, RefParent: tc.refParent,
				Population: "F2", HighBulkSize: 20, LowBulkSize: 20}
			sm, err := SmoothAndNormalise(cfg, tc.mode, smoothInput(11, 1200))
			if err != nil {
				t.Fatal(err)
			}
			tsv, err := os.ReadFile(filepath.Join(cfg.OutputDir, "stats", "GoBSAseq."+tc.mode+".smoothed_and_normalised.tsv"))
			if err != nil {
				t.Fatal(err)
			}
			for file, data := range map[string][]byte{tc.name + ".tsv": tsv, tc.name + ".values": dumpSmoothed(sm)} {
				path := filepath.Join("testdata", "smooth_golden", file)
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
