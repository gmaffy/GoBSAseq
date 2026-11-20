package stats

import (
	"io"
	"os"
	"testing"

	"github.com/fatih/color"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/utils"
)

// benchSites tiles goldenSites to n sites with increasing positions.
func benchSites(nSamples, n int) []filter.Site {
	base := goldenSites(nSamples, 1)
	out := make([]filter.Site, n)
	for i := range out {
		out[i] = base[i%len(base)]
		out[i].Pos = int64(i + 1)
		out[i].Chrom = "chr1"
	}
	return out
}

// BenchmarkRawStats: go test ./stats -run '^$' -bench RawStats -benchmem
func BenchmarkRawStats(b *testing.B) {
	for _, m := range []struct {
		mode string
		n    int
	}{{"2p2b", 4}, {"hphb", 2}} {
		b.Run(m.mode, func(b *testing.B) {
			sites := benchSites(m.n, 2_000_000)
			idxs := make([]int, m.n)
			for i := range idxs {
				idxs[i] = i
			}
			cfg := utils.AnalysisConfig{OutputDir: b.TempDir(), Population: "F2", RefParent: "high"}
			color.Output, os.Stderr = io.Discard, devNull(b) // silence logs and progress bar
			b.ResetTimer()
			for range b.N {
				if _, err := RawStats(cfg, m.mode, idxs, sites); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

func devNull(b *testing.B) *os.File {
	f, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
	if err != nil {
		b.Fatal(err)
	}
	b.Cleanup(func() { f.Close() })
	return f
}
