package filter

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/bgzf"
)

func TestReg2Bin(t *testing.T) {
	cases := []struct {
		beg, end int
		want     uint32
	}{
		{0, 1, 4681},          // smallest 16 kb bin
		{16384, 16385, 4682},  // next 16 kb bin
		{0, 16385, 585},       // spans two 16 kb bins -> 128 kb bin
		{0, 1 << 26, 1},       // first 64 Mb bin
		{0, (1 << 26) + 1, 0}, // whole-reference bin
		{(1 << 26), (1 << 26) + 1, 4681 + (1<<26)>>14},
	}
	for _, c := range cases {
		if got := reg2bin(c.beg, c.end); got != c.want {
			t.Errorf("reg2bin(%d,%d) = %d, want %d", c.beg, c.end, got, c.want)
		}
	}
}

// TestBuildTabixIndexRefs guards against the biogo/hts v1.4.5 tabix bug where
// every record became its own reference: the index must list each contig once.
func TestBuildTabixIndexRefs(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "t.vcf.gz")

	f, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	w := bgzf.NewWriter(f, 1)
	fmt.Fprintln(w, "##fileformat=VCFv4.2")
	fmt.Fprintln(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	contigs := []string{"chr1", "chr2", "chr3"}
	for _, c := range contigs {
		// Enough records to span several BGZF blocks.
		for i := 0; i < 3000; i++ {
			fmt.Fprintf(w, "%s\t%d\t.\tA\tG\t50\tPASS\t%s\n", c, 1+i*1000, strings.Repeat("X", 20))
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	f.Close()

	if err := buildTabixIndex(path); err != nil {
		t.Fatalf("buildTabixIndex: %v", err)
	}

	tf, err := os.Open(path + ".tbi")
	if err != nil {
		t.Fatal(err)
	}
	defer tf.Close()
	r, err := bgzf.NewReader(tf, 1)
	if err != nil {
		t.Fatal(err)
	}
	raw, err := io.ReadAll(r)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.HasPrefix(raw, []byte("TBI\x01")) {
		t.Fatalf("bad magic %q", raw[:4])
	}
	nRef := int32(binary.LittleEndian.Uint32(raw[4:8]))
	if nRef != int32(len(contigs)) {
		t.Fatalf("n_ref = %d, want %d", nRef, len(contigs))
	}
	lNm := int32(binary.LittleEndian.Uint32(raw[32:36]))
	names := strings.Split(strings.TrimSuffix(string(raw[36:36+lNm]), "\x00"), "\x00")
	if strings.Join(names, ",") != strings.Join(contigs, ",") {
		t.Fatalf("names = %v, want %v", names, contigs)
	}
}
