package filter

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"os"

	"github.com/biogo/hts/bgzf"
)

// TBI writer following the SAMtools tabix specification. It is written here
// rather than with biogo/hts/tabix because that package (v1.4.5) never records
// reference names in its name map (every record becomes a new reference) and
// its linear index panics on records that end on a 16 kb tile boundary.

const (
	tbiTileShift = 14      // 16 kb linear-index tiles
	tbiMaxPos    = 1 << 29 // TBI cannot address positions beyond 2^29
)

// reg2bin is the UCSC binning scheme used by BAI/TBI (end is exclusive).
func reg2bin(beg, end int) uint32 {
	end--
	switch {
	case beg>>14 == end>>14:
		return uint32(((1<<15)-1)/7 + (beg >> 14))
	case beg>>17 == end>>17:
		return uint32(((1<<12)-1)/7 + (beg >> 17))
	case beg>>20 == end>>20:
		return uint32(((1<<9)-1)/7 + (beg >> 20))
	case beg>>23 == end>>23:
		return uint32(((1<<6)-1)/7 + (beg >> 23))
	case beg>>26 == end>>26:
		return uint32(((1<<3)-1)/7 + (beg >> 26))
	}
	return 0
}

type tbiChunk struct{ beg, end uint64 }

type tbiRef struct {
	name     string
	binOrder []uint32
	bins     map[uint32][]tbiChunk
	linear   []uint64
}

func voff(file int64, block uint16) uint64 { return uint64(file)<<16 | uint64(block) }

// writeTabix builds the .tbi from exact virtual offsets. Unsorted input or
// positions beyond 2^29 return an error; the caller then skips the index
// rather than failing the run.
func writeTabix(path string, bw *bgzfWriter, entries []indexEntry) error {
	var refs []*tbiRef
	seen := map[string]bool{}
	var cur *tbiRef
	lastBeg := -1
	for _, e := range entries {
		if cur == nil || cur.name != e.chrom {
			if seen[e.chrom] {
				return fmt.Errorf("input is not sorted: %s appears in more than one block", e.chrom)
			}
			seen[e.chrom] = true
			cur = &tbiRef{name: e.chrom, bins: map[uint32][]tbiChunk{}}
			refs = append(refs, cur)
			lastBeg = -1
		}
		if e.beg < lastBeg {
			return fmt.Errorf("input is not sorted at %s:%d", e.chrom, e.beg+1)
		}
		if e.end > tbiMaxPos {
			return fmt.Errorf("%s:%d is beyond the 2^29 bp TBI limit (index with `tabix -C` for CSI)", e.chrom, e.beg+1)
		}
		lastBeg = e.beg
		end := e.end
		if end <= e.beg {
			end = e.beg + 1
		}
		bf, bo := bw.virtualOffset(e.begV)
		ef, eo := bw.virtualOffset(e.endV)
		c := tbiChunk{voff(bf, bo), voff(ef, eo)}

		bin := reg2bin(e.beg, end)
		chunks, ok := cur.bins[bin]
		if !ok {
			cur.binOrder = append(cur.binOrder, bin)
		}
		if n := len(chunks); n > 0 && chunks[n-1].end >= c.beg {
			if c.end > chunks[n-1].end {
				chunks[n-1].end = c.end
			}
		} else {
			chunks = append(chunks, c)
		}
		cur.bins[bin] = chunks

		for t := e.beg >> tbiTileShift; t <= (end-1)>>tbiTileShift; t++ {
			for len(cur.linear) <= t {
				cur.linear = append(cur.linear, 0)
			}
			if cur.linear[t] == 0 {
				cur.linear[t] = c.beg
			}
		}
	}
	// Tiles without records point at the previous tile's offset (htslib does
	// the same), so a query starting in an empty tile still scans forward.
	for _, r := range refs {
		for t := 1; t < len(r.linear); t++ {
			if r.linear[t] == 0 {
				r.linear[t] = r.linear[t-1]
			}
		}
	}

	var b bytes.Buffer
	le := binary.LittleEndian
	w32 := func(v int32) { binary.Write(&b, le, v) }
	b.WriteString("TBI\x01")
	w32(int32(len(refs)))
	w32(2)   // format: VCF
	w32(1)   // col_seq
	w32(2)   // col_beg
	w32(0)   // col_end
	w32('#') // meta char
	w32(0)   // skip
	var names bytes.Buffer
	for _, r := range refs {
		names.WriteString(r.name)
		names.WriteByte(0)
	}
	w32(int32(names.Len()))
	b.Write(names.Bytes())
	for _, r := range refs {
		w32(int32(len(r.binOrder)))
		for _, bin := range r.binOrder {
			binary.Write(&b, le, bin)
			cs := r.bins[bin]
			w32(int32(len(cs)))
			for _, c := range cs {
				binary.Write(&b, le, c.beg)
				binary.Write(&b, le, c.end)
			}
		}
		w32(int32(len(r.linear)))
		for _, o := range r.linear {
			binary.Write(&b, le, o)
		}
	}

	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	zw := bgzf.NewWriter(f, 1)
	if _, err := zw.Write(b.Bytes()); err != nil {
		zw.Close()
		return err
	}
	return zw.Close()
}
