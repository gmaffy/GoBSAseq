package filter

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bgzf"
)

// This file builds a tabix (.tbi) index for a BGZF-compressed VCF.
//
// It deliberately does not use github.com/biogo/hts/tabix: in v1.4.5,
// tabix.Index.Add never records a new reference name in its name map, so every
// record is registered as a brand-new reference and the resulting index is
// unusable for region queries (n_ref == number of records).

const (
	tbiMinShift   = 14 // 16 kb linear-index windows
	tbiMaxBinMeta = 37450
)

type tbiChunk struct{ beg, end uint64 }

type tbiRef struct {
	name    string
	bins    map[uint32][]tbiChunk
	binList []uint32 // insertion order, for deterministic output
	linear  []uint64
	hasLin  []bool
}

// vOffset converts a bgzf.Offset to a 64-bit BGZF virtual offset.
func vOffset(o bgzf.Offset) uint64 {
	return uint64(o.File)<<16 | uint64(o.Block)
}

// reg2bin computes the UCSC/SAM binning-scheme bin for the 0-based half-open
// interval [beg, end).
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

func (r *tbiRef) add(beg, end int, c tbiChunk) {
	if end <= beg {
		end = beg + 1
	}
	bin := reg2bin(beg, end)
	chunks, ok := r.bins[bin]
	if !ok {
		r.binList = append(r.binList, bin)
	}
	// Records arrive in file order, so extend the last chunk when contiguous.
	if n := len(chunks); n > 0 && chunks[n-1].end >= c.beg {
		if c.end > chunks[n-1].end {
			chunks[n-1].end = c.end
		}
	} else {
		chunks = append(chunks, c)
	}
	r.bins[bin] = chunks

	for w := beg >> tbiMinShift; w <= (end-1)>>tbiMinShift; w++ {
		for len(r.linear) <= w {
			r.linear = append(r.linear, 0)
			r.hasLin = append(r.hasLin, false)
		}
		if !r.hasLin[w] || c.beg < r.linear[w] {
			r.linear[w] = c.beg
			r.hasLin[w] = true
		}
	}
}

// buildTabixIndex reads a finished BGZF-compressed VCF back, records the exact
// virtual-offset chunk of every data line, and writes a BGZF-compressed .tbi
// next to it. Reading the closed file avoids depending on the writer's
// (asynchronous) compressed byte position.
func buildTabixIndex(vcfGzPath string) error {
	f, err := os.Open(vcfGzPath)
	if err != nil {
		return err
	}
	defer f.Close()

	r, err := bgzf.NewReader(f, 1)
	if err != nil {
		return err
	}
	defer r.Close()

	var (
		refs    []*tbiRef
		byName  = map[string]*tbiRef{}
		line    []byte
		begin   bgzf.Offset
		started bool
		lineNo  int
	)
	for {
		b, rerr := r.ReadByte()
		if rerr != nil && rerr != io.EOF {
			return rerr
		}
		if rerr == io.EOF && !started {
			break
		}
		if rerr == nil {
			if !started {
				begin = r.LastChunk().Begin
				started = true
			}
			if b != '\n' {
				line = append(line, b)
				continue
			}
		}

		// End of a line (or a final unterminated line at EOF).
		lineNo++
		end := r.LastChunk().End
		if len(line) > 0 && line[0] != '#' {
			fields := strings.SplitN(string(line), "\t", 5)
			if len(fields) < 4 {
				return fmt.Errorf("index %s: malformed VCF line %d", vcfGzPath, lineNo)
			}
			pos, perr := strconv.Atoi(fields[1])
			if perr != nil || pos < 1 {
				return fmt.Errorf("index %s: bad POS on line %d", vcfGzPath, lineNo)
			}
			ref, ok := byName[fields[0]]
			if !ok {
				ref = &tbiRef{name: fields[0], bins: map[uint32][]tbiChunk{}}
				byName[fields[0]] = ref
				refs = append(refs, ref)
			} else if ref != refs[len(refs)-1] {
				return fmt.Errorf("index %s: contig %s is not contiguous (line %d); file must be sorted", vcfGzPath, fields[0], lineNo)
			}
			ref.add(pos-1, pos-1+len(fields[3]), tbiChunk{vOffset(begin), vOffset(end)})
		}
		line = line[:0]
		started = false
		if rerr == io.EOF {
			break
		}
	}

	out, err := os.Create(vcfGzPath + ".tbi")
	if err != nil {
		return err
	}
	defer out.Close()
	gz := bgzf.NewWriter(out, 1)
	bw := bufio.NewWriter(gz)
	if err := writeTBI(bw, refs); err != nil {
		gz.Close()
		return err
	}
	if err := bw.Flush(); err != nil {
		gz.Close()
		return err
	}
	if err := gz.Close(); err != nil {
		return err
	}
	return out.Close()
}

// writeTBI serialises the index in the tabix binary format (uncompressed; the
// caller wraps it in BGZF).
func writeTBI(w io.Writer, refs []*tbiRef) error {
	le := binary.LittleEndian
	put := func(v any) error { return binary.Write(w, le, v) }

	var names []byte
	for _, r := range refs {
		names = append(names, r.name...)
		names = append(names, 0)
	}
	header := []any{
		[4]byte{'T', 'B', 'I', 1},
		int32(len(refs)),
		int32(2), // format: VCF
		int32(1), // col_seq
		int32(2), // col_beg
		int32(0), // col_end
		int32('#'),
		int32(0), // skip
		int32(len(names)),
	}
	for _, v := range header {
		if err := put(v); err != nil {
			return err
		}
	}
	if _, err := w.Write(names); err != nil {
		return err
	}

	for _, r := range refs {
		if err := put(int32(len(r.binList))); err != nil {
			return err
		}
		for _, bin := range r.binList {
			if bin == tbiMaxBinMeta {
				continue // reserved pseudo-bin; unreachable for real coordinates
			}
			chunks := r.bins[bin]
			if err := put(bin); err != nil {
				return err
			}
			if err := put(int32(len(chunks))); err != nil {
				return err
			}
			for _, c := range chunks {
				if err := put(c.beg); err != nil {
					return err
				}
				if err := put(c.end); err != nil {
					return err
				}
			}
		}
		// Linear index: windows with no record inherit the previous offset.
		var last uint64
		for i := range r.linear {
			if r.hasLin[i] {
				last = r.linear[i]
			} else {
				r.linear[i] = last
			}
		}
		if err := put(int32(len(r.linear))); err != nil {
			return err
		}
		for _, off := range r.linear {
			if err := put(off); err != nil {
				return err
			}
		}
	}
	return nil
}
