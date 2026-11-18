package filter

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/biogo/hts/bgzf"
)

// vcfFixedColumns is the number of columns before the first sample column
// (CHROM POS ID REF ALT QUAL FILTER INFO FORMAT).
const vcfFixedColumns = 9

// openSubsetVCF opens a VCF (plain, gzip or BGZF) and returns a stream in
// which only the given sample columns remain, in the given order. Header
// lines pass through unchanged except #CHROM, which is subset the same way.
//
// Why: vcfgo fully parses every FORMAT field of every sample on each record.
// GoBSAseq uses at most four samples, so on a multisample VCF nearly all of
// the parsing (time and garbage) was spent on samples that were then thrown
// away. Removing them first makes filtering cost independent of how many
// samples the VCF holds.
//
// BGZF input (the norm for bcftools/GATK/DeepVariant output) is decompressed
// on all cores; plain gzip falls back to the standard single-threaded reader.
func openSubsetVCF(path string, keep []int) (io.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}

	src, closeSrc, err := decompressor(f)
	if err != nil {
		f.Close()
		return nil, nil, fmt.Errorf("open %s: %w", path, err)
	}

	pr, pw := io.Pipe()
	go func() {
		pw.CloseWithError(subsetColumns(src, pw, keep))
	}()

	cleanup := func() {
		pr.Close() // unblocks the subsetting goroutine if the reader stops early
		closeSrc()
		f.Close()
	}
	return pr, cleanup, nil
}

// decompressor picks a reader for f by sniffing its first bytes.
func decompressor(f *os.File) (io.Reader, func(), error) {
	br := bufio.NewReaderSize(f, 1<<16)
	head, err := br.Peek(18)
	if err != nil && !errors.Is(err, io.EOF) {
		return nil, nil, err
	}
	isGzip := len(head) >= 2 && head[0] == 0x1f && head[1] == 0x8b
	if !isGzip {
		return br, func() {}, nil
	}
	// BGZF: gzip with FEXTRA set and a "BC" subfield first.
	isBGZF := len(head) >= 14 && head[3]&0x04 != 0 && head[12] == 'B' && head[13] == 'C'
	if isBGZF {
		bg, err := bgzf.NewReader(br, runtime.GOMAXPROCS(0))
		if err != nil {
			return nil, nil, err
		}
		return bg, func() { bg.Close() }, nil
	}
	gz, err := gzip.NewReader(br)
	if err != nil {
		return nil, nil, err
	}
	return gz, func() { gz.Close() }, nil
}

// subsetColumns copies src to dst line by line, keeping the 9 fixed columns
// and only the sample columns listed in keep (0-based sample indices).
func subsetColumns(src io.Reader, dst io.Writer, keep []int) error {
	r := bufio.NewReaderSize(src, 1<<20)
	w := bufio.NewWriterSize(dst, 1<<20)

	// Reused buffers: line accumulator (for lines longer than the reader's
	// buffer) and tab offsets.
	var long []byte
	starts := make([]int, 0, 64)
	ends := make([]int, 0, 64)
	lineNo := 0

	for {
		line, err := r.ReadSlice('\n')
		if errors.Is(err, bufio.ErrBufferFull) {
			long = append(long[:0], line...)
			for errors.Is(err, bufio.ErrBufferFull) {
				line, err = r.ReadSlice('\n')
				long = append(long, line...)
			}
			line = long
		}
		if err != nil && !errors.Is(err, io.EOF) {
			return err
		}
		if len(line) == 0 {
			break
		}
		lineNo++

		// Strip the line terminator (\n or \r\n); we always write \n.
		body := bytes.TrimRight(line, "\r\n")

		if bytes.HasPrefix(body, []byte("##")) {
			w.Write(body)
			w.WriteByte('\n')
		} else if len(body) > 0 {
			// Field boundaries.
			starts, ends = starts[:0], ends[:0]
			s := 0
			for i, c := range body {
				if c == '\t' {
					starts = append(starts, s)
					ends = append(ends, i)
					s = i + 1
				}
			}
			starts = append(starts, s)
			ends = append(ends, len(body))

			if len(starts) <= vcfFixedColumns {
				// Sites-only line: nothing to subset.
				w.Write(body)
			} else {
				w.Write(body[:ends[vcfFixedColumns-1]])
				for _, k := range keep {
					col := vcfFixedColumns + k
					if k < 0 || col >= len(starts) {
						return fmt.Errorf("VCF line %d: sample column %d not present (%d columns)", lineNo, k, len(starts))
					}
					w.WriteByte('\t')
					w.Write(body[starts[col]:ends[col]])
				}
			}
			w.WriteByte('\n')
		}

		if errors.Is(err, io.EOF) {
			break
		}
	}
	return w.Flush()
}
