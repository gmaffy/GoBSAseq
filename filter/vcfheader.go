package filter

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"runtime"
	"strings"

	"github.com/biogo/hts/bgzf"
)

// FieldDef is one ##INFO or ##FORMAT header definition.
type FieldDef struct {
	Number string // "1", "A", "R", "G", "." …
	Type   string // "Integer", "Float", "String", "Flag" …
}

// Header is the part of a VCF header the filter needs. Meta keeps every "##"
// line verbatim so the filtered VCF carries the original header.
type Header struct {
	Meta    []string // "##…" lines, without trailing newline
	Samples []string // sample names from the #CHROM line
	Infos   map[string]FieldDef
	Formats map[string]FieldDef
}

// ParseHeader reads header lines from br up to and including #CHROM, leaving br
// positioned at the first record.
func ParseHeader(br *bufio.Reader) (*Header, error) {
	h := &Header{Infos: map[string]FieldDef{}, Formats: map[string]FieldDef{}}
	for {
		line, err := br.ReadString('\n')
		line = strings.TrimRight(line, "\r\n")
		switch {
		case strings.HasPrefix(line, "##"):
			h.Meta = append(h.Meta, line)
			if strings.HasPrefix(line, "##INFO=<") {
				if id, def, ok := parseFieldDef(line[len("##INFO=<"):]); ok {
					h.Infos[id] = def
				}
			} else if strings.HasPrefix(line, "##FORMAT=<") {
				if id, def, ok := parseFieldDef(line[len("##FORMAT=<"):]); ok {
					h.Formats[id] = def
				}
			}
		case strings.HasPrefix(line, "#CHROM"):
			cols := strings.Split(line, "\t")
			if len(cols) > 9 {
				h.Samples = cols[9:]
			}
			return h, nil
		case line != "":
			return nil, fmt.Errorf("VCF header has no #CHROM line before the first record")
		}
		if err != nil {
			if err == io.EOF {
				return nil, fmt.Errorf("VCF ended before the #CHROM header line")
			}
			return nil, err
		}
	}
}

// parseFieldDef extracts ID, Number and Type from the inside of an
// ##INFO=<…> / ##FORMAT=<…> line. Description may contain commas inside quotes.
func parseFieldDef(body string) (string, FieldDef, bool) {
	var id string
	var def FieldDef
	inQuote := false
	start := 0
	take := func(kv string) {
		eq := strings.IndexByte(kv, '=')
		if eq < 0 {
			return
		}
		switch kv[:eq] {
		case "ID":
			id = kv[eq+1:]
		case "Number":
			def.Number = kv[eq+1:]
		case "Type":
			def.Type = kv[eq+1:]
		}
	}
	for i := 0; i < len(body); i++ {
		switch body[i] {
		case '"':
			inQuote = !inQuote
		case ',', '>':
			if !inQuote {
				take(body[start:i])
				start = i + 1
			}
		}
	}
	return id, def, id != ""
}

// HasInfo reports whether an INFO key is declared in the header.
func (h *Header) HasInfo(id string) bool { _, ok := h.Infos[id]; return ok }

// openVCFStream opens a .vcf, .vcf.gz (bgzip) or plain-gzip VCF. BGZF input is
// decompressed by several goroutines; plain gzip falls back to compress/gzip,
// which is single-threaded (bgzip the file to make it faster).
func openVCFStream(path string) (*bufio.Reader, func(), error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	peek := make([]byte, 18)
	n, _ := io.ReadFull(f, peek)
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		f.Close()
		return nil, nil, err
	}
	peek = peek[:n]

	isGzip := n >= 2 && peek[0] == 0x1f && peek[1] == 0x8b
	isBGZF := isGzip && n >= 16 && peek[3]&4 != 0 && peek[12] == 'B' && peek[13] == 'C'

	const bufSize = 4 << 20
	switch {
	case isBGZF:
		threads := runtime.GOMAXPROCS(0)
		if threads > 8 {
			threads = 8 // decompression saturates well before this
		}
		bz, err := bgzf.NewReader(f, threads)
		if err != nil {
			f.Close()
			return nil, nil, fmt.Errorf("open bgzf: %w", err)
		}
		return bufio.NewReaderSize(bz, bufSize), func() { bz.Close(); f.Close() }, nil
	case isGzip:
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, nil, fmt.Errorf("open gzip: %w", err)
		}
		return bufio.NewReaderSize(gz, bufSize), func() { gz.Close(); f.Close() }, nil
	default:
		return bufio.NewReaderSize(f, bufSize), func() { f.Close() }, nil
	}
}

// ReadVCFHeader reads only the header of a VCF (used for sample selection).
func ReadVCFHeader(path string) (*Header, error) {
	br, closer, err := openVCFStream(path)
	if err != nil {
		return nil, err
	}
	defer closer()
	return ParseHeader(br)
}

// readLineInto appends the next line (without '\n') to arena and returns the
// extended arena and the line's start offset. It reuses the bufio buffer via
// ReadSlice, so a line costs one copy into the batch arena and no allocation.
func readLineInto(br *bufio.Reader, arena []byte) ([]byte, int, error) {
	start := len(arena)
	for {
		chunk, err := br.ReadSlice('\n')
		arena = append(arena, chunk...)
		if err == bufio.ErrBufferFull {
			continue
		}
		if n := len(arena); n > start && arena[n-1] == '\n' {
			arena = arena[:n-1]
			if n := len(arena); n > start && arena[n-1] == '\r' {
				arena = arena[:n-1]
			}
		}
		return arena, start, err
	}
}

// isBGZF reports whether data starts with a BGZF block header (used by tests).
func isBGZF(data []byte) bool {
	return len(data) >= 16 && bytes.Equal(data[:2], []byte{0x1f, 0x8b}) && data[3]&4 != 0 && data[12] == 'B' && data[13] == 'C'
}
