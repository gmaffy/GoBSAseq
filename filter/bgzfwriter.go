package filter

import (
	"bytes"
	"compress/flate"
	"encoding/binary"
	"hash/crc32"
	"io"
	"runtime"
	"sync"
)

// bgzfMaxBlock is the largest uncompressed payload per BGZF block (htslib's
// BGZF_BLOCK_SIZE, 0xff00), leaving room for deflate overhead in 64 KiB.
const bgzfMaxBlock = 0xff00

// bgzfEOF is the standard empty BGZF block that marks end-of-file.
var bgzfEOF = []byte{
	0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
	0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00,
}

// vpos is a BGZF virtual position expressed as (block number, offset in block).
// The block number becomes a compressed file offset once the block is written,
// which is what lets index entries be exact even though compression is async.
type vpos struct {
	block uint32
	off   uint16
}

type bgzfBlock struct {
	data []byte
	out  bytes.Buffer
	done chan struct{}
}

// bgzfWriter compresses BGZF blocks on several goroutines and writes them in
// order. Only one goroutine may call Write/Pos/Close.
type bgzfWriter struct {
	w       io.Writer
	level   int
	cur     []byte
	blockNo uint32

	jobs    chan *bgzfBlock
	ordered chan *bgzfBlock
	pool    sync.Pool
	wg      sync.WaitGroup
	outDone chan struct{}

	offsets []int64 // compressed file offset of each block (filled by the output goroutine)
	written int64
	err     error
}

func newBGZFWriter(w io.Writer, level int) *bgzfWriter {
	workers := runtime.GOMAXPROCS(0)
	bw := &bgzfWriter{
		w:       w,
		level:   level,
		cur:     make([]byte, 0, bgzfMaxBlock),
		jobs:    make(chan *bgzfBlock, workers*2),
		ordered: make(chan *bgzfBlock, workers*4),
		outDone: make(chan struct{}),
	}
	bw.pool.New = func() any { return make([]byte, 0, bgzfMaxBlock) }
	for range workers {
		bw.wg.Add(1)
		go bw.compressLoop()
	}
	go bw.outputLoop()
	return bw
}

func (bw *bgzfWriter) compressLoop() {
	defer bw.wg.Done()
	var fw *flate.Writer
	for b := range bw.jobs {
		b.out.Grow(len(b.data) + 64)
		hdr := [18]byte{0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0}
		b.out.Write(hdr[:])
		level := bw.level
		for {
			if fw == nil {
				fw, _ = flate.NewWriter(&b.out, level)
			} else {
				fw.Reset(&b.out)
			}
			fw.Write(b.data)
			fw.Close()
			// A BGZF block must fit in 64 KiB; incompressible data at a
			// compressing level can overshoot, so store it instead.
			if b.out.Len()+8 <= 1<<16 || level == flate.NoCompression {
				break
			}
			b.out.Truncate(18)
			level = flate.NoCompression
			fw = nil
		}
		if level != bw.level {
			fw = nil
		}
		var tail [8]byte
		binary.LittleEndian.PutUint32(tail[0:], crc32.ChecksumIEEE(b.data))
		binary.LittleEndian.PutUint32(tail[4:], uint32(len(b.data)))
		b.out.Write(tail[:])
		raw := b.out.Bytes()
		binary.LittleEndian.PutUint16(raw[16:], uint16(len(raw)-1)) // BSIZE
		close(b.done)
	}
}

func (bw *bgzfWriter) outputLoop() {
	defer close(bw.outDone)
	for b := range bw.ordered {
		<-b.done
		bw.offsets = append(bw.offsets, bw.written)
		if bw.err == nil {
			n, err := bw.w.Write(b.out.Bytes())
			bw.written += int64(n)
			bw.err = err
		}
		bw.pool.Put(b.data[:0])
	}
}

func (bw *bgzfWriter) flushBlock() {
	if len(bw.cur) == 0 {
		return
	}
	b := &bgzfBlock{data: bw.cur, done: make(chan struct{})}
	bw.ordered <- b
	bw.jobs <- b
	bw.cur = bw.pool.Get().([]byte)
	bw.blockNo++
}

// Pos returns the virtual position of the next byte to be written.
func (bw *bgzfWriter) Pos() vpos { return vpos{bw.blockNo, uint16(len(bw.cur))} }

// Reserve starts a new block if a record of n bytes would not fit in the
// current one. Call it before Pos so a record's start offset is the block it is
// actually written to.
func (bw *bgzfWriter) Reserve(n int) {
	if len(bw.cur) > 0 && len(bw.cur)+n > bgzfMaxBlock && n <= bgzfMaxBlock {
		bw.flushBlock()
	}
}

// Write buffers p, cutting blocks at bgzfMaxBlock. A record that does not fit
// in the current block starts a new one (if it fits in an empty block), so most
// records live in a single block.
func (bw *bgzfWriter) Write(p []byte) (int, error) {
	n := len(p)
	bw.Reserve(len(p))
	for len(p) > 0 {
		room := bgzfMaxBlock - len(bw.cur)
		if room > len(p) {
			room = len(p)
		}
		bw.cur = append(bw.cur, p[:room]...)
		p = p[room:]
		if len(bw.cur) == bgzfMaxBlock {
			bw.flushBlock()
		}
	}
	return n, nil
}

// Close flushes the last block, writes the EOF marker and waits for all output.
// After Close, FileOffset resolves virtual positions.
func (bw *bgzfWriter) Close() error {
	bw.flushBlock()
	close(bw.ordered)
	close(bw.jobs)
	bw.wg.Wait()
	<-bw.outDone
	if bw.err != nil {
		return bw.err
	}
	// A virtual position pointing just past the last block refers to the EOF block.
	bw.offsets = append(bw.offsets, bw.written)
	_, err := bw.w.Write(bgzfEOF)
	return err
}

// virtualOffset converts a vpos to a BGZF virtual file offset (valid after Close).
func (bw *bgzfWriter) virtualOffset(p vpos) (file int64, block uint16) {
	return bw.offsets[p.block], p.off
}
