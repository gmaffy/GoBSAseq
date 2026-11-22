package stats

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/gob"
	"encoding/hex"
	"hash"
	"math"
	"os"
	"path/filepath"

	"github.com/gmaffy/GoBSAseq/utils"
)

// Cache of the genome-wide null (CalculateThresholds). The null depends only
// on the raw statistics' markers, depths and per-site signals (which set ρ)
// and on the settings below, so a re-run that changes only QTL calling,
// plotting or gene-space options can reuse it.

// nullCacheVersion changes whenever the null simulation or anything it
// computes changes, so an older cache is never reused.
const nullCacheVersion = 1

type nullCacheFile struct {
	Version int
	Key     string
	Th      Thresholds
}

// nullCacheKey hashes everything genomeNullThresholds reads.
func nullCacheKey(cfg utils.AnalysisConfig, bsaType string, raw []BSAstats) string {
	h := sha256.New()
	str := func(s string) { writeU64(h, uint64(len(s))); h.Write([]byte(s)) }
	i64 := func(v int64) { writeU64(h, uint64(v)) }
	f64 := func(v float64) { writeU64(h, math.Float64bits(v)) }
	b := func(v bool) {
		if v {
			i64(1)
		} else {
			i64(0)
		}
	}
	i64(nullCacheVersion)
	str(bsaType)
	str(cfg.Population)
	i64(int64(cfg.HighBulkSize))
	i64(int64(cfg.LowBulkSize))
	i64(int64(cfg.WindowSize))
	i64(int64(cfg.StepSize))
	f64(cfg.CMPerMb)
	i64(int64(cfg.NullReps))
	f64(cfg.BrmAlpha)
	b(cfg.SignedDeltaSI(bsaType))
	f64(hmmBinMorgans)
	i64(nullSeed)
	i64(int64(len(raw)))
	for i := range raw {
		s := &raw[i]
		str(s.CHROM)
		i64(s.POS)
		b(s.HighIsAlt)
		i64(int64(s.HighDepth))
		i64(int64(s.LowDepth))
		i64(int64(s.Depth))
		f64(s.HighSI)
		f64(s.LowSI)
		f64(s.OneBulkP0)
		f64(s.OneBulkAFDev)
	}
	return hex.EncodeToString(h.Sum(nil))
}

func writeU64(h hash.Hash, v uint64) {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], v)
	h.Write(buf[:])
}

// loadNullCache returns the cached thresholds at path if they were computed
// for key.
func loadNullCache(path, key string) (Thresholds, bool) {
	f, err := os.Open(path)
	if err != nil {
		return Thresholds{}, false
	}
	defer f.Close()
	var c nullCacheFile
	if err := gob.NewDecoder(f).Decode(&c); err != nil || c.Version != nullCacheVersion || c.Key != key {
		return Thresholds{}, false
	}
	return c.Th, true
}

// saveNullCache writes th for key to path (via a temporary file, so an
// interrupted write never leaves a truncated cache).
func saveNullCache(path, key string, th Thresholds) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(path), filepath.Base(path)+".*.tmp")
	if err != nil {
		return err
	}
	defer os.Remove(tmp.Name())
	if err := gob.NewEncoder(tmp).Encode(nullCacheFile{nullCacheVersion, key, th}); err != nil {
		tmp.Close()
		return err
	}
	if err := tmp.Close(); err != nil {
		return err
	}
	return os.Rename(tmp.Name(), path)
}
