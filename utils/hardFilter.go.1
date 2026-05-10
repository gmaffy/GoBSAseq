package utils

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/tabix"
	"github.com/brentp/vcfgo"
)

// filterThresholds holds the hard-filter cut-offs for a variant class.
// Every field is a (value, failWhen) pair encoded as a plain float64; the
// sentinel math.NaN() means "not applicable for this class".
type filterThresholds struct {
	qdMin             float64 // QD   < qdMin        → fail
	qualMin           float64 // QUAL < qualMin       → fail
	fsMax             float64 // FS   > fsMax         → fail
	sorMax            float64 // SOR  > sorMax        → fail  (SNP only)
	mqMin             float64 // MQ   < mqMin         → fail  (SNP only)
	mqRankSumMin      float64 // MQRankSum  < min     → fail  (SNP only)
	readPosRankSumMin float64 // ReadPosRankSum < min → fail
}

var (
	snpThresholds = filterThresholds{
		qdMin:             2.0,
		qualMin:           30.0,
		fsMax:             60.0,
		sorMax:            3.0,
		mqMin:             40.0,
		mqRankSumMin:      -12.5,
		readPosRankSumMin: -8.0,
	}

	indelThresholds = filterThresholds{
		qdMin:             2.0,
		qualMin:           30.0,
		fsMax:             200.0,
		sorMax:            -1, // not applied for INDELs (sentinel)
		mqMin:             -1, // not applied for INDELs (sentinel)
		mqRankSumMin:      -1, // not applied for INDELs (sentinel)
		readPosRankSumMin: -20.0,
	}
)

// isINDEL returns true when at least one ALT allele is an insertion or
// deletion (i.e. ref and alt lengths differ).
func isINDEL(v *vcfgo.Variant) bool {
	ref := v.Ref()
	for _, alt := range v.Alt() {
		if len(alt) != len(ref) {
			return true
		}
	}
	return false
}

// getFloat safely retrieves a float-typed INFO field.  vcfgo can return the
// value as float32 or float64 depending on the header declaration, so we
// handle both.
func getFloat(v *vcfgo.Variant, key string) (float64, bool) {
	raw, err := v.Info().Get(key)
	if err != nil || raw == nil {
		return 0, false
	}
	switch val := raw.(type) {
	case float32:
		return float64(val), true
	case float64:
		return val, true
	case int:
		return float64(val), true
	default:
		return 0, false
	}
}

// passesHardFilters applies GATK best-practice hard filters to a single
// variant, choosing SNP or INDEL thresholds automatically.
// Returns true when the variant should be kept (PASS).
func passesHardFilters(v *vcfgo.Variant) bool {
	t := snpThresholds
	if isINDEL(v) {
		t = indelThresholds
	}

	// QUAL is a top-level field on the Variant struct (float32).
	if float64(v.Quality) < t.qualMin {
		return false
	}

	if qd, ok := getFloat(v, "QD"); ok && qd < t.qdMin {
		return false
	}
	if fs, ok := getFloat(v, "FS"); ok && t.fsMax > 0 && fs > t.fsMax {
		return false
	}
	if t.sorMax > 0 {
		if sor, ok := getFloat(v, "SOR"); ok && sor > t.sorMax {
			return false
		}
	}
	if t.mqMin > 0 {
		if mq, ok := getFloat(v, "MQ"); ok && mq < t.mqMin {
			return false
		}
	}
	if t.mqRankSumMin < 0 && t.mqRankSumMin != -1 {
		if mqrs, ok := getFloat(v, "MQRankSum"); ok && mqrs < t.mqRankSumMin {
			return false
		}
	} else if t.mqRankSumMin > 0 { // positive sentinel means check
		if mqrs, ok := getFloat(v, "MQRankSum"); ok && mqrs < t.mqRankSumMin {
			return false
		}
	}
	if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < t.readPosRankSumMin {
		return false
	}

	return true
}

// vcfRecord wraps a VCF line so it satisfies tabix.Record for index building.
type vcfRecord struct {
	chrom string
	start int // 0-based
	end   int // 0-based half-open
}

func (r vcfRecord) RefName() string { return r.chrom }
func (r vcfRecord) Start() int      { return r.start }
func (r vcfRecord) End() int        { return r.end }

type countingWriter struct {
	w io.Writer
	n int64
}

func (w *countingWriter) Write(p []byte) (int, error) {
	n, err := w.w.Write(p)
	w.n += int64(n)
	return n, err
}

// HardFilter reads the input VCF (plain or bgzipped), applies GATK
// best-practice hard filters to every SNP and INDEL, and writes the passing
// variants to outPath as a bgzipped VCF.  A tabix index (.tbi) is written
// alongside the output file.
func HardFilter(vcfPath, outPath string, verbose bool) error {
	// ── open input ────────────────────────────────────────────────────────────
	inFile, err := os.Open(vcfPath)
	if err != nil {
		return fmt.Errorf("opening input: %w", err)
	}
	defer inFile.Close()

	var vcfReader io.Reader = inFile
	if strings.HasSuffix(vcfPath, ".gz") {
		gz, err := gzip.NewReader(inFile)
		if err != nil {
			return fmt.Errorf("opening gzip input: %w", err)
		}
		defer gz.Close()
		vcfReader = gz
	}

	rdr, err := vcfgo.NewReader(vcfReader, false)
	if err != nil {
		return fmt.Errorf("creating vcf reader: %w", err)
	}

	// ── open output (bgzf) ────────────────────────────────────────────────────
	outFile, err := os.Create(outPath)
	if err != nil {
		return fmt.Errorf("creating output file: %w", err)
	}
	defer outFile.Close()

	// bgzf.NewWriter(w, concurrency).  Use 1 for deterministic block boundaries
	// (required for correct tabix offsets).
	countingOut := &countingWriter{w: outFile}
	bgzfWriter := bgzf.NewWriter(countingOut, 1)

	vcfWriter, err := vcfgo.NewWriter(bgzfWriter, rdr.Header)
	if err != nil {
		return fmt.Errorf("creating vcf writer: %w", err)
	}
	if err := bgzfWriter.Flush(); err != nil {
		return fmt.Errorf("flushing vcf header: %w", err)
	}
	if err := bgzfWriter.Wait(); err != nil {
		return fmt.Errorf("writing vcf header: %w", err)
	}

	// ── tabix index ───────────────────────────────────────────────────────────
	tbx := tabix.New()
	tbx.Format = 2
	tbx.NameColumn = 1
	tbx.BeginColumn = 2
	tbx.MetaChar = '#'

	var kept, total int
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		total++

		if !passesHardFilters(v) {
			continue
		}

		// Capture the virtual offset *before* writing this record.
		startOffset := bgzf.Offset{File: countingOut.n}

		vcfWriter.WriteVariant(v)

		if err := bgzfWriter.Flush(); err != nil {
			return fmt.Errorf("flushing variant: %w", err)
		}
		if err := bgzfWriter.Wait(); err != nil {
			return fmt.Errorf("writing variant: %w", err)
		}
		endOffset := bgzf.Offset{File: countingOut.n}

		// tabix uses 0-based half-open coordinates.
		rec := vcfRecord{
			chrom: v.Chromosome,
			start: int(v.Pos) - 1,                // VCF POS is 1-based
			end:   int(v.Pos) - 1 + len(v.Ref()), // approximate; correct for SNPs & simple indels
		}

		if err := tbx.Add(rec, bgzf.Chunk{Begin: startOffset, End: endOffset}, true, true); err != nil {
			if verbose {
				fmt.Fprintf(os.Stderr, "tabix add warning for %s:%d : %v\n", rec.chrom, v.Pos, err)
			}
		}

		kept++
		if verbose {
			fmt.Fprintf(os.Stderr, "kept %s:%d\n", v.Chromosome, v.Pos)
		}
	}

	if err := rdr.Error(); err != nil && err != io.EOF {
		return fmt.Errorf("vcf read error: %w", err)
	}

	if err := bgzfWriter.Close(); err != nil {
		return fmt.Errorf("closing bgzf writer: %w", err)
	}

	fmt.Printf("%s created (%d/%d variants passed hard filters)\n", outPath, kept, total)

	// ── write tabix index ─────────────────────────────────────────────────────
	tbiPath := outPath + ".tbi"
	tbiFile, err := os.Create(tbiPath)
	if err != nil {
		return fmt.Errorf("creating tbi file: %w", err)
	}
	defer tbiFile.Close()

	// The tabix spec requires the .tbi itself to be bgzipped.
	tbiGz := bgzf.NewWriter(tbiFile, 1)
	if err := tabix.WriteTo(tbiGz, tbx); err != nil {
		return fmt.Errorf("writing tabix index: %w", err)
	}
	if err := tbiGz.Close(); err != nil {
		return fmt.Errorf("closing tbi bgzf: %w", err)
	}

	fmt.Printf("%s created\n", tbiPath)
	return nil
}
