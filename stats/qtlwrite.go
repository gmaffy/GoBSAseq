package stats

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/xuri/excelize/v2"
)

// QTL table output: every list is written as Excel and as TSV (writeTable).

// WriteIndividualQTLs writes the individual-statistic QTLs to filename
// (.xlsx) and the same path with a .tsv extension.
func WriteIndividualQTLs(qtls []ConsolidatedQTL, bsaType string, filename string) error {
	slots := modeSlots(bsaType)
	header := []string{"Chrom", "Start", "End", "Threshold", "Direction", "N_Statistics", "Lead_Stat", "Lead_Pos"}
	for _, s := range slots {
		header = append(header, slotColumns[s]+"_Pos", slotColumns[s]+"_Value")
	}
	rows := make([][]interface{}, 0, len(qtls))
	for _, q := range qtls {
		row := []interface{}{q.Chrom, bp(q.Start), bp(q.End), levelLabel(q.Threshold), q.Direction, q.NStats, q.LeadStat, q.LeadPos}
		for _, s := range slots {
			if p := q.Peaks[s]; p != nil {
				row = append(row, p.Pos, p.Value)
			} else {
				row = append(row, "", "")
			}
		}
		rows = append(rows, row)
	}
	return writeTable(filename, header, rows)
}

// WriteCompositeZQTLs writes the CompositeZ QTLs to filename (.xlsx) and the
// same path with a .tsv extension.
func WriteCompositeZQTLs(peaks []PeakIntersection, filename string) error {
	header := []string{"Chrom", "Start", "End", "Threshold", "CompositeZ_Pos", "CompositeZ_Value"}
	rows := make([][]interface{}, 0, len(peaks))
	for _, p := range peaks {
		rows = append(rows, []interface{}{p.Chrom, bp(p.Start), bp(p.End), levelLabel(p.ThresholdLevel), p.PeakPos, p.PeakValue})
	}
	return writeTable(filename, header, rows)
}

// finalHeader and finalRow define the final QTL table (Excel and TSV).
var finalHeader = []string{"Chrom", "Start", "Stop", "Span_Start", "Span_Stop", "Consensus", "Support", "Direction",
	"Lead_Method", "Lead_Pos",
	"HMMLOD_Pos", "HMMLOD_Value", "HMMLOD_Threshold", "HMMLOD_Support_Start", "HMMLOD_Support_Stop",
	"BRM_Pos", "BRM_Peak", "BRM_Threshold", "BRM_Z",
	"CompositeZ_Pos", "CompositeZ_Value", "CompositeZ_Threshold",
	"Statistics"}

func finalRow(q FinalQTL) []interface{} {
	cons := "no"
	if q.Consensus {
		cons = "yes"
	}
	row := []interface{}{q.Chrom, q.Start, q.Stop, q.SpanStart, q.SpanStop, cons, q.Support, q.Direction, q.LeadMethod, q.LeadPos}
	peak := func(p *MethodPeak, withZ bool) {
		switch {
		case p == nil && withZ:
			row = append(row, "", "", "", "")
		case p == nil:
			row = append(row, "", "", "")
		case withZ:
			row = append(row, p.Pos, p.Value, p.Threshold, nanOrFloat(p.Z))
		default:
			row = append(row, p.Pos, p.Value, p.Threshold)
		}
	}
	peak(q.HMM, false)
	if q.HMM != nil {
		row = append(row, q.HMMSupportStart, q.HMMSupportStop)
	} else {
		row = append(row, "", "")
	}
	peak(q.BRM, true)
	peak(q.CZ, false)
	return append(row, q.Statistics)
}

// WriteFinalQTLs writes the final QTL table to filename (.xlsx) and to the
// same path with a .tsv extension.
func WriteFinalQTLs(qtls []FinalQTL, filename string) error {
	rows := make([][]interface{}, len(qtls))
	for i, q := range qtls {
		rows[i] = finalRow(q)
	}
	return writeTable(filename, finalHeader, rows)
}

// writeTable writes header and rows to xlsxPath and, as tab-separated text,
// to the same path with a .tsv extension; it returns the first error.
func writeTable(xlsxPath string, header []string, rows [][]interface{}) error {
	if err := writeSheet(xlsxPath, header, rows); err != nil {
		return fmt.Errorf("writing %s: %w", xlsxPath, err)
	}
	tsv := strings.TrimSuffix(xlsxPath, filepath.Ext(xlsxPath)) + ".tsv"
	if err := writeTSV(tsv, header, rows); err != nil {
		return fmt.Errorf("writing %s: %w", tsv, err)
	}
	return nil
}

// writeSheet writes header and rows to Sheet1 of a new workbook at filename
// and returns the first error (cell writes, save or close).
func writeSheet(filename string, header []string, rows [][]interface{}) (err error) {
	f := excelize.NewFile()
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	hdr := make([]interface{}, len(header))
	for i, h := range header {
		hdr[i] = h
	}
	if err := f.SetSheetRow("Sheet1", "A1", &hdr); err != nil {
		return err
	}
	for i := range rows {
		cell, err := excelize.CoordinatesToCellName(1, i+2)
		if err != nil {
			return err
		}
		if err := f.SetSheetRow("Sheet1", cell, &rows[i]); err != nil {
			return err
		}
	}
	return f.SaveAs(filename)
}

// writeTSV writes a header and rows as tab-separated text, returning any
// write, flush or close error.
func writeTSV(path string, header []string, rows [][]interface{}) (err error) {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer func() {
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}()
	w := bufio.NewWriter(f)
	w.WriteString(strings.Join(header, "\t") + "\n")
	for _, r := range rows {
		for i, v := range r {
			if i > 0 {
				w.WriteByte('\t')
			}
			switch x := v.(type) {
			case float64:
				w.WriteString(strconv.FormatFloat(x, 'g', 6, 64))
			default:
				fmt.Fprint(w, x)
			}
		}
		w.WriteByte('\n')
	}
	return w.Flush()
}

// levelLabel is the threshold level as written: p95 regions (used only when
// a statistic has no p99 region in the genome) are suggestive.
func levelLabel(level string) string {
	if level == "p95" {
		return "p95 (suggestive)"
	}
	return level
}

func nanOrFloat(v float64) interface{} {
	if math.IsNaN(v) {
		return ""
	}
	return v
}

// bp is an interpolated position as a whole base pair.
func bp(x float64) int64 { return int64(math.Round(x)) }
