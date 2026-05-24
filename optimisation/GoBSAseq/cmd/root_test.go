package cmd

import (
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestReadVCFSamples(t *testing.T) {
	dir := t.TempDir()
	vcf := filepath.Join(dir, "input.vcf")
	data := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tB1\tB2\n"
	if err := os.WriteFile(vcf, []byte(data), 0o644); err != nil {
		t.Fatal(err)
	}

	got, err := readVCFSamples(vcf)
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"P1", "P2", "B1", "B2"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("readVCFSamples() = %#v, want %#v", got, want)
	}
}

func TestSplitCSVSkipsNoneAndEmpty(t *testing.T) {
	got := splitCSV("ParentA, none, ,ParentB")
	want := []string{"ParentA", "ParentB"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("splitCSV() = %#v, want %#v", got, want)
	}
}
