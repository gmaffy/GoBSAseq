package main

import (
	"fmt"
	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/filter"
	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
	"os"
	"compress/gzip"
)

func main() {
	f, err := os.Open("Pepo.CMV238.ALL.vcf.gz")
	if err != nil {
		panic(err)
	}
	defer f.Close()

	gr, err := gzip.NewReader(f)
	if err != nil {
		panic(err)
	}
	defer gr.Close()

	rdr, err := vcfgo.NewReader(gr, false)
	if err != nil {
		panic(err)
	}

	cfg := utils.AnalysisConfig{
		VCF:             "Pepo.CMV238.ALL.vcf.gz",
		HighParentName:  "SQM238",
		LowParentName:   "SQM097",
		HighBulkName:    "238R_Pool_1",
		Population:      "F2",
		HighParentDepth: 10,
		LowParentDepth:  10,
		HighBulkDepth:   10,
	}

	sampleNames := rdr.Header.SampleNames
	hpIdx, lpIdx, hbIdx, lbIdx := -1, -1, -1, -1
	for i, name := range sampleNames {
		if name == cfg.HighParentName {
			hpIdx = i
		}
		if name == cfg.LowParentName {
			lpIdx = i
		}
		if name == cfg.HighBulkName {
			hbIdx = i
		}
	}
	cfg.HighParentIdx = hpIdx
	cfg.LowParentIdx = lpIdx
	cfg.HighBulkIdx = hbIdx
	cfg.LowBulkIdx = lbIdx
	cfg.Rdr = rdr
	
	fmt.Printf("Indices: hp=%d, lp=%d, hb=%d\n", hpIdx, lpIdx, hbIdx)

	var passed []*vcfgo.Variant
	count := 0
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		count++
		if count%50000 == 0 {
			fmt.Printf("Processed %d variants...\n", count)
		}
		if filter.BsaSeqFilter(v, cfg, "2phb") {
			passed = append(passed, v)
		}
	}

	fmt.Printf("Total variants: %d, Passed: %d\n", count, len(passed))

	rawStats, err := stats.RawStats(cfg, "2phb", []int{hpIdx, lpIdx, hbIdx}, passed)
	if err != nil {
		panic(err)
	}
	fmt.Printf("RawStats generated successfully: %d items\n", len(rawStats))
}
