package utils

import (
	"fmt"
	"os"

	"github.com/biogo/hts/bgzf"
	"github.com/brentp/vcfgo"
)

func isSNP(v *vcfgo.Variant) bool {
	if len(v.Ref()) != 1 {
		return false
	}
	for _, alt := range v.Alt() {
		if len(alt) != 1 || alt == "." || alt == "*" {
			return false
		}
	}
	return true
}

func isIndel(v *vcfgo.Variant) bool {
	refLen := len(v.Ref())
	for _, alt := range v.Alt() {
		if alt == "." || alt == "*" {
			continue
		}
		if len(alt) != refLen {
			return true
		}
	}
	return false
}

func failsSnpHardFilters(v *vcfgo.Variant, cfg HardFilterConfig) bool {
	if float64(v.Quality) < cfg.SNP_QUAL_Min {
		return true
	}
	if qd, ok := getFloat(v, "QD"); ok && qd < cfg.SNP_QD_Min {
		return true
	}
	if fs, ok := getFloat(v, "FS"); ok && fs > cfg.SNP_FS_Max {
		return true
	}
	if sor, ok := getFloat(v, "SOR"); ok && sor > cfg.SNP_SOR_Max {
		return true
	}
	if mq, ok := getFloat(v, "MQ"); ok && mq < cfg.SNP_MQ_Min {
		return true
	}
	if mqRankSum, ok := getFloat(v, "MQRankSum"); ok && mqRankSum < cfg.SNP_MQRankSum_Min {
		return true
	}
	if readPosRankSum, ok := getFloat(v, "ReadPosRankSum"); ok && readPosRankSum < cfg.SNP_ReadPosRankSum_Min {
		return true
	}
	return false
}

func failsIndelHardFilters(v *vcfgo.Variant, cfg HardFilterConfig) bool {
	if float64(v.Quality) < cfg.INDEL_QUAL_Min {
		return true
	}
	if qd, ok := getFloat(v, "QD"); ok && qd < cfg.INDEL_QD_Min {
		return true
	}
	if fs, ok := getFloat(v, "FS"); ok && fs > cfg.INDEL_FS_Max {
		return true
	}
	if readPosRankSum, ok := getFloat(v, "ReadPosRankSum"); ok && readPosRankSum < cfg.INDEL_ReadPosRankSum_Min {
		return true
	}
	return false
}

func HardFilterVcf(rdr *vcfgo.Reader, out string, cfg HardFilterConfig) ([]*vcfgo.Variant, error) {
	f, err := os.Create(out)
	if err != nil {
		return nil, fmt.Errorf("failed to create file: %w", err)
	}
	defer f.Close()

	bgzfWriter := bgzf.NewWriter(f, 1)
	writer, err := vcfgo.NewWriter(bgzfWriter, rdr.Header)
	if err != nil {
		bgzfWriter.Close()
		return nil, fmt.Errorf("failed to create VCF writer: %w", err)
	}

	var hardFilteredVariants []*vcfgo.Variant
	for {
		v := rdr.Read()
		if v == nil {
			break
		}

		keep := false
		if isSNP(v) {
			keep = !failsSnpHardFilters(v, cfg)
		} else if isIndel(v) {
			keep = !failsIndelHardFilters(v, cfg)
		}
		if !keep {
			continue
		}

		hardFilteredVariants = append(hardFilteredVariants, v)
		writer.WriteVariant(v)
	}

	if err := rdr.Error(); err != nil {
		bgzfWriter.Close()
		return nil, fmt.Errorf("failed to read VCF: %w", err)
	}
	if err := bgzfWriter.Close(); err != nil {
		return nil, fmt.Errorf("failed to close VCF writer: %w", err)
	}

	return hardFilteredVariants, nil
}
