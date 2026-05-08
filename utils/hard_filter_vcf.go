package utils

import "github.com/brentp/vcfgo"

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
            return true  // insertion or deletion
        }
    }
    return false
}

func HardFilterVcf(rdr *vcfgo.Reader, out string, cfg HardFilterConfig) error {
	// 1. Create the output file
	f, err := os.Create(outPath)
	if err != nil {
		return fmt.Errorf("failed to create file: %w", err)
	}
	defer f.Close()

	// 2. Wrap in gzip writer for .vcf.gz
	// Note: For Tabix compatibility, BGZF is preferred over standard Gzip.
	// If you don't have a BGZF library, standard gzip works with newer tabix versions.
	bgzfWriter := bgzf.NewWriter(f, 1) // 1 is compression level
	defer bgzfWriter.Close()


	writer, err := vcfgo.NewWriter(gw, rdr.Header())
	if err != nil {
		return fmt.Errorf("failed to create VCF writer: %w", err)
	}

	var hardFilteredVariants []*vcfgo.Variant
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		if isSNP(v) {
			if v.Qual() < cfg.MinSnpQual &&
			v.Info().Get("QD").(float64) < cfg.MinSnpQD &&
			v.Info().Get("FS").(float64) > cfg.MaxSnpFS &&
			v.Info().Get("SOR").(float64) > cfg.MaxSnpSOR &&
			v.Info().Get("MQ").(float64) < cfg.MinSnpMQ &&
			v.Info().Get("MQRankSum").(float64) < cfg.MinSnpMQRankSum &&
			v.Info().Get("ReadPosRankSum").(float64) < cfg.MinSnpReadPosRankSum {
				hardFilteredVariants = append(hardFilteredVariants, v)

		} else if isIndel(v) {
			if v.Qual() < cfg.MinIndelQual &&
			v.Info().Get("QD").(float64) < cfg.MinIndelQD &&
			v.Info().Get("FS").(float64) > cfg.MaxIndelFS &&
			v.Info().Get("ReadPosRankSum").(float64) < cfg.MinIndelReadPosRankSum {
				hardFilteredVariants = append(hardFilteredVariants, v)
			}
		}


	return nil
}
