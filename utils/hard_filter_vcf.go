package utils

import "github.com/brentp/vcfgo"

func HardFilterVcf(rdr *vcfgo.Reader, out string, cfg HardFilterConfig) error {
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		if indelThresholds.readPo

	return nil
}
