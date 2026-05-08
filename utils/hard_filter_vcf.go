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
	for {
		v := rdr.Read()
		if v == nil {
			break
		}
		if isSNP(v) {
			if v.Qual() < cfg.MinSnpQual {
				v.AddFilter("LowQual")
			}


	return nil
}
