package filter

import (
	"github.com/brentp/vcfgo"
	"github.com/gmaffy/GoBSAseq/utils"
)

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

func classifyVariant(v *vcfgo.Variant) (isSNP, isIndel bool) {
	refLen := len(v.Ref())
	isSNP = refLen == 1
	for _, alt := range v.Alt() {
		if alt == "." || alt == "*" || (len(alt) > 0 && alt[0] == '<') {
			continue
		}
		if len(alt) != 1 {
			isSNP = false
		}
		if len(alt) != refLen {
			isIndel = true
		}
	}
	return
}

func PassesHardFilter(v *vcfgo.Variant, hfcfg utils.HardFilterConfig) bool {
	isSNP, isIndel := classifyVariant(v)

	switch {
	case isSNP:
		if float64(v.Quality) < hfcfg.SNP_QUAL_Min {
			return false
		}
		if qd, ok := getFloat(v, "QD"); ok && qd < hfcfg.SNP_QD_Min {
			return false
		}
		if fs, ok := getFloat(v, "FS"); ok && fs > hfcfg.SNP_FS_Max {
			return false
		}
		if sor, ok := getFloat(v, "SOR"); ok && sor > hfcfg.SNP_SOR_Max {
			return false
		}
		if mq, ok := getFloat(v, "MQ"); ok && mq < hfcfg.SNP_MQ_Min {
			return false
		}
		if mqrs, ok := getFloat(v, "MQRankSum"); ok && mqrs < hfcfg.SNP_MQRankSum_Min {
			return false
		}
		if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.SNP_ReadPosRankSum_Min {
			return false
		}
		return true

	case isIndel:
		if float64(v.Quality) < hfcfg.INDEL_QUAL_Min {
			return false
		}
		if qd, ok := getFloat(v, "QD"); ok && qd < hfcfg.INDEL_QD_Min {
			return false
		}
		if fs, ok := getFloat(v, "FS"); ok && fs > hfcfg.INDEL_FS_Max {
			return false
		}
		if sor, ok := getFloat(v, "SOR"); ok && sor > hfcfg.INDEL_SOR_Max {
			return false
		}
		if rprs, ok := getFloat(v, "ReadPosRankSum"); ok && rprs < hfcfg.INDEL_ReadPosRankSum_Min {
			return false
		}
		return true

	default:
		return false // mixed / MNP / symbolic-only
	}
}
