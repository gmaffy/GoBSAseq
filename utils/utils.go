package utils

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"time"

	"github.com/brentp/vcfgo"
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

func CreateResultsDir(outputDir string) (string, error) {

	baseDir := filepath.Join(outputDir, "goBSAseqResults")
	bErr := os.MkdirAll(filepath.Join(outputDir, "goBSAseqResults"), 0755)
	if bErr != nil {
		log.Fatalf("Error creating results directory: %s\n", bErr)
		return "", bErr
	}

	now := time.Now()
	resultsDir := filepath.Join(baseDir, fmt.Sprintf("%02d_%02d_%04d_%02d_%02d_%02d", now.Day(), now.Month(), now.Year(), now.Hour(), now.Minute(), now.Second()))

	err := os.MkdirAll(resultsDir, 0755)
	if err != nil {
		log.Fatalf("Error creating results directory: %s\n", err)
		return "", err

	}
	fmt.Printf("Created results directory at %s ..\n\n", resultsDir)

	return resultsDir, nil
}
