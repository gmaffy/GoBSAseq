package utils

import (
	"fmt"
	"os"
	"path/filepath"
)

func CreateResultsDir(outputDir string) (string, error) {

	baseDir := filepath.Join(outputDir, "goBSAseqResults")
	//if err := os.MkdirAll(baseDir, 0755); err != nil {
	//	return "", fmt.Errorf("error creating results directory: %w", err)
	//}
	//
	//now := time.Now()
	//resultsDir := filepath.Join(baseDir, fmt.Sprintf("%02d_%02d_%04d_%02d_%02d_%02d", now.Day(), now.Month(), now.Year(), now.Hour(), now.Minute(), now.Second()))
	resultsDir, _ := filepath.Abs(baseDir)
	err := os.MkdirAll(resultsDir, 0755)
	if err != nil {
		return "", fmt.Errorf("error creating results directory: %w", err)
	}
	fmt.Printf("Created results directory at %s ..\n\n", resultsDir)

	return resultsDir, nil
}
