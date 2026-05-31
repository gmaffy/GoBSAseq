package utils

import (
	"fmt"
	"os"
	"path/filepath"
	"time"
)

// PrepareRunDir creates a timestamped subdirectory under baseOut for this analysis run.
// Returns the path to use as cfg.OutputDir.
func PrepareRunDir(baseOut string) (string, error) {
	root := filepath.Join(baseOut, "GoBSAseqResults")
	if err := os.MkdirAll(root, 0o755); err != nil {
		return "", fmt.Errorf("create results root: %w", err)
	}
	now := time.Now()
	runDir := filepath.Join(root, fmt.Sprintf("%02d_%02d_%04d_%02d_%02d_%02d",
		now.Day(), now.Month(), now.Year(), now.Hour(), now.Minute(), now.Second()))
	if err := os.MkdirAll(runDir, 0o755); err != nil {
		return "", fmt.Errorf("create run directory: %w", err)
	}
	return runDir, nil
}
