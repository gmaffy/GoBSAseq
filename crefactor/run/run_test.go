package run

import (
	"testing"

	"github.com/gmaffy/GoBSAseq/crefactor/utils"
)

func TestResolveAnalysisModeBulksOnly(t *testing.T) {
	cfg := utils.AnalysisConfig{HighBulkName: "high", LowBulkName: "low"}
	if got := resolveAnalysisMode(cfg, 0, 0, 0, 0); got != modeBulksOnly {
		t.Fatalf("resolveAnalysisMode() = %v, want %v", got, modeBulksOnly)
	}
}

func TestResolveAnalysisModeSingleBulkDoesNotFallThroughToTwoBulk(t *testing.T) {
	cfg := utils.AnalysisConfig{
		HighParentName: "res",
		LowParentName:  "sus",
		OneBulkName:    "bulk",
	}
	if got := resolveAnalysisMode(cfg, 0, 0, 0, 0); got != modeOneBulkLow {
		t.Fatalf("resolveAnalysisMode() = %v, want %v", got, modeOneBulkLow)
	}
}

func TestResolveAnalysisModeRejectsParentsWithoutBulk(t *testing.T) {
	cfg := utils.AnalysisConfig{HighParentName: "res", LowParentName: "sus"}
	if got := resolveAnalysisMode(cfg, 0, 0, 0, 0); got != modeUnsupported {
		t.Fatalf("resolveAnalysisMode() = %v, want %v", got, modeUnsupported)
	}
}
