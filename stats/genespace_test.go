package stats

import (
	"testing"

	"github.com/gmaffy/GoBSAseq/utils"
)

func TestGeneSpaceGTLines(t *testing.T) {
	cfg := utils.AnalysisConfig{
		HighParentName: "R_parent",
		LowParentName:  "S_parent",
		HighBulkName:   "R_bulk",
		LowBulkName:    "S_bulk",
	}

	res, sus := geneSpaceGTLines(cfg, "2p2b")
	if len(res) != 2 || res[0] != "R_parent.GT" || res[1] != "R_bulk.GT" {
		t.Fatalf("2p2b res = %v, want [R_parent.GT R_bulk.GT]", res)
	}
	if len(sus) != 2 || sus[0] != "S_parent.GT" || sus[1] != "S_bulk.GT" {
		t.Fatalf("2p2b sus = %v, want [S_parent.GT S_bulk.GT]", sus)
	}

	res, sus = geneSpaceGTLines(cfg, "2b")
	if len(res) != 1 || res[0] != "R_bulk.GT" {
		t.Fatalf("2b res = %v, want [R_bulk.GT]", res)
	}
	if len(sus) != 1 || sus[0] != "S_bulk.GT" {
		t.Fatalf("2b sus = %v, want [S_bulk.GT]", sus)
	}
}
