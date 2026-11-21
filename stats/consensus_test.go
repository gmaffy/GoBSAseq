package stats

import (
	"math"
	"testing"
)

// consensusInput: one chromosome, markers every 10 kb; HMM LOD above 5 at
// markers lo..hi, CompositeZ above 10 at czLo..czHi.
func consensusInput(hmmLo, hmmHi, czLo, czHi int) []SmoothedStats {
	sm := make([]SmoothedStats, 200)
	for i := range sm {
		sm[i] = SmoothedStats{CHROM: "c", POS: int64(10_000 * (i + 1))}
		if i >= hmmLo && i <= hmmHi {
			sm[i].HMMLOD = 8 + float64(i-hmmLo)
		}
		if i >= czLo && i <= czHi {
			sm[i].CompositeZ = 12
		}
	}
	return sm
}

func consensusThresholds() Thresholds {
	var th Thresholds
	th.HMM.LODP99, th.Z.CompositeZP99 = 5, 10
	return th
}

func TestConsensusAllAgree(t *testing.T) {
	sm := consensusInput(40, 60, 35, 55) // HMM 410–610 kb, CompositeZ 360–560 kb
	brm := []BRMBlock{{Chrom: "c", Start: 450_000, Stop: 700_000, PeakPos: 500_000, Peak: 0.4, Threshold: 0.2, PeakZ: 6}}
	ind := []PeakIntersection{{Stat: "Gstat", Chrom: "c", Start: 400_000, End: 500_000}, {Stat: "ED4", Chrom: "c", Start: 900_000, End: 950_000}}
	q := ConsensusQTLs(sm, consensusThresholds(), brm, ind, true, true, 0, 1)
	if len(q) != 1 {
		t.Fatalf("loci %+v", q)
	}
	f := q[0]
	// core: HMM ∩ CompositeZ ∩ BRM ≈ [450 kb, 560 kb]; span ≈ [360 kb, 700 kb]
	if !f.Consensus || f.Support != "HMMLOD+BRM+CompositeZ" || f.Direction != "+" || f.LeadMethod != "HMMLOD" || f.HMMSupportStart != 600_000 || f.HMMSupportStop != 610_000 || f.LeadPos != 605_000 ||
		f.Start != 450_000 || math.Abs(float64(f.Stop)-560_000) > 10_000 || f.SpanStop != 700_000 || f.SpanStart > 360_000 ||
		f.BRM.Z != 6 || f.Statistics != "Gstat" {
		t.Fatalf("locus %+v", f)
	}
}

func TestConsensusPartialAndNoHMM(t *testing.T) {
	sm := consensusInput(-1, -1, 35, 55) // no HMM region
	brm := []BRMBlock{{Chrom: "c", Start: 400_000, Stop: 520_000, PeakPos: 450_000, Peak: -0.4, PeakZ: -5}}
	q := ConsensusQTLs(sm, consensusThresholds(), brm, nil, true, true, 0, 1)
	if len(q) != 1 || q[0].Consensus || q[0].Support != "BRM+CompositeZ" || q[0].Direction != "−" || q[0].LeadMethod != "BRM" {
		t.Fatalf("HMM available but absent: %+v", q)
	}
	// 2b without --ref-parent: no HMM, polarity-free BRM: BRM + CompositeZ is a consensus, direction none
	q = ConsensusQTLs(sm, consensusThresholds(), brm, nil, false, false, 0, 1)
	if len(q) != 1 || !q[0].Consensus || q[0].Support != "BRM+CompositeZ" || q[0].Direction != "none" {
		t.Fatalf("no HMM: %+v", q)
	}
}

func TestConsensusSeparateLociAndMixed(t *testing.T) {
	sm := consensusInput(20, 30, 20, 30)
	for i := 120; i <= 130; i++ { // a second locus, CompositeZ only
		sm[i].CompositeZ = 15
	}
	brm := []BRMBlock{
		{Chrom: "c", Start: 200_000, Stop: 250_000, PeakPos: 220_000, Peak: 0.3},
		{Chrom: "c", Start: 260_000, Stop: 300_000, PeakPos: 280_000, Peak: -0.5},
	}
	q := ConsensusQTLs(sm, consensusThresholds(), brm, nil, true, true, 0, 1)
	if len(q) != 2 || !q[0].Consensus || q[0].Direction != "mixed" || q[0].BRM.Value != -0.5 || q[1].Consensus || q[1].Support != "CompositeZ" {
		t.Fatalf("loci %+v", q)
	}
}

// Regions joined only through the merge gap have no common core: the core
// is the lead (HMM) region, whose ends are interpolated threshold crossings.
func TestConsensusGapCore(t *testing.T) {
	sm := consensusInput(40, 50, 60, 70) // HMM 410–510 kb, CompositeZ 610–710 kb
	brm := []BRMBlock{{Chrom: "c", Start: 520_000, Stop: 600_000, PeakPos: 560_000, Peak: 0.3}}
	if q := ConsensusQTLs(sm, consensusThresholds(), brm, nil, true, true, 0, 1); len(q) != 3 {
		t.Fatalf("gap 0: want 3 separate loci, got %+v", q)
	}
	q := ConsensusQTLs(sm, consensusThresholds(), brm, nil, true, true, 50_000, 1)
	if len(q) != 1 || !q[0].Consensus || q[0].LeadMethod != "HMMLOD" || q[0].Start > 410_000 || q[0].Stop < 510_000 || q[0].Stop > 520_000 {
		t.Fatalf("gap 50 kb: %+v", q)
	}
}

// --min-qtl-markers applies to the HMM and CompositeZ regions here too.
func TestConsensusMinMarkers(t *testing.T) {
	sm := consensusInput(40, 42, 40, 42) // 3 markers each
	brm := []BRMBlock{{Chrom: "c", Start: 400_000, Stop: 440_000, PeakPos: 420_000, Peak: 0.3}}
	if q := ConsensusQTLs(sm, consensusThresholds(), brm, nil, true, true, 0, 5); len(q) != 1 || q[0].Support != "BRM" {
		t.Fatalf("min 5: %+v", q)
	}
}

// On a flat LOD plateau the support interval covers the plateau and the lead
// position is its midpoint, not the (arbitrary) maximum.
func TestLODSupportPlateau(t *testing.T) {
	sm := make([]SmoothedStats, 100)
	for i := range sm {
		sm[i] = SmoothedStats{CHROM: "c", POS: int64(1000 * (i + 1))}
		if i >= 20 && i < 80 {
			sm[i].HMMLOD = 11.5 + 0.3*math.Sin(float64(i)) // plateau with noise, all within 1.5 of the maximum
		}
	}
	sm[30].HMMLOD = 12.4 // arbitrary maximum near one end
	s, e := lodSupport(sm, 31_000, 1.5)
	if s != 21_000 || e != 80_000 {
		t.Fatalf("support %d–%d, want 21000–80000", s, e)
	}
}
