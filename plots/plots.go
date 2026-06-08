package plots

import (
	"sort"

	"github.com/gmaffy/GoBSAseq/stats"
	"github.com/gmaffy/GoBSAseq/utils"
	"github.com/go-echarts/go-echarts/v2/components"
)

func Plot(allSmoothed []stats.SmoothedStats, cfg utils.AnalysisConfig) error {
	if len(allSmoothed) == 0 {
		return nil
	}

	byChr := make(map[string][]stats.SmoothedStats)
	chroms := make([]string, 0)

	for _, s := range allSmoothed {
		if _, ok := byChr[s.CHROM]; !ok {
			chroms = append(chroms, s.CHROM)
		}
		byChr[s.CHROM] = append(byChr[s.CHROM], s)
	}

	sort.Strings(chroms)

	individualPage := components.NewPage()
	individualPage.SetLayout(components.PageFlexLayout)
	individualPage.PageTitle = "GoBSAseq — Individual Smoothed Statistics"

	normalisedPage := components.NewPage()
	normalisedPage.SetLayout(components.PageFlexLayout)
	normalisedPage.PageTitle = "GoBSAseq — Normalised Smoothed Stats"

	compositePage := components.NewPage()
	compositePage.SetLayout(components.PageFlexLayout)
	compositePage.PageTitle = "GoBSAseq — GoStat (Composite statistics)"

	for _, chrom := range chroms {
		stat := byChr[chrom]
		n := len(stat)
		x := make([]int64, n)
		hi := make([]float64, n)
		li := make([]float64, n)
		dsi := make([]float64, n)
		gs := make([]float64, n)
		ed := make([]float64, n)
		lod := make([]float64, n)
		bbl := make([]float64, n)
		hiT99, hiTM99 := make([]float64, n), make([]float64, n)
		hiT95, hiTM95 := make([]float64, n), make([]float64, n)
		liT99, liTM99 := make([]float64, n), make([]float64, n)
		liT95, liTM95 := make([]float64, n), make([]float64, n)
		dsiT99, dsiTM99 := make([]float64, n), make([]float64, n)
		dsiT95, dsiTM95 := make([]float64, n), make([]float64, n)
		gsT99, gsT95 := make([]float64, n), make([]float64, n)
		edT99, edT95 := make([]float64, n), make([]float64, n)
		lodT99, lodT95 := make([]float64, n), make([]float64, n)
		bblT99, bblT95 := make([]float64, n), make([]float64, n)
		if len(stat) == 0 {
			continue
		}
		for i, s := range stat {
			x[i] = s.POS
			hi[i], li[i], dsi[i] = s.SmHighSI, s.SmLowSI, s.SmDeltaSI
			gs[i], ed[i], lod[i], bbl[i] = s.SmGstat, s.SmED, s.SmLOD, s.SmBBLogBF
			t := s.Thresholds
			hiT99[i], hiTM99[i] = t.HighP99, t.HighMp99
			hiT95[i], hiTM95[i] = t.HighP95, t.HighMp95
			liT99[i], liTM99[i] = t.LowP99, t.LowMp99
			liT95[i], liTM95[i] = t.LowP95, t.LowMp95
			dsiT99[i], dsiTM99[i] = t.DsiP99, t.DsiMp99
			dsiT95[i], dsiTM95[i] = t.DsiP95, t.DsiMp95
			gsT99[i], gsT95[i] = t.GsP99, t.GsP95
			edT99[i], edT95[i] = t.EdP99, t.EdP95
			lodT99[i], lodT95[i] = t.LodP99, t.LodP95
			bblT99[i], bblT95[i] = t.BbP99, t.BbP95
		}
	}

	return nil
}
