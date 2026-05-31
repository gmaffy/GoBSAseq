# GoBSAseq refactor (`mrefactor/`)

## Research synthesis (one-bulk vs two-bulk)

| Design | Literature terms | Null at unlinked loci | Primary statistics |
|--------|------------------|----------------------|-------------------|
| **Two-bulk** | Bidirectional QTL-seq, classical BSA-seq | Equal allele frequency in both pools | Δ(SNP index), G′ (Magwene), ED⁴, LOD, beta-binomial BF, Fisher 2×2 |
| **One-bulk** | Unidirectional QTL-seq, MutMap-style | Mendelian expectation in single pool (0.5 F₂; 0.75/0.25 BC) | SNP index, \|SI−p₀\|, one-bulk G/LOD/BF |

Reference pipelines: QTL-seq, PyBSASeq (Fisher + significant-SNP windows), BRM (block regression + CIs), MULTIPOOL, QTLseqr.

## Comparison: `genome-whisperer/bsaseq` vs root `GoBSAseq` vs `mrefactor`

| Aspect | genome-whisperer `bsaseq` | Root `GoBSAseq` | `mrefactor` |
|--------|---------------------------|-----------------|-------------|
| Input | Filtered TSV from external VCF step | Streaming VCF (`vcfgo`) | Streaming VCF |
| Statistics | SNP index, Δ-index, G only | + ED⁴, LOD, BF, BRM, robust Z | Shared `stats/` package + Fisher + SigSNP window ratio |
| One-bulk high pool | Supported in filters | Not wired in `run` | `RunTwoParentsHighBulk` |
| Output layout | Timestamped `goBSAseqResults/` | Flat `--out` | Timestamped `GoBSAseqResults/` |
| Code shape | Monolithic package | `cmd` / `run` / `utils` / `twobulk` / `onebulk` | Same layout + `stats/` |

## Architecture

```
mrefactor/
  main.go
  cmd/root.go          # Cobra CLI (same flags as root + fisher-alpha)
  run/run.go           # Mode dispatch + timestamped output dir
  utils/               # VCF hard filter, config, SimulateAF, PrepareRunDir
  stats/               # Pure statistics (literature-aligned)
  twobulk/pipeline.go  # Two-bulk VCF → stats → smooth → QTL → HTML
  onebulk/pipeline.go  # One-bulk (low or high pool via bsaFilterType)
```

## Analysis modes

1. **Two parents + two bulks** — `twobulk.RunTwoBulkTwoParents`
2. **Bulks only** — `twobulk.RunTwoBulksOnly`
3. **Two parents + low bulk** — `onebulk.RunTwoParentsLowBulk` (filter type 2)
4. **Two parents + high bulk** — `onebulk.RunTwoParentsHighBulk` (filter type 3)

## Improvements implemented

- **`stats` package**: central G-test, ED⁴, LOD, beta-binomial BF, Fisher exact, population-aware null AF.
- **PyBSASeq-style `SigSNPRatio`**: per sliding window, fraction of SNPs with Fisher p ≤ `--fisher-alpha`.
- **Fisher p-value** on raw marker TSV (`FisherP` column).
- **BC recurrent flag** (`--recurrent`) for 0.75 vs 0.25 null in backcross one-bulk stats.
- **Timestamped run directories** under `GoBSAseqResults/`.
- **Complete mode routing** in `run.resolveAnalysisMode` (CLI and interactive).

## Build

```powershell
cd mrefactor
go mod tidy
go build -o GoBSAseq.exe .
```

## Future work

- LOESS block regression CIs exactly as Huang et al. BRM (currently BRM-style blocks on smoothed curves).
- Optional TSV-only path compatible with genome-whisperer preprocessing.
- Unit tests for `stats` hypergeometric/Fisher and threshold simulation.
