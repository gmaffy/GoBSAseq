# GoBSAseq (crefactor)

Refactored BSA-seq pipeline with the same CLI shape as the parent [GoBSAseq](../README.md) repo, incorporating literature methods and lessons from [genome-whisperer/bsaseq](https://github.com/gmaffy/genome-whisperer/tree/main/bsaseq).

The refactor keeps direct Go VCF processing and the existing output workflow, but expands the analysis core around QTL-seq/MutMap SNP-index logic, Magwene G/ED statistics, PyBSASeq-style Fisher windows, BRM-style interval support, robust Z-score overlays, Wilson confidence intervals, and FDR utilities.

See [IMPLEMENTATION.md](./IMPLEMENTATION.md) for design notes and comparison tables.

## Quick start

```powershell
cd crefactor
go build -o GoBSAseq.exe .
.\GoBSAseq.exe --variant your.vcf.gz --parents Res,Sus --bulks High,Low --out results
```

Results are written to `results/GoBSAseqResults/<timestamp>/`.
