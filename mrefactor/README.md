# GoBSAseq (mrefactor)

Refactored BSA-seq pipeline with the same CLI shape as the parent [GoBSAseq](../README.md) repo, incorporating literature methods and lessons from [genome-whisperer/bsaseq](https://github.com/gmaffy/genome-whisperer/tree/main/bsaseq).

See [IMPLEMENTATION.md](./IMPLEMENTATION.md) for design notes and comparison tables.

## Quick start

```powershell
cd mrefactor
go build -o GoBSAseq.exe .
.\GoBSAseq.exe --variant your.vcf.gz --parents Res,Sus --bulks High,Low --out results
```

Results are written to `results/GoBSAseqResults/<timestamp>/`.
