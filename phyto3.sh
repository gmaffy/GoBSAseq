#!/bin/bash

go run main.go \
	--parents-bams /mnt/v/DATA/moschata/2025/HM_PHYBCF1_BSASEQ/reference_genomes/v1/bams/HM.RGMD_bqsr.cram,/mnt/v/DATA/moschata/2026/MOB139/reference_genomes/v1/bams/MOB139.rgmd_bqsr.cram \
	--bulks-bams /mnt/v/DATA/moschata/2026/Phyto3_HR_BULK_BC2F2/reference_genomes/v1/bams/Phyto3_HR_BULK_BC2F2.sorted.RGMD_bqsr.cram,/mnt/v/DATA/moschata/2026/Phyto3_SUS_breakoff_BULK_BC2F2/reference_genomes/v1/bams/Phyto3_SUS_breakoff_BULK_BC2F2.sorted.RGMD_bqsr.cram \
	--snpEffDB Moschatav1 \
	-r /mnt/z/genomes/moschata/v1/assembly/Cmoschata_genome_v1.fa \
	--prg ~/tools/genome-whisperer/variant_annotation/PRG_FILES/MOSCHATA_PRG_TOP_HITS.txt \
	--gene-descriptions ~/tools/genome-whisperer/variant_annotation/GENE_DESCRIPTIONS/MOSCHATA_gene_descriptions.txt \
	--gff /mnt/z/genomes/moschata/v1/annotation/Cmoschata_v1.gff3 \
	--population BC2H \
	-o /mnt/f/2026_Phyto3 \
	--caller DeepVariant \
	--merger glnexus

