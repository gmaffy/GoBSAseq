#!/bin/bash

go run main.go \
	--parents-bams /mnt/v/DATA/moschata/2025/HM_PHYBCF1_BSASEQ/reference_genomes/v1/bams/HM.RGMD_bqsr.cram,/mnt/v/DATA/moschata/2025/MO991_PM/reference_genomes/v1/bams/MO991_PM.RGMD_bqsr.cram \
	--bulks-bams /mnt/v/DATA/moschata/2026/Phyto2_RUS_BULK_F2/reference_genomes/v1/bams/Phyto2_RUS_BULK_F2.sorted.RGMD_bqsr.cram,/mnt/v/DATA/moschata/2026/Phyto2_SUS_BULK_F2/reference_genomes/v1/bams/Phyto2_SUS_BULK_F2.sorted.RGMD_bqsr.cram \
	-r /mnt/z/genomes/moschata/v1/assembly/Cmoschata_genome_v1.fa \
	--caller DeepVariant \
	--merger glnexus \
	--snpEffDB Moschatav1 \
	--gff /mnt/z/genomes/moschata/v1/annotation/Cmoschata_v1.gff3 \
	--prg /home/godwin/tools/genome-whisperer/variant_annotation/PRG_FILES/MOSCHATA_PRG_TOP_HITS.txt \
	--gene-descriptions /home/godwin/tools/genome-whisperer/variant_annotation/GENE_DESCRIPTIONS/MOSCHATA_gene_descriptions.txt \
	-o /mnt/e/PHY2_F2 \
	--verbose
