#!/bin/bash -l

confbed_na12878=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed
confbed_na24385=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA24385/HG002_GRCh38_1_22_v4.1_draft_benchmark.bed
outdir=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/intersecting_bedfile

module load bedtools
bedtools intersect -a $confbed_na12878 -b $confbed_na24385 >> $outdir/na12878_na24385_intersect_confbed.bed
