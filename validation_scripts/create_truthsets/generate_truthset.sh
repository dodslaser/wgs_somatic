#!/bin/bash -l

bedfile=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/intersecting_bedfile/na12878_na24385_intersect_confbed.bed
rtg=/apps/bio/dependencies/wgs_somatic/rtg-tools-3.11/rtg
#rtg=/apps/bio/software/rtg-tools/rtg-tools/dist/rtg-tools-3.9.1-6dde278/rtg
#hg38sdf=/apps/bio/dependencies/wgs_somatic/hg38/rtgeval/rtg/GRCh38.sdf
hg38sdf=/apps/bio/dependencies/wgs_somatic/hg38/rtgeval/rtg/GRCh38_hs38d1.sdf

na24385_truth=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA24385/HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf.gz
na12878_truth=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
#na24385_truth=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA24385/HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf
#na12878_truth=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf

output=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA12878_NA24385_truthset
$rtg vcfeval --squash-ploidy -e $bedfile -t $hg38sdf -c $na24385_truth -b $na12878_truth -o $output/eval

output=/medstore/Development/WGS_Somatic/artificial_coriell/truthset/hg38/NA24385_NA12878_truthset
$rtg vcfeval --squash-ploidy -e $bedfile -t $hg38sdf -b $na24385_truth -c $na12878_truth -o $output/eval
