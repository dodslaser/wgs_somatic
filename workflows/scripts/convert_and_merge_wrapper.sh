#!/bin/bash -l

for dir in $(ls /seqstore/webfolders/wgs/barncancer); do
    maindir=/seqstore/webfolders/wgs/barncancer/$dir
    tumorsnv=$(find $maindir/tumor -name "*somatic_refseq3kfilt.vcf" | head -n1)
    tumorcnv=$(find $maindir/tumor -name "*CNV_somatic.vcf" | head -n1)
    normalsnv=$(find $maindir/normal -name "*germline_refseq3kfilt.vcf" | head -n1)
    normalcnv=$(find $maindir/normal -name "*CNV_germline.vcf" | head -n1)
    if [ -z "$(ls $maindir | grep SNV_CNV | grep somatic)" ];
    then
        ./convertCNV_and_mergeSNV.sh $tumorsnv $tumorcnv $maindir
    fi
    if [ -z "$(ls $maindir | grep SNV_CNV | grep germline)" ];
    then
    ./convertCNV_and_mergeSNV.sh $normalsnv $normalcnv $maindir
    fi
done
