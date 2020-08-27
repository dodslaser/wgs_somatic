#!/bin/bash -l

snv=$(echo "$1")
cnv=$(echo "$2")
outputdir=$(echo "$3")
if [ -z "$snv" ] || [ -z "$cnv" ] || [ -z "$outputdir" ]; then
    echo "arguments missing, supply -- snv cnv outputdir "
    echo "snv $snv"
    echo "cnv $cnv"
    echo "outputdir $outputdir"    
    exit
fi
outputvcf=$(basename "$1" | sed 's/.vcf$//g')_SNV_CNV.vcf.gz
alissacnv=$(echo "$cnv" | sed 's/.vcf$//g')_alissaformatted.vcf
if [ ! -f "$alissacnv" ]; then
    ./canvas_to_interpreter/canvasvcf_to_interpreter.py -l -q $cnv $alissacnv
fi
module load bcftools
if [ ! -f "$outputdir/$outputvcf" ]; then
    bcftools view $snv -Ob -o ${snv}.bgz
    bcftools view $snv -Ob -o ${alissacnv}.bgz
    bcftools index ${snv}.bgz
    bcftools index ${alissacnv}.bgz
    bcftools concat --allow-overlaps ${snv}.bgz ${alissacnv}.bgz -Oz -o $outputdir/$outputvcf
    rm -f ${snv}.bgz*
    rm -f ${alissacnv}.bgz* 
fi
