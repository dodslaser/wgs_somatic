#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 40

SAMPLE=$(echo "$1")
VCF=$(echo "$2")
SIZE=$(echo "$3")
REF=$(echo "$4")

conda deactivate
module load anaconda2
source activate wopr_alissa

python workflows/scripts/Alissa_API_tools/alissa_API_tools.py -a $SAMPLE -v $VCF -s $SIZE -n $SAMPLE -ref $REF

