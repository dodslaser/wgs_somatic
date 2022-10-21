#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 40

IGVDIR="$1"
SAMPLENAME="$2"


find $IGVDIR -name "*$SAMPLENAME*" -name "*_REALIGNED.bam" -exec petasuite -c -X --numthreads 40 -m bqfilt --validate full {} \;
#find $IGVDIR -name "*$SAMPLENAME*" -name "*_REALIGNED.bam" -exec echo {} \;
