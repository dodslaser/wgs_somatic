#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1
#$ -pe mpi 40

export LD_PRELOAD=/usr/lib/petalink.so
export PETASUITE_REFPATH=/seqstore/software/petagene/corpus:/opt/petagene/petasuite/species


OUTPUTDIR=$(echo "$1")

find $OUTPUTDIR -name "*_REALIGNED.bam" -exec petasuite -c -X --numthreads 40 -m bqfilt --validate full {} \;
find $OUTPUTDIR -name "*_realignedTNscope.bam" -exec petasuite -c -X --numthreads 40 -m bqfilt --validate full {} \;

find $OUTPUTDIR -name "*.bam" -not -name "*_REALIGNED.bam" -not -name "*_realignedTNscope.bam" -exec rm {} \;
find $OUTPUTDIR -name "*.bam.bai" -not -name "*_REALIGNED.bam.bai" -not -name "*_realignedTNscope.bam.bai" -exec rm {} \;
