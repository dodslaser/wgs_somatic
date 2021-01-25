#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -l excl=1

export LD_PRELOAD=/usr/lib/petalink.so
export PETASUITE_REFPATH=/seqstore/software/petagene/corpus:/opt/petagene/petasuite/species

find args.outputdir -name "*.bam” -exec petasuite -c -X --numthreads 40 -m bqfilt --validate full {} \;
