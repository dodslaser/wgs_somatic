#!/bin/bash -l

outputdir=$(echo "$1")
dag=$(echo "$2")

date=$(date '+%Y-%m-%d')
if [ ! -d "$outputdir" ]; then
    mkdir $outputdir
fi
if [ ! -d "$outputdir" ]; then
    echo "$outputdir does not exist and could not be created, exiting"
    exit
fi
logdir=$(echo "$outputdir/qsublogs")
if [ ! -d "$logdir" ]; then
    mkdir $logdir
fi

configs=$(echo "$outputdir/configs")
if [ ! -d "$configs" ]; then
    mkdir $configs
fi
for config in $(ls configs); do
    cp -f configs/$config $configs
done

module load anaconda2
source activate wgs_somatic

if [ -z "$dag" ];
then
    snakemake -s pipeline.snakefile --cluster-config configs/cluster.yaml --cluster "qsub -S /bin/bash -pe mpi {cluster.threads} -q {cluster.queue} -N {cluster.name} -o $logdir/{cluster.output} -e $logdir/{cluster.error} -l {cluster.excl}" --jobs 999 --latency-wait 60 --directory $outputdir &>> $logdir/log_$date
else
    snakemake -s pipeline.snakefile --dag | dot -Tsvg > /seqstore/webfolders/wgs/admin/barncancer/dags/dag_${date}.svg
fi
