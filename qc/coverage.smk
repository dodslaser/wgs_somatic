# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule wgs_coverage:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["wgs_coverage"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome
    output:
        "{stype}/reports/{sname}_WGScov.tsv"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -i {input} -r {params.referencegenome} --interval 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT --algo WgsMetricsAlgo {output}")

rule y_coverage:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["y_coverage"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome
    output:
        "{stype}/reports/{sname}_Ycov.tsv"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -i {input} -r {params.referencegenome} --interval Y --algo WgsMetricsAlgo {output}")
