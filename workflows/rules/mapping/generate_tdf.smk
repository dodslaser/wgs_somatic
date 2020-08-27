# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule generate_tdf:
    input:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    params:
        igvjar = pipeconfig["rules"]["generate_tdf"]["tdfgen"]
    output:
        "{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf"
    run:
        shell("nohup java -Xmx32768m -jar {params.igvjar} count {input} {output} hg19")
