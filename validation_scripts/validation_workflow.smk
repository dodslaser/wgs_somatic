# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule validation_wf:
    input:
        "{workingdir}/rtgeval/summary.txt"
    output:
        "{workingdir}/reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")


rule rtgtools_eval:
    input:
        "{workingdir}/{stype}/tnscope/{sname}_somatic.vcf"
    params:
        rtg =
        sdf = 
        truthset = 
        bedfile =  
    output:
        "{workingdir}/rtgeval/summary.txt"
    run:
