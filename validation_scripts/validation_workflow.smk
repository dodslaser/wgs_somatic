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
        expand("{workingdir}/{stype}/tnscope/{sname}_somatic.vcf", workingdir=workingdir, stype="tumor", sname=tumorid)
    params:
        rtg = config["rtg"]["tools"],
        sdf = config["rtg"]["sdf"],
        truthset = config["data"]["bed"],
        bedfile = config["data"]["tset"]
    output:
        "{workingdir}/rtgeval/summary.txt"
    run:
        shell("{params.rtg} bgzip {input}")
        shell("{params.rtg} index {input}.gz")
        shell("{params.rtg} vcfeval --squash-ploidy -e {params.bed} -b {params.truthset} -c {input}.gz -o {wildcards.workingdir}/rtgeval/") 
