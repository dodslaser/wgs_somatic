# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import shutil

rule validation_wf:
    input:
        expand("{workingdir}/rtgeval/{tnsetting}/{sname}_summary.txt", workingdir=workingdir, tnsetting=tnscopesetting_list, sname=tumorid),
        expand("{workingdir}/{stype}/tnscope_given/{sname}_TNscope_tn_given.vcf", workingdir=workingdir, stype="tumor", sname=tumorid)
    output:
        "{workingdir}/reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")

rule rtgtools_eval:
    input:
        #expand("{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_somatic_{tnsetting}.vcf", tnsetting=tnscopesetting_list, workingdir=workingdir, stype="tumor", sname=tumorid)
        "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_somatic_{tnsetting}.vcf"
    params:
        rtg = config["rtg"]["tools"],
        sdf = config["rtg"]["sdf"],
        truthset = config["data"]["tset"],
        bedfile = config["data"]["bed"]
    output:
        "{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_summary.txt"
    run:
        #if not os.path.isdir("{wildcards.workingdir}/rtgeval/"):
        #    os.mkdir(f"{wildcards.workingdir}/rtgeval/")
        vcfbase = os.path.basename(f"{input}")
        shutil.copyfile(f"{input}", f"{wildcards.workingdir}/rtgeval/{vcfbase}")
        vcfbase = os.path.basename(f"{input}")
        #if not os.path.isfile(f"{wildcards.workingdir}/rtgeval/{vcfbase}.gz"):
        shell(f"{params.rtg} bgzip {wildcards.workingdir}/rtgeval/{vcfbase}")
        if os.path.isdir(f"{wildcards.workingdir}/rtgeval/{wildcards.tnsetting}"):
            shutil.rmtree(f"{wildcards.workingdir}/rtgeval/{wildcards.tnsetting}")            
        shell("{params.rtg} index {wildcards.workingdir}/rtgeval/{vcfbase}.gz")
        shell("{params.rtg} vcfeval --squash-ploidy -t {params.sdf} -e {params.bedfile} -b {params.truthset} -c {wildcards.workingdir}/rtgeval/{vcfbase}.gz -o {wildcards.workingdir}/rtgeval/{wildcards.tnsetting}/") 
        shell("mv {wildcards.workingdir}/rtgeval/{wildcards.tnsetting}/eval/summary.txt {wildcards.workingdir}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_summary.txt")
