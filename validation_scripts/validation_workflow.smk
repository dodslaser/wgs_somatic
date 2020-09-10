# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import shutil

rule validation_wf:
    input:
        expand("{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_summary.txt", workingdir=workingdir, stype="tumor", tnsetting=tnscopesetting_list, sname=tumorid),
        expand("{workingdir}/{stype}/tnscope_given/{sname}_TNscope_tn_given.vcf", workingdir=workingdir, stype="tumor", sname=tumorid)
    output:
        "{workingdir}/reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")

rule rtgtools_eval:
    input:
        "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_somatic_{tnsetting}.vcf"
        #"{workingdir}/tumor/tnscope/{tnsetting}/{sname}_somatic_{tnsetting}.vcf"
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
        shutil.copyfile(f"{input}", f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}")
        vcfbase = os.path.basename(f"{input}")
        #if not os.path.isfile(f"{wildcards.workingdir}/rtgeval/{vcfbase}.gz"):
        shell(f"{params.rtg} bgzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}")
        if os.path.isdir(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}"):
            shutil.rmtree(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}")            
        shell("{params.rtg} index {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}.gz")
        shell("{params.rtg} vcfeval --squash-ploidy -t {params.sdf} -e {params.bedfile} -b {params.truthset} -c {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}.gz -o {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/") 
        shell("mv {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/summary.txt {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_summary.txt")
