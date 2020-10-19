# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
import shutil

rule validation_wf:
    input:
#        expand("{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_summary.txt", workingdir=workingdir, stype="tumor", tnsetting=tnscopesetting_list, sname=tumorid),
        expand("{workingdir}/{stype}/shared_to_igv.txt", workingdir=workingdir, stype="tumor")
    output:
        "{workingdir}/reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")

rule share_eval_to_igv:
    input:
        expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype="tumor"),
        expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype="normal"),
        expand("{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_{tnsetting}_FP.vcf", workingdir=workingdir, sname=tumorid, stype="tumor", tnsetting=tnscopesetting_list),
        expand("{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_{tnsetting}_FN.vcf", workingdir=workingdir, sname=tumorid, stype="tumor", tnsetting=tnscopesetting_list),
        expand("{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_{tnsetting}_somatic_w_normal_w_filters.vcf", workingdir=workingdir,  sname=tumorid, stype="tumor", tnsetting=tnscopesetting_list),
        expand("{workingdir}/{stype}/rtgeval/{tnsetting}_freqrange/{sname}_{tnsetting}_FP_freqrange.vcf", workingdir=workingdir, sname=tumorid, stype="tumor", tnsetting=tnscopesetting_list),
        expand("{workingdir}/{stype}/rtgeval/{tnsetting}_freqrange/{sname}_{tnsetting}_FN_freqrange.vcf", workingdir=workingdir, sname=tumorid, stype="tumor", tnsetting=tnscopesetting_list)
    params:
        igvdatadir = pipeconfig["rules"]["share_to_igv"]["igvdatadir"]
    output:
        "{workingdir}/{stype}/shared_to_igv.txt"
    run:
        igvsharedir = f"{params.igvdatadir}/mathias.johansson/"
        for sharefile in input:
            link_sharefile = os.path.abspath(sharefile)
            shell("ln -sf {link_sharefile} {igvsharedir}")
            if sharefile.endswith("REALIGNED.bam"):
                shell("ln -sf {link_sharefile}.bai {igvsharedir}")
        shell("echo {input} >> {output}")


rule filter_truthset:
    input: 
        expand("{workingdir}/{stype}/tnscope_given/{sname}_TNscope_tn_given.vcf", workingdir=workingdir, stype="tumor", sname=tumorid)
    params:
        rtg = config["rtg"]["tools"],
        bcftools = pipeconfig["rules"]["tnscope_vcffilter"]["bcftools"],
        high_thresh = config["data"]["thresholds"]["high"],
        low_thresh = config["data"]["thresholds"]["low"]
    output:
        "{workingdir}/{stype}/tnscope_given/truthset.filtered.ADabove2.vcf.gz",
        "{workingdir}/{stype}/tnscope_given/truthset.filtered.ADabove2.freqthreshold.vcf.gz"
    run:
        shell("{params.bcftools} filter -i 'FORMAT/AD[0:1] > 2' {input} > {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.vcf")
        shell("{params.bcftools} filter -i '(FORMAT/AD[0:1]*100)/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= {params.low_thresh} & (FORMAT/AD[0:1]*100)/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) <= {params.high_thresh}' {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.vcf > {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.freqthreshold.vcf")
        shell("{params.rtg} bgzip {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.vcf")
        shell("{params.rtg} index {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.vcf.gz")
        shell("{params.rtg} bgzip {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.freqthreshold.vcf")
        shell("{params.rtg} index {wildcards.workingdir}/{wildcards.stype}/tnscope_given/truthset.filtered.ADabove2.freqthreshold.vcf.gz")

rule rtgtools_eval:
    input:
        filtered_truthset = "{workingdir}/{stype}/tnscope_given/truthset.filtered.ADabove2.vcf.gz",
        calls = "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_somatic_{tnsetting}.vcf"
    params:
        rtg = config["rtg"]["tools"],
        sdf = config["rtg"]["sdf"],
        truthset = config["data"]["tset"],
        bedfile = config["data"]["bed"],
        tumorname = config["data"]["tumorname"]
    output:
        "{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_summary.txt",
        "{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_{tnsetting}_FP.vcf",
        "{workingdir}/{stype}/rtgeval/{tnsetting}/{sname}_{tnsetting}_FN.vcf"
    run:
        vcfbase = os.path.basename(f"{input.calls}")
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}"):
            shutil.copyfile(f"{input.calls}", f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}")
            shell(f"{params.rtg} bgzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}")
            shell("{params.rtg} index {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}.gz")
        if os.path.isdir(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}"):
            shutil.rmtree(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}")            
        shell("{params.rtg} vcfeval --squash-ploidy --sample={params.tumorname} -t {params.sdf} -e {params.bedfile} -b {input.filtered_truthset} -c {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}.gz -o {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/") 
        shell("mv {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/summary.txt {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_summary.txt")
        shell("cp {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/fp.vcf.gz {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_{wildcards.tnsetting}_FP.vcf.gz")
        shell("cp {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/fn.vcf.gz {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_{wildcards.tnsetting}_FN.vcf.gz")
        shell("gunzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_{wildcards.tnsetting}_FN.vcf.gz")
        shell("gunzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}/{wildcards.sname}_{wildcards.tnsetting}_FP.vcf.gz")


rule rtgtools_eval_freq:
    input:
        filtered_truthset = "{workingdir}/{stype}/tnscope_given/truthset.filtered.ADabove2.freqthreshold.vcf.gz",
        calls = "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_{tnsetting}_somatic_freqrange.vcf"
    params:
        rtg = config["rtg"]["tools"],
        sdf = config["rtg"]["sdf"],
        truthset = config["data"]["tset"],
        bedfile = config["data"]["bed"],
        tumorname = config["data"]["tumorname"]
    output:
        "{workingdir}/{stype}/rtgeval/{tnsetting}_freqrange/{sname}_summary.txt",
        "{workingdir}/{stype}/rtgeval/{tnsetting}_freqrange/{sname}_{tnsetting}_FP_freqrange.vcf",
        "{workingdir}/{stype}/rtgeval/{tnsetting}_freqrange/{sname}_{tnsetting}_FN_freqrange.vcf"
    run:
        vcfbase = os.path.basename(f"{input.calls}")
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}"):
            shutil.copyfile(f"{input.calls}", f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}")
            shell("{params.rtg} bgzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}")
            shell("{params.rtg} index {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}.gz")
        if os.path.isdir(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange"):
            shutil.rmtree(f"{wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange")
        shell("{params.rtg} vcfeval --squash-ploidy --sample={params.tumorname} -t {params.sdf} -e {params.bedfile} -b {input.filtered_truthset} -c {wildcards.workingdir}/{wildcards.stype}/rtgeval/{vcfbase}.gz -o {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/")
        shell("mv {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/summary.txt {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/{wildcards.sname}_summary.txt")
        shell("cp {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/fp.vcf.gz {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/{wildcards.sname}_{wildcards.tnsetting}_FP_freqrange.vcf.gz")
        shell("cp {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/fn.vcf.gz {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/{wildcards.sname}_{wildcards.tnsetting}_FN_freqrange.vcf.gz")
        shell("gunzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/{wildcards.sname}_{wildcards.tnsetting}_FN_freqrange.vcf.gz")
        shell("gunzip {wildcards.workingdir}/{wildcards.stype}/rtgeval/{wildcards.tnsetting}_freqrange/{wildcards.sname}_{wildcards.tnsetting}_FP_freqrange.vcf.gz")


