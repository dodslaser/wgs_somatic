# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from workflows.scripts.gender import calc_gender
from workflows.scripts.create_segfile import create_seg
from workflows.scripts.fix_sexploidyfile import mod_sex_vcf

if tumorid:
    rule filter_canvas_somatic:
        input:
            expand("{workingdir}/{stype}/canvas/{sname}_somatic_CNV.vcf.gz", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"])
        params:
            annotate = pipeconfig["rules"]["canvas"]["annotate"],
            annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"]
        output:
            "{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx",
            "{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf"
        run:
            shell("gunzip {input}")
            shell("grep -v 'Canvas:REF' {wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_somatic_CNV.vcf > {wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf")
            shell("{params.annotate} -v {wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf -g {params.annotate_ref} -o {wildcards.workingdir}/{wildcards.stype}/canvas/")
            os.rename(f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf.xlsx", f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf.xlsx")
            os.rename(f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic_noref.vcf", f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf")

rule filter_canvas_germline:
    input:
        "{workingdir}/{stype}/canvas/{sname}_germline_CNV.vcf.gz"
    params:
        annotate = pipeconfig["rules"]["canvas"]["annotate"],
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"]
    output:
        "{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf.xlsx",
        "{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf"
    run:
        shell("gunzip {input}")
        shell("grep -v 'Canvas:REF' {wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_germline_CNV.vcf > {wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf")
        shell("{params.annotate} -v {wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf -g {params.annotate_ref} -o {wildcards.workingdir}/{wildcards.stype}/canvas/")
        os.rename(f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf.xlsx", f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf.xlsx")
        os.rename(f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline_noref.vcf", f"{wildcards.workingdir}/{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf")


if tumorid:
    rule canvas_somatic:
        input:
            germline_snv_vcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            somatic_vcf = expand("{workingdir}/{stype}/tnscope/{sname}_somatic.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            bamfile = "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam",
            normal_wgscov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            normal_ycov = expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
        params:
            genomeversion = config["reference"],
            dll = pipeconfig["singularities"]["canvas"]["dll"],
            annotate = pipeconfig["rules"]["canvas"]["annotate"],
            annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
            genomedir = pipeconfig["singularities"]["canvas"]["reference"],
            kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
            run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
            filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
            samplename = sampleconfig["tumorname"]
        singularity:
            pipeconfig["singularities"]["canvas"]["sing"]
        output:
            "{workingdir}/{stype}/canvas/{sname}_somatic_CNV.vcf.gz",
            "{workingdir}/{stype}/canvas/{sname}_somatic_CNV_observed.seg",
            "{workingdir}/{stype}/canvas/{sname}_somatic_CNV_called.seg"
        shell:
            "echo $HOSTNAME;"
            "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bamfile} --normal_vcf {input.germline_snv_vcf} --o {wildcards.workingdir}/{wildcards.stype}/canvas/ -t TN --samplename {wildcards.sname} --wgscovfile {input.normal_wgscov} --ycovfile {input.normal_ycov} --somatic_vcf {input.somatic_vcf} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"
 
rule canvas_germline:
    input:
        germline_snv_vcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        bamfile = "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam",
        normal_wgscov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        normal_ycov = expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        genomeversion = config["reference"],
        dll = pipeconfig["singularities"]["canvas"]["dll"],
        annotate = pipeconfig["rules"]["canvas"]["annotate"],
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
        genomedir = pipeconfig["singularities"]["canvas"]["reference"],
        kmerfile = pipeconfig["singularities"]["canvas"]["kmerfile"],
        run_py = pipeconfig["singularities"]["canvas"]["tool_path"],
        filter13 = pipeconfig["singularities"]["canvas"]["filter13"],
        samplename = sampleconfig["normalname"]
    singularity:
        pipeconfig["singularities"]["canvas"]["sing"]
    output:
        "{workingdir}/{stype}/canvas/{sname}_germline_CNV.vcf.gz",
        "{workingdir}/{stype}/canvas/{sname}_germline_CNV_observed.seg",
        "{workingdir}/{stype}/canvas/{sname}_germline_CNV_called.seg"
    shell:
        "echo $HOSTNAME;"
        "{params.run_py} --genomeversion {params.genomeversion} --bam {input.bamfile} --normal_vcf {input.germline_snv_vcf} --o {wildcards.workingdir}/{wildcards.stype}/canvas/ -t germline --samplename {wildcards.sname} --wgscovfile {input.normal_wgscov} --ycovfile {input.normal_ycov} --referencedir {params.genomedir} --kmerfile {params.kmerfile} --canvasdll {params.dll} --filterfile {params.filter13}"


rule convert_to_alissaformat:
    input:
        germline_cnv_vcf = expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        converter = pipeconfig["rules"]["convert_to_alissaformat"]["converter"],
        referencegenome = pipeconfig["referencegenome"]
    output:
        "{workingdir}/{stype}/canvas/{sname}_CNV_germline_alissaformat.vcf"
    run:
        shell(f"{params.python} {params.converter} -l -q {input} {output} {referencegenome}")
