# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile


if tumorid:
    rule share_to_resultdir:
        input:
            expand("{workingdir}/{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
            expand("{workingdir}/{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
            expand("{workingdir}/qc_report/{tumorname}_qc_stats.xlsx", workingdir=workingdir, tumorname=tumorname),
            expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference)            
        output:
            "{workingdir}/reporting/shared_result_files.txt"
        run:
            for resultfile in input:
                filebase = os.path.basename(f"{resultfile}")
                copyfile(f"{resultfile}", f"{wildcards.workingdir}/{filebase}")
            shell("echo {input} >> {output}")



else:
    rule share_to_resultdir:
        input:
            expand("{workingdir}/{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
            expand("{workingdir}/qc_report/{normalname}_qc_stats.xlsx", workingdir=workingdir, normalname=normalname),
            expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference)
        output:
            "{workingdir}/reporting/shared_result_files.txt"
        run:
            for resultfile in input:
                filebase = os.path.basename(f"{resultfile}")
                copyfile(f"{resultfile}", f"{wildcards.workingdir}/{filebase}")
            shell("echo {input} >> {output}")
