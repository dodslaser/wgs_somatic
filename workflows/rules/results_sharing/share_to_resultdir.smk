# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile

rule share_to_resultdir:
    input:
        expand("{workingdir}/qc_report/{tumorname}_qc_stats.xlsx", workingdir=workingdir, tumorname=tumorname),
        #expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        #expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        #expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"])
    output:
        "{workingdir}/reporting/shared_result_files.txt"
    run:
        for resultfile in input:
            filebase = os.path.basename(f"{resultfile}")
            copyfile(f"{resultfile}", f"{filebase}")
        shell("echo {input} >> {output}")
