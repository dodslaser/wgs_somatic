# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile

rule share_to_resultdir:
    input:
        expand("qc_report/{tumorname}_qc_stats.xlsx", tumorname=tumorname),
        expand("{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorid, stype=sampleconfig[tumorname]["stype"])
    output:
        "reporting/shared_result_files.txt"
    run:
        for resultfile in input:
            filebase = os.path.basename(f"{resultfile}")
            copyfile(f"{resultfile}", f"{filebase}")
        shell("echo {input} >> {output}")
