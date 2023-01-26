# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from shutil import copyfile


if tumorid:
    if normalid:
        rule share_to_resultdir:
            input:
                expand("{workingdir}/{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("{workingdir}/{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], caller="dnascope", sname=normalid, vcftype="germline"),
                expand("{workingdir}/qc_report/{tumorname}_qc_stats.xlsx", workingdir=workingdir, tumorname=tumorname),
                expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference),            
                expand("{workingdir}/{stype}/pindel/{sname}_pindel.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_somatic_CNV_observed.seg", workingdir=workingdir, vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_somatic_CNV_called.seg", workingdir=workingdir, vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam.bai", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam.bai", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir,  sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gz", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gzi.csi", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gz", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gzi.csi", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf.gzi.csi", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.gz.csi", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
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
                expand("{workingdir}/{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], caller="tnscope", sname=tumorid, vcftype="somatic"),
                expand("{workingdir}/qc_report/{tumorname}_qc_stats.xlsx", workingdir=workingdir, tumorname=tumorname),
                expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/pindel/{sname}_pindel.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam.bai", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gz", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaBNDs.vcf.gz.csi", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gz", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
                expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf.gzi.csi", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"])
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
            expand("{workingdir}/{stype}/dnascope/{sname}_{hgX}_SNV_CNV_germline.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"], hgX=reference),
            expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam.bai", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir,  sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf.gz.csi", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.gz", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf.csi", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
        output:
            "{workingdir}/reporting/shared_result_files.txt"
        run:
            for resultfile in input:
                filebase = os.path.basename(f"{resultfile}")
                copyfile(f"{resultfile}", f"{wildcards.workingdir}/{filebase}")
            shell("echo {input} >> {output}")
