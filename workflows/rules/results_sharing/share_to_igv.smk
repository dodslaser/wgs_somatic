# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os

if tumorid:
    rule share_to_igv:
        input:
            expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf", workingdir=workingdir, sname=normalid, vartype="germline", stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf", workingdir=workingdir, sname=tumorid, vartype="somatic", stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_somatic_CNV_observed.seg", workingdir=workingdir, vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_somatic_CNV_called.seg", workingdir=workingdir, vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir,  sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaBNDs.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/tnscope/{sname}_REALIGNED_realignedTNscope.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"])
        params:
            updateigv = pipeconfig["rules"]["share_to_igv"]["updateigv"],
            igvdatadir = pipeconfig["rules"]["share_to_igv"]["igvdatadir"]
        output:
            "{workingdir}/reporting/shared_igv_files.txt" 
        run:
            igvsharedir = f"{params.igvdatadir}/{igvuser}/"
            for sharefile in input:
                print(sharefile)
                link_sharefile = os.path.abspath(sharefile)
                shell("ln -sf {link_sharefile} {igvsharedir}")
                if sharefile.endswith("REALIGNED.bam"):
                    shell("ln -sf {link_sharefile}.bai {igvsharedir}")
            shell("{params.updateigv}")
            shell("echo {input} >> {output}")

else:

    rule share_to_igv:
        input:
            expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf", workingdir=workingdir, sname=normalid, vartype="germline", stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/canvas/{sname}_germline_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir,  sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
            expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
        params:
            updateigv = pipeconfig["rules"]["share_to_igv"]["updateigv"],
            igvdatadir = pipeconfig["rules"]["share_to_igv"]["igvdatadir"]
        output:
            "{workingdir}/reporting/shared_igv_files.txt"
        run:
            igvsharedir = f"{params.igvdatadir}/{igvuser}/"
            for sharefile in input:
                print(sharefile)
                link_sharefile = os.path.abspath(sharefile)
                shell("ln -sf {link_sharefile} {igvsharedir}")
                if sharefile.endswith("REALIGNED.bam"):
                    shell("ln -sf {link_sharefile}.bai {igvsharedir}")
            shell("{params.updateigv}")
            shell("echo {input} >> {output}")
