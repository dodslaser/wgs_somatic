# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os

rule share_to_igv:
    input:
        expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_CNV_somatic.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", vartype="somatic", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_called.seg", vartype="somatic", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", vartype="germline", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_called.seg", vartype="germline", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_REALIGNED.bam.tdf",  sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/reports/{sname}_baf.igv", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/manta/{sname}_somatic_MantaBNDs.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/manta/{sname}_germline_MantaBNDs.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/manta/{sname}_germline_MantaNOBNDs.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"])
    params:
        updateigv = pipeconfig["rules"]["share_to_igv"]["updateigv"],
        igvdatadir = pipeconfig["rules"]["share_to_igv"]["igvdatadir"]
    output:
        "reporting/shared_igv_files.txt" 
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
