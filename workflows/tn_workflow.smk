# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule tn_workflow:
    input:
        expand("{stype}/tnscope/{sname}_somatic.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/dnascope/{sname}_germline.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/dnascope/{sname}_germline.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_baf.igv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_WGScov.tsv", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_CNV_somatic.vcf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_called.seg", vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_called.seg", vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_REALIGNED.bam.tdf",  sname=normalid, stype=sampleconfig[normalname]["stype"]),
        "reporting/shared_result_files.txt"
    output:
        "reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")
