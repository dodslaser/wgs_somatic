# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule tn_workflow:
    input:
        expand("{stype}/tnscope/{sname}_somatic.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/dnascope/{sname}_germline.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/dnascope/{sname}_germline.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_baf.igv", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_WGScov.tsv", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_CNV_germline.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_CNV_somatic.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", vartype="somatic", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_called.seg", vartype="somatic", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", vartype="germline", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/canvas/{sname}_{vartype}_CNV_called.seg", vartype="germline", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        expand("{stype}/reports/{sname}_REALIGNED.bam.tdf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        expand("{stype}/reports/{sname}_REALIGNED.bam.tdf",  sname=normalname, stype=sampleconfig[normalname]["stype"]),
        "reporting/shared_result_files.txt"
    output:
        "reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")
