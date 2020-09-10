# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule tn_workflow:
    input:
        expand("{workingdir}/{stype}/tnscope/{sname}_somatic.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/dnascope/{sname}_germline.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{workingdir}/{stype}/dnascope/{sname}_germline.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", workingdir=workingdir, vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_{vartype}_CNV_called.seg", workingdir=workingdir, vartype="somatic", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_{vartype}_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir,  sname=normalid, stype=sampleconfig[normalname]["stype"]),
        "{workingdir}/reporting/shared_result_files.txt"
    output:
        "{workingdir}/reporting/workflow_finished.txt"
    run:
        shell("echo {input} >> {output}")
