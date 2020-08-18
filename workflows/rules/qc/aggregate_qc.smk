# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule excel_qc:
    input:
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalid),
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP.txt", stype=sampleconfig[normalname]["stype"], sname=normalid),
        tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid),
        normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalid),
        tumorcanvas = expand("{stype}/canvas/{sname}_CNV_somatic.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorid)
    params:
        qcsumscript = pipeconfig["rules"]["excel_qc"]["qcsumscript"]
    output:
        "qc_report/{tumorname}_qc_stats.xlsx"
    run:
        shell("{params.qcsumscript} -tc {input.tumorcov} -nc {input.normalcov} -td {input.tumordedup} -nd {input.normaldedup} -tv {input.tumorvcf} -nv {input.normalvcf} -cv {input.tumorcanvas} -o {output}") 
