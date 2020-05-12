# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule excel_qc:
    input:
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorname),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalname),
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP_score.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorname),
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP_score.txt", stype=sampleconfig[normalname]["stype"], sname=normalname),
        tumorvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorname),
        normalvcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", stype=sampleconfig[normalname]["stype"], sname=normalname),
        tumorcanvas = expand("{stype}/canvas/{sname}_CNV_somatic.vcf", stype=sampleconfig[tumorname]["stype"], sname=tumorname)
    params:
        qcsumscript = pipeconfig["rules"]["excel_qc"]["qcsumscript"]
    output:
        "qc_report/{tumorname}_qc_stats.xlsx"
    run:
        shell("{params.qcsumscript} -tc {input.tumorcov} -nc {input.normalcov} -td {input.tumordedup} -nd {input.normaldedup} -tv {input.tumorvcf} -nv {input.normalvcf} -cv {input.tumorcanvas} -o {output}") 
