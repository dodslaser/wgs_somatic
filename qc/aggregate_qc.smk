# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


rule excel_qc:
    input:
        tumorcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[tumorname]["stype"], sname=tumorname),
        normalcov = expand("{stype}/reports/{sname}_WGScov.tsv", stype=sampleconfig[normalname]["stype"], sname=normalname),
        tumordedup = expand("{stype}/dedup/{sname}_DEDUP_score.txt", stype=sampleconfig[tumorname]["stype"], sname=tumorname),
        normaldedup = expand("{stype}/dedup/{sname}_DEDUP_score.txt", stype=sampleconfig[normalname]["stype"], sname=normalname)
    params:
        qcsumscript = pipeconfig["rules"]["excel_qc"]["qcsumscript"]
    output:
        "qc_report/{tumorname}_qc_stats.xlsx"
    run:
        shell("{params.qcsumscript} -tc {input.tumorcov} -nc {input.normalcov} -td {input.tumordedup} -nd {input.normaldedup} -o {output}") 
