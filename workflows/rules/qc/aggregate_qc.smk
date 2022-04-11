# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main

if tumorid:
    if normalid:
        rule excel_qc:
            input:
                tumorcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normalcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumordedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normaldedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumorvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normalvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumorcanvas = expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid)
            output:
                "{workingdir}/qc_report/{tumorname}_qc_stats.xlsx"
            run:
                create_excel_main(f"{input.tumorcov}", f"{input.normalcov}", f"{input.tumordedup}", f"{input.normaldedup}", f"{input.tumorvcf}", f"{input.normalvcf}", f"{input.tumorcanvas}", f"{output}") 
    else:
        rule excel_qc:
            input:
                tumorcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumordedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumorvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                #tumorcanvas = expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid)
            output:
                "{workingdir}/qc_report/{tumorname}_qc_stats.xlsx"
            run:
                #create_excel_main(tumorcov = f"{input.tumorcov}", tumordedup = f"{input.tumordedup}", tumorvcf = f"{input.tumorvcf}", f"{input.tumorcanvas}", output = f"{output}")
                create_excel_main(tumorcov = f"{input.tumorcov}", tumordedup = f"{input.tumordedup}", tumorvcf = f"{input.tumorvcf}", output = f"{output}")

else:
    rule excel_qc:
        input:
            normalcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
            normaldedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
            normalvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
        output:
            "{workingdir}/qc_report/{normalname}_qc_stats.xlsx"
        run:
            create_excel_main(normalcov = f"{input.normalcov}", normaldedup = f"{input.normaldedup}", normalvcf =f"{input.normalvcf}", output = f"{output}")
