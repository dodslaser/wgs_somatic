# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.create_qc_excel import create_excel_main
import os

if tumorid:
    if normalid:
        # excel_qc rule for paired analysis
        rule excel_qc:
            input:
                tumorcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                ycov = expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                normalcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumordedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normaldedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                tumorvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                normalvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
                canvasvcf = expand("{workingdir}/{stype}/canvas/{sname}_CNV_somatic.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                insilicofile = expand("{workingdir}/{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], insiliconame=sampleconfig["insilico"], sname=normalid),
            output:
                "{workingdir}/qc_report/{tumorname}_qc_stats.xlsx"
            run:
                my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0]
                insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] # Somehow this stuff works
                create_excel_main(tumorcov = f"{input.tumorcov}", ycov = f"{input.ycov}" , normalcov = f"{input.normalcov}", tumordedup = f"{input.tumordedup}", normaldedup = f"{input.normaldedup}", tumorvcf = f"{input.tumorvcf}", normalvcf = f"{input.normalvcf}", canvasvcf = f"{input.canvasvcf}", output = f"{output}", insilicodir = f"{insilicodir}") 
    else:
        # excel_qc rule for tumor only
        rule excel_qc:
            input:
                tumorcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                ycov = expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumordedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                tumorvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], sname=tumorid),
                insilicofile = expand("{workingdir}/{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", workingdir=workingdir, stype=sampleconfig[tumorname]["stype"], insiliconame=sampleconfig["insilico"], sname=tumorid),
            output:
                "{workingdir}/qc_report/{tumorname}_qc_stats.xlsx"
            run:
                my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0]
                insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] # Somehow this stuff works
                create_excel_main(tumorcov = f"{input.tumorcov}", ycov = f"{input.ycov}", tumordedup = f"{input.tumordedup}", tumorvcf = f"{input.tumorvcf}", output = f"{output}", insilicodir = f"{insilicodir}")

else:
    # excel_qc rule for normal only
    rule excel_qc:
        input:
            normalcov = expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
            ycov = expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
            normaldedup = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.txt", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
            normalvcf = expand("{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], sname=normalid),
            insilicofile = expand("{workingdir}/{stype}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx", workingdir=workingdir, stype=sampleconfig[normalname]["stype"], insiliconame=sampleconfig["insilico"], sname=normalid),
        output:
            "{workingdir}/qc_report/{normalname}_qc_stats.xlsx"
        run:
            my_insilicofile = f"{input.insilicofile}".split(" ", 1)[0]
            insilicodir = os.path.dirname(f"{my_insilicofile}").rsplit("/", 1)[0] # Somehow this stuff works
            create_excel_main(normalcov = f"{input.normalcov}", ycov = f"{input.ycov}", normaldedup = f"{input.normaldedup}", normalvcf =f"{input.normalvcf}", output = f"{output}", insilicodir = f"{insilicodir}")
