# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from misc_tools.gender import calc_gender
from misc_tools.create_segfile import create_seg
from misc_tools.fix_sexploidyfile import mod_sex_vcf

rule canvas_somatic:
    input:
        germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        somatic_vcf = "{stype}/tnscope/{sname}_somatic.vcf",
        bamfile = "{stype}/realign/{sname}_REALIGNED.bam",
        normal_wgscov = expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        normal_ycov = expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"])
    params:
        annotate = pipeconfig["rules"]["canvas"]["annotate"],
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
        canvasdll = pipeconfig["rules"]["canvas"]["canvasdll"],
        dotnet = pipeconfig["rules"]["canvas"]["dotnet"],
        genomedir = pipeconfig["rules"]["canvas"]["genomedir"],
        reference = pipeconfig["rules"]["canvas"]["reference"],
        filter13 = pipeconfig["rules"]["canvas"]["filter13"],
        malevcf = pipeconfig["rules"]["canvas"]["malevcf"],
        femalevcf = pipeconfig["rules"]["canvas"]["femalevcf"]
    output:
        "{stype}/canvas/{sname}_CNV_somatic.vcf",
        "{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx"
    run:
        calculated_gender = calc_gender(input.normal_wgscov[0], input.normal_ycov[0])
        dotnetdir = os.path.dirname(params.dotnet)
        if calculated_gender == "male":
            ploidyvcf = params.malevcf
            ploidyvcf = mod_sex_vcf(ploidyvcf, f"{wildcards.sname}", f"{wildcards.stype}/canvas/")
        else:
            ploidyvcf = params.femalevcf
            ploidyvcf = mod_sex_vcf(ploidyvcf, f"{wildcards.sname}", f"{wildcards.stype}/canvas/")
        shell("export PATH={dotnetdir}:$PATH ; {params.dotnet} {params.canvasdll} Somatic-WGS --bam={input.bamfile} --somatic-vcf={input.somatic_vcf} --sample-b-allele-vcf={input.germline_snv_vcf} --sample-name={wildcards.sname} --ploidy-vcf={ploidyvcf} -g {params.genomedir} -r {params.reference} -f {params.filter13} -o {wildcards.stype}/canvas/")
        if os.path.isfile(f"{wildcards.stype}/canvas/CNV.vcf.gz"):
            shell("gunzip {wildcards.stype}/canvas/CNV.vcf.gz")
        shell("grep -v 'Canvas:REF' {wildcards.stype}/canvas/CNV.vcf > {wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf")
        shell("{params.annotate} -v {wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf -g {params.annotate_ref} -o {wildcards.stype}/canvas/")
        #os.rename(f"{wildcards.stype}/canvas/CNV.vcf", "{wildcards.stype}/canvas/{wildcards.sname}_CNV_somatic.vcf")
#        create_segfile(f"{wildcards.stype}/canvas", analysistime, f"{wildcards.sname}")
 
rule canvas_germline:
    input:
        germline_snv_vcf = expand("{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        bamfile = "{stype}/realign/{sname}_REALIGNED.bam",
        normal_wgscov = expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        normal_ycov = expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"])
    params:
        annotate = pipeconfig["rules"]["canvas"]["annotate"],
        annotate_ref = pipeconfig["rules"]["canvas"]["annotate_ref"],
        canvasdll = pipeconfig["rules"]["canvas"]["canvasdll"],
        dotnet = pipeconfig["rules"]["canvas"]["dotnet"],
        genomedir = pipeconfig["rules"]["canvas"]["genomedir"],
        reference = pipeconfig["rules"]["canvas"]["reference"],
        filter13 = pipeconfig["rules"]["canvas"]["filter13"],
        malevcf = pipeconfig["rules"]["canvas"]["malevcf"],
        femalevcf = pipeconfig["rules"]["canvas"]["femalevcf"]
    output:
        "{stype}/canvas/{sname}_CNV_germline.vcf",
        "{stype}/canvas/{sname}_CNV_germline.vcf.xlsx"
    run:
        calculated_gender = calc_gender(input.normal_wgscov[0], input.normal_ycov[0])
        dotnetdir = os.path.dirname(params.dotnet)
        if calculated_gender == "male":
            ploidyvcf = params.malevcf
            ploidyvcf = mod_sex_vcf(ploidyvcf, f"{wildcards.sname}", f"{wildcards.stype}/canvas/")
        else:
            ploidyvcf = params.femalevcf
            ploidyvcf = mod_sex_vcf(ploidyvcf, f"{wildcards.sname}", f"{wildcards.stype}/canvas/")
        shell("export PATH={dotnetdir}:$PATH ; {params.dotnet} {params.canvasdll} SmallPedigree-WGS --bam={input.bamfile} --sample-b-allele-vcf={input.germline_snv_vcf} --ploidy-vcf={ploidyvcf} -g {params.genomedir} -r {params.reference} -f {params.filter13} -o {wildcards.stype}/canvas/")
        if os.path.isfile(f"{wildcards.stype}/canvas/CNV.vcf.gz"):
            shell("gunzip {wildcards.stype}/canvas/CNV.vcf.gz")
        shell("grep -v 'Canvas:REF' {wildcards.stype}/canvas/CNV.vcf > {wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf")
        shell("{params.annotate} -v {wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf -g {params.annotate_ref} -o {wildcards.stype}/canvas/")
        #os.rename(f"{wildcards.stype}/canvas/CNV.vcf", f"{wildcards.stype}/canvas/{wildcards.sname}_CNV_germline.vcf")
#        create_segfile(f"{wildcards.stype}/canvas", analysistime, f"{wildcards.sname}")

rule createseg:
    input:
        "{stype}/canvas/{sname}_CNV_{vartype}.vcf"
    output:
        "{stype}/canvas/{sname}_{vartype}_CNV_observed.seg",
        "{stype}/canvas/{sname}_{vartype}_CNV_called.seg"
    run:
        #canvasname = f"{wildcards.sname}_{wildcards.vartype}"
        create_seg(f"{wildcards.stype}/canvas", f"{wildcards.sname}", f"{wildcards.vartype}")
        

