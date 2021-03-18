# vim: syntax=python tabstop=4 expandtab
# # coding: utf-8
import os

rule manta_germline:
    input:
        expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"], 
        annotate = pipeconfig["rules"]["manta"]["annotate"],
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
        swegendb = pipeconfig["rules"]["manta"]["swegendb"],
        gnomaddb = pipeconfig["rules"]["manta"]["gnomaddb"],
        localdb = pipeconfig["rules"]["manta"]["localdb"],
        bcftools = pipeconfig["rules"]["manta"]["bcftools"] 
    output:
        "{workingdir}/{stype}/manta/{sname}_germline_mantaSV.vcf",
        "{workingdir}/{stype}/manta/{sname}_germline_mantaSV.vcf.xlsx",
        "{workingdir}/{stype}/manta/{sname}_germline_MantaBNDs.vcf",
        "{workingdir}/{stype}/manta/{sname}_germline_MantaNOBNDs.vcf"
    run:
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/manta/runWorkflow.py"):
            shell("{params.mantaconf} --normalBam={input} --referenceFasta {params.reference} --runDir {wildcards.workingdir}/{wildcards.stype}/manta/") # Preparing Manta
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV.vcf.gz"):
            shell("{wildcards.workingdir}/{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV.vcf"):
            shell("gunzip {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV.vcf.gz")
        shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS.vcf --db {params.swegendb} --out_frq SWEFRQ --in_frq FRQ --out_occ SWEOCC --in_occ OCC > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen.vcf --db {params.gnomaddb} --out_frq GNOMAD_AF --in_frq AF --out_occ GNOMAD_AC --in_occ AC > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad.vcf --out_frq KGGFRQ --out_occ KGGOCC --sqdb {params.localdb} > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg.vcf")
        shell("{params.bcftools} filter -e 'INFO/GNOMAD_AF >= 0.03 | INFO/SWEFRQ >= 0.03 | INFO/KGGFRQ >= 0.05' {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf")
        shell("grep -Ev 'GL000|hs37d5' {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/diploidSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf")
        shell("grep -e '^#' -e 'MantaBND:' {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_germline_MantaBNDs.vcf")
        shell("grep -v 'MantaBND:' {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_germline_MantaNOBNDs.vcf")
        shell("{params.annotate} -v {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_germline_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.workingdir}/{wildcards.stype}/manta")

rule manta_somatic:
    input:
        tumorbam = expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"])
    params:
        reference = pipeconfig["referencegenome"],
        svdb = pipeconfig["rules"]["manta"]["svdb"],
        mantaconf = pipeconfig["rules"]["manta"]["mantaconf"],
        annotate = pipeconfig["rules"]["manta"]["annotate"],
        annotate_ref = pipeconfig["rules"]["manta"]["annotate_ref"],
        swegendb = pipeconfig["rules"]["manta"]["swegendb"],
        gnomaddb = pipeconfig["rules"]["manta"]["gnomaddb"],
        localdb = pipeconfig["rules"]["manta"]["localdb"],
        bcftools = pipeconfig["rules"]["manta"]["bcftools"]
    output:
        "{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf",
        "{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx",
        "{workingdir}/{stype}/manta/{sname}_somatic_MantaBNDs.vcf",
        "{workingdir}/{stype}/manta/{sname}_somatic_MantaNOBNDs.vcf"
    run:
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/manta/runWorkflow.py"):
            shell("{params.mantaconf} --tumorBam={input.tumorbam} --normalBam={input.normalbam} --referenceFasta {params.reference} --runDir {wildcards.workingdir}/{wildcards.stype}/manta/") # Preparing Manta
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV.vcf.gz"):
            shell("{wildcards.workingdir}/{wildcards.stype}/manta/runWorkflow.py -m local") #Running Manta
        if not os.path.isfile(f"{wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV.vcf"):
            shell("gunzip {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV.vcf.gz")
        shell("grep -e $'\t'PASS$'\t' -e '^#' {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS.vcf --db {params.swegendb} --out_frq SWEFRQ --in_frq FRQ --out_occ SWEOCC --in_occ OCC > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen.vcf --db {params.gnomaddb} --out_frq GNOMAD_AF --in_frq AF --out_occ GNOMAD_AC --in_occ AC > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad.vcf")
        shell("{params.svdb} --query --query_vcf {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad.vcf --out_frq KGGFRQ --out_occ KGGOCC --sqdb {params.localdb} > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg.vcf")
        shell("{params.bcftools} filter -e 'INFO/GNOMAD_AF >= 0.03 | INFO/SWEFRQ >= 0.03 | INFO/KGGFRQ >= 0.05' {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf")
        shell("grep -Ev 'GL000|hs37d5' {wildcards.workingdir}/{wildcards.stype}/manta/results/variants/somaticSV_PASS_swegen_gnomad_kgg_freqfiltered.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf")
        shell("grep -e '^#' -e 'MantaBND:' {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_somatic_MantaBNDs.vcf")
        shell("grep -v 'MantaBND:' {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf > {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_somatic_MantaNOBNDs.vcf")
        shell("{params.annotate} -v {wildcards.workingdir}/{wildcards.stype}/manta/{wildcards.sname}_somatic_mantaSV.vcf -g {params.annotate_ref} -o {wildcards.workingdir}/{wildcards.stype}/manta")

rule manta_summary:
    input:
        "{workingdir}/{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx"
    output:
        "{workingdir}/{stype}/manta/{sname}_somatic_mantaSV_Summary.xlsx"
    run:
        manta_summary(input, output, tumorname, normalname)


