# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from os.path import join
import glob
import time
from pathlib import Path
import yaml
import helpers

__author__ = "Rickard 'Ricksy' Rickardsson"

###########################################################
# Defining Non Cluster Rules
localrules: all, upload_to_iva, share_to_igv
##########################################################


pipeconfig = helpers.read_config()
#inputdict = helpers.read_inputfile()
clusterconf = helpers.read_clusterconf()

shell.executable("/bin/bash")

analysistime = time.strftime("%Y-%m-%d-%H-%M-%S")

normalfastqs = config["normalfastqs"]
normalname = config["normalname"]
tumorfastqs = config["tumorfastqs"]
tumorname = config["tumorname"]
igvuser = config["igvuser"]
ivauser = config["ivauser"]

sampleconfig = {}
sampleconfig[normalname] = {}
sampleconfig[normalname]["stype"] = "normal"
sampleconfig[tumorname] = {}
sampleconfig[tumorname]["stype"] = "tumor"
sampleconfig["normal"] = normalname
sampleconfig["tumor"] = tumorname

normal_fastqpairs, = glob_wildcards(normalfastqs + "/{normal_fastqpairs}_R1_001.fastq.gz")
tumor_fastqpairs, = glob_wildcards(tumorfastqs + "/{tumor_fastqpairs}_R1_001.fastq.gz")

if not normal_fastqpairs:
    normal_fastqpairs, = glob_wildcards(normalfastqs + "/{normal_fastqpairs}_1.fastq.gz")
    n_pattern_r1 = '_1.fastq.gz'
    n_pattern_r2 = '_2.fastq.gz'
else:
    n_pattern_r1 = '_R1_001.fastq.gz'
    n_pattern_r2 = '_R2_001.fastq.gz'

if not tumor_fastqpairs:
    tumor_fastqpairs, = glob_wildcards(tumorfastqs + "/{tumor_fastqpairs}_1.fastq.gz")
    t_pattern_r1 = '_1.fastq.gz'
    t_pattern_r2 = '_2.fastq.gz'
else:
    t_pattern_r1 = '_R1_001.fastq.gz'
    t_pattern_r2 = '_R2_001.fastq.gz'

sentieon = pipeconfig["sentieon"]
referencegenome = pipeconfig["referencegenome"]

if igvuser and not ivauser:
    rule all:
        input:
            expand("{stype}/tnscope/{sname}_somatic.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/dnascope/{sname}_germline.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/dnascope/{sname}_germline.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/reports/{sname}_WGScov.tsv", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            "reporting/shared_igv_files.txt",
            expand("{stype}/manta/{sname}_germline_mantaSV.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_germline_mantaSV.vcf.xlsx", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_somatic_mantaSV.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalname, stype=sampleconfig[normalname]["stype"])

elif igvuser and ivauser:
    rule all:
        input:
            expand("{stype}/tnscope/{sname}_somatic.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/dnascope/{sname}_germline.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/dnascope/{sname}_germline.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/reports/{sname}_WGScov.tsv", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/reports/{sname}_WGScov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/reports/{sname}_Ycov.tsv", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            "reporting/shared_igv_files.txt",
            expand("reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", sname=normalname, stype=sampleconfig[normalname]["stype"], caller="dnascope", vcftype="germline"),
            expand("reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", sname=tumorname, stype=sampleconfig[tumorname]["stype"], caller="tnscope", vcftype="somatic"),
            expand("{stype}/manta/{sname}_germline_mantaSV.vcf", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_germline_mantaSV.vcf.xlsx", sname=normalname, stype=sampleconfig[normalname]["stype"]),
            expand("{stype}/manta/{sname}_somatic_mantaSV.vcf", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/manta/{sname}_somatic_mantaSV.vcf.xlsx", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/canvas/{sname}_CNV_somatic.vcf.xlsx", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
            expand("{stype}/canvas/{sname}_CNV_germline.vcf.xlsx", sname=normalname, stype=sampleconfig[normalname]["stype"])

else:
    rule all:
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
            expand("{stype}/reports/{sname}_REALIGNED.bam.tdf",  sname=normalname, stype=sampleconfig[normalname]["stype"])


########################################
# Mapping 
include:        "mapping/mapping.smk"
include:        "mapping/generate_tdf.smk"

########################################
# VariantCalling
include:        "variantcalling/tnscope.smk"
include:        "variantcalling/dnascope.smk"
include:        "misc_tools/ballele.smk"
include:        "variantcalling/canvas.smk"
include:        "variantcalling/manta.smk"

#########################################
# QC
include:        "qc/coverage.smk"

#########################################
# ResultSharing:
include:        "results_sharing/share_to_igv.smk"
include:        "results_sharing/upload_to_iva.smk"
