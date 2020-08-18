# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
from os.path import join
import glob
import time
from pathlib import Path
import yaml
import helpers
import os

__author__ = "Rickard 'Ricksy' Rickardsson"

normalfastqs = config["normalfastqs"]
normalname = config["normalname"]
normalid = config["normalid"]

tumorfastqs = config["tumorfastqs"]
tumorname = config["tumorname"]
tumorid = config["tumorid"]

igvuser = config["igvuser"]
ivauser = config["ivauser"]
reference = config["reference"]

workingdir = config["workingdir"]

##########################################
# PetaGene EnvVariables       ( Does not appear to do anything...)
#os.environ["PETASUITE_REFPATH"] = "/usr/lib/petalink.so"
#os.environ["LD_PRELOAD"] = "/seqstore/software/petagene/corpus:/opt/petagene/petasuite/species"
##########################################

if reference == "hg38":
    configfilepath = "configs/config_hg38.json"
else:
    configfilepath = "configs/config_hg19.json"

pipeconfig = helpers.read_config(configfilepath)
clusterconf = helpers.read_clusterconf()

shell.executable("/bin/bash")

analysistime = time.strftime("%Y-%m-%d-%H-%M-%S")

sampleconfig = {}
sampleconfig[normalname] = {}
sampleconfig[normalname]["stype"] = "normal"
sampleconfig[tumorname] = {}
sampleconfig[tumorname]["stype"] = "tumor"
sampleconfig["normal"] = normalid
sampleconfig["tumor"] = tumorid
sampleconfig["normalname"] = normalname
sampleconfig["tumorname"] = tumorname

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

print(tumor_fastqpairs)
print(normal_fastqpairs)

sentieon = pipeconfig["sentieon"]
referencegenome = pipeconfig["referencegenome"]

wildcard_constraints:
    sname="[^_]*_[^_]*_[^_]*"

###########################################################
# Defining Non Cluster Rules
localrules: all, upload_to_iva, share_to_igv, tn_workflow, share_to_resultdir, excel_qc
###########################################################

########################################
# Workflows
include:        "workflows/tn_workflow.smk"

########################################
# Mapping
include:        "workflows/rules/mapping/generate_tdf.smk"

#########################################
# VariantCalling
include:        "workflows/rules/variantcalling/tnscope.smk"
include:        "workflows/rules/variantcalling/dnascope.smk"
include:        "workflows/rules/misc_tools/ballele.smk"
include:        "workflows/rules/variantcalling/canvas.smk"

#########################################
# QC
include:        "workflows/rules/qc/coverage.smk"
include:        "workflows/rules/qc/aggregate_qc.smk"

#########################################
# ResultSharing:
include:        "workflows/rules/results_sharing/share_to_igv.smk"
include:        "workflows/rules/results_sharing/share_to_resultdir.smk"
include:        "workflows/rules/results_sharing/upload_to_iva.smk"



if reference == "hg38":
    ###########################################################
    # HG38 rules
    ###########################################################
    # Mapping
    include:    "workflows/rules/mapping/mapping_hg38.smk"
    # Variantcalling
    include:    "workflows/rules/variantcalling/manta_hg38.smk"
else:
    ###########################################################
    # HG19 rules
    ###########################################################
    # Mapping
    include:        "workflows/rules/mapping/mapping.smk"
    # VariantCalling
    include:        "workflows/rules/variantcalling/manta.smk"

def get_igv_input(wildcards):
    if igvuser:
        return expand("{workingdir}/reporting/shared_igv_files.txt", workingdir=workingdir),
    return []

def get_iva_input(wildcards):
    input_list = []
    if ivauser:
        if ivauser == "testing":
            input_list.append(expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}_test.txt", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"], caller="dnascope", vcftype="germline"))
            input_list.append(expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}_test.txt", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"], caller="tnscope", vcftype="somatic"))
        else:
            input_list.append(expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"], caller="dnascope", vcftype="germline"))
            input_list.append(expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"], caller="tnscope", vcftype="somatic"))
         return ",".join(input_list)
    return []

rule all:
    input:
        get_igv_input
        get_iva_input
        expand("{workingdir}/reporting/workflow_finished.txt", workingdir=workingdir)
