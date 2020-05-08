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
localrules: all, upload_to_iva, share_to_igv, tn_workflow, share_to_resultdir
##########################################################


pipeconfig = helpers.read_config()
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
            "reporting/shared_igv_files.txt",
            "reporting/workflow_finished.txt"

elif igvuser and ivauser:
    rule all:
        input:
            "reporting/shared_igv_files.txt",
            expand("reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", sname=normalname, stype=sampleconfig[normalname]["stype"], caller="dnascope", vcftype="germline"),
            expand("reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", sname=tumorname, stype=sampleconfig[tumorname]["stype"], caller="tnscope", vcftype="somatic"),
            "reporting/workflow_finished.txt"
else:
    rule all:
        input:
            "reporting/workflow_finished.txt"

########################################
# Workflows
include:        "workflows/tn_workflow.smk"

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
include:        "qc/aggregate_qc.smk"

#########################################
# ResultSharing:
include:        "results_sharing/share_to_igv.smk"
include:        "results_sharing/share_to_resultdir.smk"
include:        "results_sharing/upload_to_iva.smk"
