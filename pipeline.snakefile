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

normalfastqdirs = config["normalfastqs"]
normalname = config["normalname"]
normalid = config["normalid"]

tumorfastqdirs = config["tumorfastqs"]
tumorname = config["tumorname"]
tumorid = config["tumorid"]

igvuser = config["igvuser"]
ivauser = config["ivauser"]
reference = config["reference"]

workingdir = config["workingdir"]

##################################################
# Chose Config based on Reference
# ---------------------------------------------

if reference == "hg38":
    configfilepath = "configs/config_hg38.json"
else:
    configfilepath = "configs/config_hg19.json"
#----------------------------------------------


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

####################################################
# Prepare Fastq Variables 
# -------------------------------------------------

fwdpatterns = ["_1.fastq.gz", "_R1_001.fastq.gz", "_1.fasterq", "_R1_001.fasterq"] 
revpatterns = ["_2.fastq.gz", "_R2_001.fastq.gz", "_2.fasterq", "_R2_001.fasterq"]

fastq_dict = {}
fastq_dict["normal"] = {}
fastq_dict["normal"]["fastqpair_patterns"] = {}

fastq_dict["tumor"] = {}
fastq_dict["tumor"]["fastqpair_patterns"] = {}

# Prepare Normal Fastq Variables
for normalfastqdir in normalfastqdirs:
    for fwdpattern in fwdpatterns:
        normal_fwd_fastqs = glob.glob(f"{normalfastqdir}/*{fwdpattern}")
        if normal_fwd_fastqs:
            for normal_fwd_fastq in normal_fwd_fastqs:
                fastqpair_pattern = os.path.basename(normal_fwd_fastq).replace(fwdpattern, "")
                fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern] = {}
                fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern]["fwd"] = normal_fwd_fastq
    for revpattern in revpatterns:
        normal_rev_fastqs = glob.glob(f"{normalfastqdir}/*{revpattern}")
        if normal_rev_fastqs:
            for normal_rev_fastq in normal_rev_fastqs:
                fastqpair_pattern = os.path.basename(normal_rev_fastq).replace(revpattern, "")
                fastq_dict["normal"]["fastqpair_patterns"][fastqpair_pattern]["rev"] = normal_rev_fastq

# Prepare Tumor Fastq Variables
for tumorfastqdir in tumorfastqdirs:
    for fwdpattern in fwdpatterns:
        tumor_fwd_fastqs = glob.glob(f"{tumorfastqdir}/*{fwdpattern}")
        if tumor_fwd_fastqs:
            for tumor_fwd_fastq in tumor_fwd_fastqs:
                fastqpair_pattern = os.path.basename(tumor_fwd_fastq).replace(fwdpattern, "")
                fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern] = {}
                fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern]["fwd"] = tumor_fwd_fastq
    for revpattern in revpatterns:
        tumor_rev_fastqs = glob.glob(f"{tumorfastqdir}/*{revpattern}")
        if tumor_rev_fastqs:
            for tumor_rev_fastq in tumor_rev_fastqs:
                fastqpair_pattern = os.path.basename(tumor_rev_fastq).replace(revpattern, "")
                fastq_dict["tumor"]["fastqpair_patterns"][fastqpair_pattern]["rev"] = tumor_rev_fastq 

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
include:        "workflows/rules/small_tools/ballele.smk"
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
        return expand("{workingdir}/reporting/shared_igv_files.txt", workingdir=workingdir)
    return []

def upload_germline_iva(wildcards):
    if ivauser:
        if ivauser == "testing":
            return expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}_test.txt", workingdir=workingdir, sname=normalid, stype="normal", caller="dnascope", vcftype="germline")
        else:
            return expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", workingdir=workingdir, sname=normalid, stype="normal", caller="dnascope", vcftype="germline")
    return []
def upload_somatic_iva(wildcards):
    if ivauser:
        if ivauser == "testing":
            return expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}_test.txt", workingdir=workingdir, sname=tumorid, stype="tumor", caller="tnscope", vcftype="somatic")
        else:
            return expand("{workingdir}/reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt", workingdir=workingdir, sname=tumorid, stype="tumor", caller="tnscope", vcftype="somatic")
    return []

rule all:
    input:
        get_igv_input,
        upload_somatic_iva,
        upload_germline_iva,
        expand("{workingdir}/reporting/workflow_finished.txt", workingdir=workingdir)
