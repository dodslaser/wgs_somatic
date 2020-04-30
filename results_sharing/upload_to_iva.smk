# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from helpers import read_passconfig
from shutil import copyfile

rule filter_variants_in_bed:
    input:
        "{stype}/{caller}/{sname}_{vcftype}.vcf"    
    params:
        rtg = pipeconfig["rules"]["filter_variants_in_bed"]["rtg"],
        bedfile = pipeconfig["rules"]["filter_variants_in_bed"]["bedfile"],
    output:
        "{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf"
    run:
        if not os.path.isfile(f"{output}.gz"):
            shell("{params.rtg} vcffilter --include-bed={params.bedfile} --output={output} --input={input}")
        if os.path.isfile(f"{output}.gz"):
            shell("gunzip {output}.gz")


rule upload_to_iva:
    input:
        "{stype}/{caller}/{sname}_{vcftype}_refseq3kfilt.vcf"     
    params:
        clcserver = pipeconfig["rules"]["upload_to_iva"]["clcserver"],
        clcport = pipeconfig["rules"]["upload_to_iva"]["clcport"],
        clcqueue = pipeconfig["rules"]["upload_to_iva"]["clcqueue"],
        clcuser = pipeconfig["rules"]["upload_to_iva"]["clcuser"],
        clccmd = pipeconfig["rules"]["upload_to_iva"]["clccmd"],
        passwords = pipeconfig["rules"]["upload_to_iva"]["passwords"],
        clusterclcdir = pipeconfig["rules"]["upload_to_iva"]["clusterclcdir"],
        clcref = pipeconfig["rules"]["upload_to_iva"]["clcref"],
        clcivadir = pipeconfig["rules"]["upload_to_iva"]["clcivadir"],
        clcivadir_servpath = pipeconfig["rules"]["upload_to_iva"]["clcivadir_servpath"]
    output:
        "reporting/uploaded_to_iva_{stype}_{caller}_{sname}_{vcftype}.txt"
    run:
        # copy to clc accessible directory
        vcfbase = os.path.basename(f"{input}")
        passconfig = read_passconfig()
        clcpass = passconfig["clcpass"]
        ivapass = passconfig[ivauser]["ivapass"]
        ivaemail = passconfig[ivauser]["ivauser"]

        if not os.path.isfile(f"{params.clusterclcdir}/{vcfbase}"):
            copyfile(f"{input}", f"{params.clusterclcdir}/{vcfbase}")
        # upload to clc
        clcvcfname = vcfbase.replace(".vcf", f"_{wildcards.sname}")
        clcvcfpath = f"{params.clcivadir_servpath}/{clcvcfname}.clc"
        if not os.path.isfile(clcvcfpath):
            shell("echo {params.clccmd} -S {params.clcserver} -P {params.clcport} -U {params.clcuser} -W {clcpass} -G {params.clcqueue} -A import_tracks --type VCF --reference-track {params.clcref} -f clc://serverfile{params.clusterclcdir}/{vcfbase} -d {params.clcivadir}")
            shell("{params.clccmd} -S {params.clcserver} -P {params.clcport} -U {params.clcuser} -W {clcpass} -G {params.clcqueue} -A import_tracks --type VCF --reference-track {params.clcref} -f clc://serverfile{params.clusterclcdir}/{vcfbase} -d {params.clcivadir}")
        
        # upload from clc to iva
        shell("echo {params.clccmd} -S {params.clcserver} -P {params.clcport} -U {params.clcuser} -W {clcpass} -G {params.clcqueue} -A iva --reference {params.clcref} --study-person-email {ivaemail} --study-person-password {ivapass} --analysis-pipeline NONE -i {params.clcivadir}/{clcvcfname}")
        shell("{params.clccmd} -S {params.clcserver} -P {params.clcport} -U {params.clcuser} -W {clcpass} -G {params.clcqueue} -A iva --reference {params.clcref} --study-person-email {ivaemail} --study-person-password {ivapass} --analysis-pipeline NONE -i {params.clcivadir}/{clcvcfname}")
        shell("echo {input} >> {output}")
