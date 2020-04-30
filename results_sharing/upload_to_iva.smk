# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
from helpers import read_passconfig
from shutil import copyfile

rule filter_variants_in_bed:
    input:
        "{stype}/{caller}/{sname}_{vcftype}.vcf"    
    params:
        rtg = config["rules"]["filter_variants_in_bed"]["rtg"],
        bedfile = config["rules"]["filter_variants_in_bed"]["bedfile"],
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
        clcserver = config["rules"]["upload_to_iva"]["clcserver"],
        clcport = config["rules"]["upload_to_iva"]["clcport"],
        clcqueue = config["rules"]["upload_to_iva"]["clcqueue"],
        clcuser = config["rules"]["upload_to_iva"]["clcuser"],
        clccmd = config["rules"]["upload_to_iva"]["clccmd"],
        passwords = config["rules"]["upload_to_iva"]["passwords"],
        clusterclcdir = config["rules"]["upload_to_iva"]["clusterclcdir"],
        clcref = config["rules"]["upload_to_iva"]["clcref"],
        clcivadir = config["rules"]["upload_to_iva"]["clcivadir"],
        clcivadir_servpath = config["rules"]["upload_to_iva"]["clcivadir_servpath"]
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
