# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os

def get_bedfile_path(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]['bed']

def get_bedfile_version(wcs):
    return config["insilico"][f"{wcs.insiliconame}"]['version']

def get_insiliconame(wcs):
    return config["insilico"]

rule insilico_coverage:
    input:
        bamfile = "{workingdir}/{normalid}/dedup/{sname}_DEDUP.bam"
        #bedinfo = get_insilico_info
    params:
        insilico_wrapper = pipeconfig["rules"]["insilico"]["insilico_main"],
        samtools = pipeconfig["rules"]["insilico"]["samtools"],
        pythonversion = pipeconfig["rules"]["insilico"]["python"],
        bedfile_path = get_bedfile_path,
        bedfile_version = get_bedfile_version
    output:
        "{workingdir}/{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_10x.xlsx",
        "{workingdir}/{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_20x.xlsx",
        "{workingdir}/{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_genes_below10x.xlsx",
        "{workingdir}/{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}.csv",
        "{workingdir}/{normalid}/insilico/{insiliconame}/{sname}_{insiliconame}_cov.tsv"
    run:
        os.makedirs(f'{workingdir}/{normalid}/{wildcards.sname}/insilico', exist_ok=True)
        os.makedirs(f'{workingdir}/{normalid}/{wildcards.sname}/insilico/{wildcards.insiliconame}', exist_ok=True)

        insilico_level = config["insilico"][f"{wildcards.insiliconame}"]["levels"]
        shell(f"{params.pythonversion} {params.insilico_wrapper} --bamfile {input.bamfile} --bedfile {params.bedfile_path} --version {params.bedfile_version} --outputdir {wildcards.workingdir}/{wildcards.normalid}/insilico/{wildcards.insiliconame} --annotationlevel {insilico_level} --samplename {wildcards.sname} --samtools {params.samtools}")
