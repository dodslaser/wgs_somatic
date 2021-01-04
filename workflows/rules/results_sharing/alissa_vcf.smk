# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule alissa_vcf:
input:
    "{workingdir}/{sname}_somatic_refseq3kfilt.vcf"
output:
    "{workingdir}/{sname}_somatic_refseq3kfilt_Alissa.vcf”
run:
    shell("cp {input} {output}")
    shell("python vcf_conversion.py {output}")
