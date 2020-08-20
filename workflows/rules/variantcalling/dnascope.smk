# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule dnascope:
    input:
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.bam"
    params:
        threads = clusterconf["dnascope"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        model = pipeconfig["singularities"]["sentieon"]["dnascope_m"],
        callsettings = pipeconfig["rules"]["dnascope"]["settings"]
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/dnascope/{sname}_DNAscope.vcf"
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} "
            "-i {input} --algo DNAscope -d {params.dbsnp} "
            "--var_type snp,indel --model {params.model} {params.callsettings} {output}"
        
rule dnascope_modelfilter:
    input:
        "{workingdir}/{stype}/dnascope/{sname}_DNAscope.vcf"
    params:
        threads = clusterconf["dnascope_modelfilter"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        model = pipeconfig["singularities"]["sentieon"]["dnascope_m"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf"
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} --algo DNAModelApply --model {params.model} -v {input} {output}"

rule dnascope_vcffilter:
    input:
        "{workingdir}/{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf"
    params:
        threads = clusterconf["dnascope_vcffilter"]["threads"],
        bcftools = pipeconfig["rules"]["dnascope_vcffilter"]["bcftools"],
        vcftools = pipeconfig["rules"]["dnascope_vcffilter"]["vcftools"],
        passfilter = "'FILTER=\"PASS\"'"
    output:
        "{workingdir}/{stype}/dnascope/{sname}_germline.vcf",
        "{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf"
    run:
        shell("{params.bcftools} filter -s 'ML_FAIL' -i 'INFO/ML_PROB <= 0.95' -m x {wildcards.workingdir}/{wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered.vcf > {wildcards.workingdir}/{wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered_0.95.vcf")
        shell("{params.bcftools} filter -i {params.passfilter} -m x {wildcards.workingdir}/{wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered_0.95.vcf > {wildcards.workingdir}/{wildcards.stype}/dnascope/{wildcards.sname}_germline.vcf")
        shell("{params.vcftools} --vcf {wildcards.workingdir}/{wildcards.stype}/dnascope/{wildcards.sname}_germline.vcf --remove-indels --recode --out {wildcards.workingdir}/{wildcards.stype}/dnascope/{wildcards.sname}_germline_SNVsOnly")
