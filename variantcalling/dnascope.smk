# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule dnascope:
    input:
        "{stype}/dedup/{sname}_DEDUP.bam"
#        bam = expand("{sampledir}/dedup/{sname}_DEDUP.bam", sampledir=tumornormaldirs, sname=tumornormalnames)
    params:
        threads = clusterconf["dnascope"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = config["rules"]["dnascope"]["dbsnp"],
        model = config["rules"]["dnascope"]["modelpath"],
        callsettings = config["rules"]["dnascope"]["settings"]
    output:
        "{stype}/dnascope/{sname}_DNAscope.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo DNAscope -d {params.dbsnp} --var_type snp,indel --model {params.model} {params.callsettings} {output}")
        
rule dnascope_modelfilter:
    input:
        "{stype}/dnascope/{sname}_DNAscope.vcf"
    params:
        threads = clusterconf["dnascope_modelfilter"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        model = config["rules"]["dnascope"]["modelpath"]
    output:
        "{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} --algo DNAModelApply --model {params.model} -v {input} {output}")

rule dnascope_vcffilter:
    input:
        "{stype}/dnascope/{sname}_DNAscope_modelfiltered.vcf"
    params:
        threads = clusterconf["dnascope_vcffilter"]["threads"],
        bcftools = config["rules"]["dnascope_vcffilter"]["bcftools"],
        vcftools = config["rules"]["dnascope_vcffilter"]["vcftools"],
        passfilter = "'FILTER=\"PASS\"'"
    output:
        "{stype}/dnascope/{sname}_germline.vcf",
        "{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf"
    run:
        shell("{params.bcftools} filter -s 'ML_FAIL' -i 'INFO/ML_PROB <= 0.95' -m x {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered.vcf > {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered_0.95.vcf")
        shell("{params.bcftools} filter -i {params.passfilter} -m x {wildcards.stype}/dnascope/{wildcards.sname}_DNAscope_modelfiltered_0.95.vcf > {wildcards.stype}/dnascope/{wildcards.sname}_germline.vcf")
        shell("{params.vcftools} --vcf {wildcards.stype}/dnascope/{wildcards.sname}_germline.vcf --remove-indels --recode --out {wildcards.stype}/dnascope/{wildcards.sname}_germline_SNVsOnly")
