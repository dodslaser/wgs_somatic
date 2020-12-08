# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule dnascope:
    input:
        normalbam = expand("normal/dedup/{normalname}_DEDUP.bam", normalname=normalname),
    params:
        threads = clusterconf["dnascope"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = config["rules"]["dnascope"]["dbsnp"],
        model = config["rules"]["dnascope"]["modelpath"],
        callsettings = config["rules"]["dnascope"]["settings"]
    output:
        "dnascope/{normalname}_DNAscope.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input.normalbam} --algo DNAscope -d {params.dbsnp} --var_type snp,indel --model {params.model} {params.callsettings} {output}")
        
rule dnascope_modelfilter:
    input:
        "dnascope/{normalname}_DNAscope.vcf"
    params:
        threads = clusterconf["dnascope_modelfilter"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        model = config["rules"]["dnascope"]["modelpath"],
    output:
        "dnascope/{normalname}_DNAscope_modelfiltered.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} --algo DNAModelApply --model {params.model} -v {input} {output}")

rule dnascope_vcffilter:
    input:
        "dnascope/{normalname}_DNAscope_modelfiltered.vcf"
    params:
        threads = clusterconf["dnascope_vcffilter"]["threads"],
        bcftools = config["rules"]["dnascope_vcffilter"]["bcftools"]
    output:
        "dnascope/{normalname}_germline.vcf"
    run:
        shell("{params.bcftools} filter -s 'ML_FAIL' -i 'INFO/ML_PROB <= 0.95' -m x dnascope/{normalname}_DNAscope_modelfiltered.vcf >> dnascope/{normalname}_DNAscope_modelfiltered_0.95.vcf")
        shell("{params.bcftools} filter -i 'FILTER=PASS' -m x dnascope/{normalname}_DNAscope_modelfiltered_0.95.vcf >> dnascope/{normalname}_germline.vcf")
