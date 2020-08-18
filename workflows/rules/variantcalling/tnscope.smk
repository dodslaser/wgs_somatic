# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule tnscope:
    input:
        tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalid, stype=sampleconfig[normalname]["stype"]),
        tumortable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normaltable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        threads = clusterconf["tnscope"]["threads"],
        tumorname = tumorname,
        normalname = normalname,
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = pipeconfig["rules"]["tnscope"]["dbsnp"],
        callsettings = pipeconfig["rules"]["tnscope"]["settings"],
    output:
        "{stype}/tnscope/{sname}_TNscope_tn.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input.tumorbam} -q {input.tumortable} -i {input.normalbam} -q {input.normaltable} --algo TNscope --tumor_sample {params.tumorname} --normal_sample {params.normalname} {params.callsettings} {output}")

rule tnscope_modelfilter:
    input:
        tnscopevcf = "{stype}/tnscope/{sname}_TNscope_tn.vcf"
    params:
        threads = clusterconf["tnscope_modelfilter"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        modelpath = pipeconfig["rules"]["tnscope_modelfilter"]["modelpath"],
    output:
        "{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} --algo TNModelApply -m {params.modelpath} -v {input.tnscopevcf} {output}")

rule tnscope_vcffilter:
    input:
        tnscopevcf_ml = "{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"
    params:
        threads = clusterconf["tnscope_vcffilter"]["threads"],
        outputdir = pipeconfig["rules"]["tnscope_vcffilter"]["outputdir"],
        bcftools = pipeconfig["rules"]["tnscope_vcffilter"]["bcftools"]
    output:
        somatic_n = "{stype}/tnscope/{sname}_somatic_w_normal.vcf",
        somatic = "{stype}/tnscope/{sname}_somatic.vcf"
    run:
        vcfname = os.path.basename(f"{input.tnscopevcf_ml}")
        vcfname = vcfname.replace(".vcf", "")
        shell_command = [f"{params.bcftools} filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml}", f"> {params.outputdir}/{vcfname}_lowqual1.vcf ;",
                        f"{params.bcftools} annotate -x FILTER/triallelic_site {params.outputdir}/{vcfname}_lowqual1.vcf", f"> {params.outputdir}/{vcfname}_triallelic2.vcf ;",
                        f"{params.bcftools} annotate -x FILTER/alt_allele_in_normal {params.outputdir}/{vcfname}_triallelic2.vcf", f"> {params.outputdir}/{vcfname}_altalleleinnormal4.vcf ;",
                        f"{params.bcftools} filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {params.outputdir}/{vcfname}_altalleleinnormal4.vcf", f"> {params.outputdir}/{vcfname}_uncertainaf6.vcf ;",
                        f"{params.bcftools} filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AD[1:1]>1' -m + {params.outputdir}/{vcfname}_uncertainaf6.vcf", f"> {params.outputdir}/{vcfname}_likelyartifact7.vcf ;",
                        f"{params.bcftools} filter -s MLrejected -e 'INFO/ML_PROB<0.37' -m + {params.outputdir}/{vcfname}_likelyartifact7.vcf", f"> {params.outputdir}/{vcfname}_mladjusted8.vcf ;",
                        f"{params.bcftools} filter -i 'FILTER=\"PASS\"' {params.outputdir}/{vcfname}_likelyartifact7.vcf > {output.somatic_n} ;",
                        f"{params.bcftools} view -s {tumorname} {output.somatic_n} > {output.somatic}"]
        shell_command = " ".join(shell_command)
        print(shell_command)      
        shell(shell_command)
