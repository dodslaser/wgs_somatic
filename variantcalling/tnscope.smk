# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule tnscope:
    input:
        tumorbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{stype}/realign/{sname}_REALIGNED.bam", sname=normalname, stype=sampleconfig[normalname]["stype"]),
        tumortable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=tumorname, stype=sampleconfig[tumorname]["stype"]),
        normaltable = expand("{stype}/recal/{sname}_RECAL_DATA.TABLE", sname=normalname, stype=sampleconfig[normalname]["stype"]),
    params:
        threads = clusterconf["tnscope"]["threads"],
        tumorname = tumorname,
        normalname = normalname,
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = config["rules"]["tnscope"]["dbsnp"],
        callsettings = config["rules"]["tnscope"]["settings"],
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
        modelpath = config["rules"]["tnscope_modelfilter"]["modelpath"],
    output:
        "{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} --algo TNModelApply -m {params.modelpath} -v {input.tnscopevcf} {output}")

rule tnscope_vcffilter:
    input:
        tnscopevcf_ml = "{stype}/tnscope/{sname}_TNscope_tn_ML.vcf"
    params:
        threads = clusterconf["tnscope_vcffilter"]["threads"],
        outputdir = config["rules"]["tnscope_vcffilter"]["outputdir"],
        bcftools = config["rules"]["tnscope_vcffilter"]["bcftools"]
    output:
        "{stype}/tnscope/{sname}_somatic.vcf"
    run:
        vcfname = os.path.basename(f"{input.tnscopevcf_ml}")
        vcfname = vcfname.replace(".vcf", "")
        shell_command = [f"{params.bcftools} filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml}", f">> {params.outputdir}/{vcfname}_lowqual1.vcf ;",
                        f"{params.bcftools} annotate -x FILTER/triallelic_site {params.outputdir}/{vcfname}_lowqual1.vcf", f">> {params.outputdir}/{vcfname}_triallelic2.vcf ;",
                        f"{params.bcftools} annotate -x FILTER/alt_allele_in_normal {params.outputdir}/{vcfname}_triallelic2.vcf", f">> {params.outputdir}/{vcfname}_altalleleinnormal4.vcf ;",
                        f"{params.bcftools} filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {params.outputdir}/{vcfname}_altalleleinnormal4.vcf", f">> {params.outputdir}/{vcfname}_uncertainaf6.vcf ;",
                        f"{params.bcftools} filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AD[1:1]>1' -m + {params.outputdir}/{vcfname}_uncertainaf6.vcf", f">> {params.outputdir}/{vcfname}_likelyartifact7.vcf ;",
                        f"{params.bcftools} filter -s MLrejected -e 'INFO/ML_PROB<0.37' -m + {params.outputdir}/{vcfname}_likelyartifact7.vcf", f">> {params.outputdir}/{vcfname}_mladjusted8.vcf ;",
                        f"cp {params.outputdir}/{vcfname}_likelyartifact7.vcf {output}"]
        shell_command = " ".join(shell_command)
        print(shell_command)      
        shell(shell_command)
