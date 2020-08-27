# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule tnscope:
    input:
        tumorbam = expand("tumor/realign/{tumorname}_REALIGNED.bam", tumorname=tumorname),
        normalbam = expand("normal/realign/{normalname}_REALIGNED.bam", normalname=normalname),
        tumortable = expand("tumor/recal/{tumorname}_RECAL_DATA.TABLE", tumorname=tumorname),
        normaltable = expand("normal/recal/{normalname}_RECAL_DATA.TABLE", normalname=normalname),
    params:
        threads = clusterconf["tnscope"]["threads"],
        normalname = normalname,
        tumorname = tumorname,
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = config["rules"]["tnscope"]["dbsnp"],
        callsettings = config["rules"]["tnscope"]["settings"],
        outputdir = config["rules"]["tnscope"]["outputdir"]
    output:
        "tnscope/{tumorname}_TNscope_tn.vcf"
    run:
        shell("echo 'export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input.tumorbam} -q {input.tumortable} -i {input.normalbam} -q {input.normaltable} --algo TNscope --tumor_sample {params.tumorname} --normal_sample {params.normalname} {params.callsettings} {output}' &>> {output}")

rule tnscope_modelfilter:
    input:
        tnscopevcf = "tnscope/{tumorname}_TNscope_tn.vcf"
    params:
        threads = clusterconf["tnscope_modelfilter"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        modelpath = config["rules"]["tnscope_modelfilter"]["modelpath"],
        outputdir = config["rules"]["tnscope_modelfilter"]["outputdir"]
    output:
        "tnscope/{tumorname}_TNscope_tn_ML.vcf"
    run:
        shell("echo 'export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} --algo TNModelApply -m {params.modelpath} -v {input.tnscopevcf} {output}' &>> {output}")

rule tnscope_vcffilter:
    input:
        tnscopevcf_ml = "tnscope/{tumorname}_TNscope_tn_ML.vcf"
    params:
        threads = clusterconf["tnscope_vcffilter"]["threads"],
        outputdir = config["rules"]["tnscope_vcffilter"]["outputdir"]
    output:
        "tnscope/{tumorname}_TNscope_tn_ML_filtered.vcf"
    run:
        vcfname = os.path.basename(f"{input.tnscopevcf_ml}")
        vcfname = vcfname.replace(".vcf", "")
        shell_command = ["module load bcftools/1.9 ;", 
                        f"bcftools filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml}", f">> {params.outputdir}/{vcfname}_lowqual1.vcf ;",
                        f"bcftools annotate -x FILTER/triallelic_site {params.outputdir}/{vcfname}_lowqual1.vcf", f">> /{params.outputdir}/{vcfname}_triallelic2.vcf ;",
                        f"bcftools annotate -x FILTER/alt_allele_in_normal {params.outputdir}/{vcfname}_triallelic2.vcf", f">> {params.outputdir}/{vcfname}_altalleleinnormal4.vcf ;",
                        f"bcftools filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {params.outputdir}/{vcfname}_altalleleinnormal4.vcf", f">> {params.outputdir}/{vcfname}_uncertainaf6.vcf ;",
                        f"bcftools filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AD[1:1]>1' -m + {params.outputdir}/{vcfname}_uncertainaf6.vcf", f">> {params.outputdir}/{vcfname}_likelyartifact7.vcf ;",
                        f"bcftools filter -s MLrejected -e 'INFO/ML_PROB<0.37' -m + {params.outputdir}/{vcfname}_likelyartifact7.vcf", f">> {params.outputdir}/{vcfname}_mladjusted8.vcf ;",
                        f"cp {params.outputdir}/{vcfname}_likelyartifact7.vcf {output}"]
        shell_command = " ".join(shell_command)
        shell("echo 'vcfilter' &>> {output}") 
#       print(shell_command)
#        shell(shell_command)
