# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_tnscope_settings(wcs):
    return tnscopesetting[f"{wcs.tnsetting}"]["settings"]

rule tnscope:
    input:
        tumorbam = expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        tumortable = expand("{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normaltable = expand("{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        threads = clusterconf["tnscope"]["threads"],
        tumorname = tumorname,
        normalname = normalname,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        callsettings = get_tnscope_settings
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_TNscope_tn_{tnsetting}.vcf"
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.reference} "
            "-i {input.tumorbam} -q {input.tumortable} -i {input.normalbam} -q {input.normaltable} "
            "--algo TNscope --tumor_sample {params.tumorname} --normal_sample {params.normalname} "
            "{params.callsettings} {output}"

rule tnscope_given:
    input:
        tumorbam = expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{workingdir}/{stype}/realign/{sname}_REALIGNED.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
        tumortable = expand("{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normaltable = expand("{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
    params:
        threads = clusterconf["tnscope"]["threads"],
        tumorname = tumorname,
        normalname = normalname,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        truthset = config["data"]["tset"]
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/tnscope_given/{sname}_TNscope_tn_given.vcf"
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.reference} "
            "-i {input.tumorbam} -q {input.tumortable} -i {input.normalbam} -q {input.normaltable} "
            "--algo TNscope --tumor_sample {params.tumorname} --normal_sample {params.normalname} "
            "--given {params.truthset} {output}"

rule tnscope_modelfilter:
    input:
        tnscopevcf = "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_TNscope_tn_{tnsetting}.vcf"
    params:
        threads = clusterconf["tnscope_modelfilter"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        reference = pipeconfig["singularities"]["sentieon"]["reference"],
        modelpath = pipeconfig["singularities"]["sentieon"]["tnscope_m"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_TNscope_tn_ML_{tnsetting}.vcf"
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.reference} "
            "--algo TNModelApply -m {params.modelpath} -v {input.tnscopevcf} {output}"

rule tnscope_vcffilter:
    input:
        tnscopevcf_ml = "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_TNscope_tn_ML_{tnsetting}.vcf"
    params:
        threads = clusterconf["tnscope_vcffilter"]["threads"],
        outputdir = pipeconfig["rules"]["tnscope_vcffilter"]["outputdir"],
        bcftools = pipeconfig["rules"]["tnscope_vcffilter"]["bcftools"]
    output:
        somatic_n = "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_somatic_w_normal_{tnsetting}.vcf",
        somatic = "{workingdir}/{stype}/tnscope/{tnsetting}/{sname}_somatic_{tnsetting}.vcf"
    run:
        vcfname = os.path.basename(f"{input.tnscopevcf_ml}")
        vcfname = vcfname.replace(".vcf", "")
        shell_command = [f"{params.bcftools} filter -s LowQual -e 'QUAL < 1' -m + {input.tnscopevcf_ml}", f"> {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_lowqual1.vcf ;",
                        f"{params.bcftools} annotate -x FILTER/triallelic_site {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_lowqual1.vcf", f"> {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_triallelic2.vcf ;",
                        f"{params.bcftools} annotate -x FILTER/alt_allele_in_normal {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_triallelic2.vcf", f"> {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_altalleleinnormal4.vcf ;",
                        f"{params.bcftools} filter -s uncertainAF -e 'FORMAT/AF[0]<0.045 && FORMAT/AD[0:1]<4' -m + {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_altalleleinnormal4.vcf", f"> {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_uncertainaf6.vcf ;",
                        f"{params.bcftools} filter -s likely_artifact -e 'FORMAT/AF[0]<0.1 && FORMAT/AD[1:1]>1' -m + {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_uncertainaf6.vcf", f"> {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_likelyartifact7.vcf ;",
                        f"{params.bcftools} filter -s MLrejected -e 'INFO/ML_PROB<0.37' -m + {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_likelyartifact7.vcf", f"> {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_mladjusted8.vcf ;",
                        f"{params.bcftools} filter -i 'FILTER=\"PASS\"' {wildcards.workingdir}/{params.outputdir}/{wildcards.tnsetting}/{vcfname}_likelyartifact7.vcf > {output.somatic_n} ;",
                        f"{params.bcftools} view -s {tumorname} {output.somatic_n} > {output.somatic}"]
        shell_command = " ".join(shell_command)
        print(shell_command)      
        shell(shell_command)
