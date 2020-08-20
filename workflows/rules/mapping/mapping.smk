# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_fastqpairs(wcs):
    if f"{wcs.stype}" == "normal":
        return normal_fastqpairs
    else:
        return tumor_fastqpairs

def get_fwd_pattern(wcs):
    if f"{wcs.stype}" == "normal":
        print(f"{normalfastqs}/{wcs.fastqpattern}{n_pattern_r1}")
        return f"{normalfastqs}/{wcs.fastqpattern}{n_pattern_r1}"
    else:
        print(f"{tumorfastqs}/{wcs.fastqpattern}{t_pattern_r1}")
        return f"{tumorfastqs}/{wcs.fastqpattern}{t_pattern_r1}"

def get_rev_pattern(wcs):
    if f"{wcs.stype}" == "normal":
        print(f"{normalfastqs}/{wcs.fastqpattern}{n_pattern_r2}")
        return f"{normalfastqs}/{wcs.fastqpattern}{n_pattern_r2}"
    else:
        print(f"{tumorfastqs}/{wcs.fastqpattern}{t_pattern_r2}")
        return f"{tumorfastqs}/{wcs.fastqpattern}{t_pattern_r2}"

def get_samplename(wcs):
    return sampleconfig[f"{wcs.stype}name"]

def get_mapping(wcs):
    if f"{wcs.stype}" == "normal":
        fastqpattern = normal_fastqpairs
    else:
        fastqpattern = tumor_fastqpairs
    return expand("{workingdir}/{stype}/mapping/{fastqpattern}.bam", workingdir=f"{wcs.workingdir}", stype=f"{wcs.stype}", fastqpattern=fastqpattern)

rule mapping:
    input:
        fwd = get_fwd_pattern,
        rev = get_rev_pattern
    params:
        threads = clusterconf["mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        samplename = get_samplename,
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"]
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/mapping/{fastqpattern}.bam"
    shell:
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' "
            "-t 16 {params.referencegenome} {input.fwd} {input.rev} "
        "| {params.sentieon} util sort -o {output} -t 20 --sam2bam -i -")

rule dedup:
    input:
        bamfiles = get_mapping
    output:
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.bam",
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.txt"
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
    shell:
        "shellbamfiles=$(echo {input.bamfiles} | sed 's/ / -i /g') ;"
        "{params.sentieon} driver -t {params.threads} "
            "-i {shellbamfiles} --algo LocusCollector "
            "--fun score_info {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt ;"
        "{params.sentieon} driver -t {params.threads} "
            "-i {shellbamfiles} --algo Dedup --rmdup "
            "--score_info {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt "
            "--metrics {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.txt "
            "{wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.bam"

rule realign_mapping:
    input:
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.bam"
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        mills = pipeconfig["rules"]["realign"]["mills"],
        tgenomes = pipeconfig["rules"]["realign"]["tgenomes"],
        outputdir = pipeconfig["rules"]["realign"]["outputdir"]
    output:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} -k {params.tgenomes} {output}")

rule baserecal:
    input:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = pipeconfig["rules"]["recal"]["dbsnp"],
        mills = pipeconfig["rules"]["recal"]["mills"],
        tgenomes = pipeconfig["rules"]["recal"]["tgenomes"],
        outputdir = pipeconfig["rules"]["recal"]["outputdir"]
    output:
        "{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.dbsnp} -k {params.tgenomes} {output}")
