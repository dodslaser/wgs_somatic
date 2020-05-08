# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def get_fastqpairs(wcs):
    if f"{wcs.stype}" == "normal":
        return normal_fastqpairs
    else:
        return tumor_fastqpairs

def get_fwd_pattern(wcs):
    if f"{wcs.stype}" == "normal":
        return f"{normalfastqs}/{wcs.fastqpattern}{n_pattern_r1}"
    else:
        return f"{tumorfastqs}/{wcs.fastqpattern}{t_pattern_r1}"

def get_rev_pattern(wcs):
    if f"{wcs.stype}" == "normal":
        return f"{normalfastqs}/{wcs.fastqpattern}{n_pattern_r2}"
    else:
        return f"{tumorfastqs}/{wcs.fastqpattern}{t_pattern_r2}"

def get_samplename(wcs):
    return sampleconfig[f"{wcs.stype}"]

def get_mapping(wcs):
    if f"{wcs.stype}" == "normal":
        fastqpattern = normal_fastqpairs
    else:
        fastqpattern = tumor_fastqpairs
    return expand("{stype}/mapping/{fastqpattern}.bam", stype=f"{wcs.stype}", fastqpattern=fastqpattern)

rule mapping:
    input:
        fwd = get_fwd_pattern,
        rev = get_rev_pattern
    params:
        threads = clusterconf["mapping"]["threads"],
        sentieon = sentieon,
        samplename = get_samplename,
        referencegenome = referencegenome,
        outputdir = pipeconfig["rules"]["mapping"]["outputdir"]
    output:
        "{stype}/mapping/{fastqpattern}.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} bwa mem -M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' -t 16 {params.referencegenome} {input.fwd} {input.rev} | {params.sentieon} util sort -o {output} -t 20 --sam2bam -i -")

rule dedup:
    input:
        bamfiles = get_mapping
    output:
        "{stype}/dedup/{sname}_DEDUP.bam",
        "{stype}/dedup/{sname}_DEDUP_score.txt"
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = sentieon,
        outputdir = pipeconfig["rules"]["mapping"]["outputdir"]
    run:
        inp_bamfiles = ""
        for bamfile in input.bamfiles:
            inp_bamfiles = f"{inp_bamfiles}-i {bamfile} "
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} {inp_bamfiles}--algo LocusCollector --fun score_info {wildcards.stype}/dedup/{params.samplename}_DEDUP_score.txt")
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} {inp_bamfiles}--algo Dedup --rmdup --score_info {wildcards.stype}/dedup/{params.samplename}_DEDUP_score.txt --metrics {wildcards.stype}/dedup/{params.samplename}_DEDUP.txt {wildcards.stype}/dedup/{params.samplename}_DEDUP.bam")

rule realign_mapping:
    input:
        "{stype}/dedup/{sname}_DEDUP.bam"
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        mills = pipeconfig["rules"]["realign"]["mills"],
        tgenomes = pipeconfig["rules"]["realign"]["tgenomes"],
        outputdir = pipeconfig["rules"]["realign"]["outputdir"]
    output:
        "{stype}/realign/{sname}_REALIGNED.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} -k {params.tgenomes} {output}")

rule baserecal:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = pipeconfig["rules"]["recal"]["dbsnp"],
        mills = pipeconfig["rules"]["recal"]["mills"],
        tgenomes = pipeconfig["rules"]["recal"]["tgenomes"],
        outputdir = pipeconfig["rules"]["recal"]["outputdir"]
    output:
        "{stype}/recal/{sname}_RECAL_DATA.TABLE"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.tgenomes} -k {params.dbsnp} {output}")
