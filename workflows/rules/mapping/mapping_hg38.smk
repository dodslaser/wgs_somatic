# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

<<<<<<< HEAD
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
=======
def get_fwd_pattern(wcs):
    return fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["fwd"] 

def get_rev_pattern(wcs):
    return fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["rev"]

def format_fwd(wcs):
    fastq = fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["rev"]
    if fastq.endswith(".fasterq"):
        fastq = fastq.replace(".fasterq", "fastq.gz")
    return fastq

def format_rev(wcs):
    fastq = fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"][f"{wcs.fastqpattern}"]["fwd"]
    if fastq.endswith(".fasterq"):
        fastq = fastq.replace(".fasterq", "fastq.gz")
    return fastq

def get_samplename(wcs):
    return sampleconfig[f"{wcs.stype}name"]

def get_mapping(wcs):
    fastqpatterns = []
    for fastqpattern in fastq_dict[f"{wcs.stype}"]["fastqpair_patterns"]:
        fastqpatterns.append(fastqpattern)
    return expand("{workingdir}/{stype}/mapping/{fastqpattern}.bam", workingdir=f"{wcs.workingdir}", stype=f"{wcs.stype}", fastqpattern=fastqpatterns)
>>>>>>> sentieon_singularity

rule mapping:
    input:
        fwd = get_fwd_pattern,
        rev = get_rev_pattern
    params:
        threads = clusterconf["mapping"]["threads"],
<<<<<<< HEAD
        sentieon = sentieon,
        samplename = get_samplename,
        referencegenome = referencegenome,
        outputdir = pipeconfig["rules"]["mapping"]["outputdir"]
    output:
        "{stype}/mapping/{fastqpattern}.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} bwa mem -M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' -t 16 {params.referencegenome} {input.fwd} {input.rev} | {params.sentieon} util sort -o {output} -t 20 --sam2bam -i -")
=======
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        samplename = get_samplename,
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        fwd_fmt = format_fwd,
        rev_fmt = format_rev
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/mapping/{fastqpattern}.bam"
    shell:
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.fastqpattern}\\tSM:{params.samplename}\\tPL:ILLUMINA' "
            "-t {params.threads} {params.referencegenome} {params.fwd_fmt} {params.rev_fmt} "
        "| {params.sentieon} util sort -o {output} -t {params.threads} --sam2bam -i -"
>>>>>>> sentieon_singularity

rule dedup:
    input:
        bamfiles = get_mapping
    output:
<<<<<<< HEAD
        "{stype}/dedup/{sname}_DEDUP.bam",
        "{stype}/dedup/{sname}_DEDUP.txt"
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
        outputdir = pipeconfig["rules"]["realign"]["outputdir"]
    output:
        "{stype}/realign/{sname}_REALIGNED.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} {output}")

rule baserecal:
    input:
        "{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = pipeconfig["rules"]["recal"]["dbsnp"],
        mills = pipeconfig["rules"]["recal"]["mills"],
        outputdir = pipeconfig["rules"]["recal"]["outputdir"]
    output:
        "{stype}/recal/{sname}_RECAL_DATA.TABLE"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.dbsnp} {output}")
=======
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.bam",
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.txt"
    params:
        threads = clusterconf["dedup"]["threads"],
        samplename = get_samplename,
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    shell:
        "shellbamfiles=$(echo {input.bamfiles} | sed 's/ / -i /g') ;"
        "{params.sentieon} driver -t {params.threads} "
            "-i $shellbamfiles "
            "--algo LocusCollector "
            "--fun score_info "
            "{wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt ;"
        "{params.sentieon} driver "
            "-t {params.threads} "
            "-i $shellbamfiles "
            "--algo Dedup "
            "--rmdup "
            "--score_info {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP_score.txt "
            "--metrics {wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.txt "
            "{wildcards.workingdir}/{wildcards.stype}/dedup/{wildcards.sname}_DEDUP.bam"

rule realign_mapping:
    input:
        "{workingdir}/{stype}/dedup/{sname}_DEDUP.bam"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    output:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    params:
        threads = clusterconf["realign_mapping"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"]
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} {output}"

rule baserecal:
    input:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    singularity:
        pipeconfig["singularities"]["sentieon"]["sing"]
    params:
        threads = clusterconf["baserecal"]["threads"],
        sentieon = pipeconfig["singularities"]["sentieon"]["tool_path"],
        referencegenome = pipeconfig["singularities"]["sentieon"]["reference"],
        dbsnp = pipeconfig["singularities"]["sentieon"]["dbsnp"],
        mills = pipeconfig["singularities"]["sentieon"]["mills"],
    output:
        "{workingdir}/{stype}/recal/{sname}_RECAL_DATA.TABLE"
    shell:
        "{params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.dbsnp} {output}"
>>>>>>> sentieon_singularity
