# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule normal_mapping:
    input:
        fwd = join(normalfastqs, n_pattern_r1),
        rev = join(normalfastqs, n_pattern_r2),
    params:
        threads = clusterconf["normal_mapping"]["threads"],
        sentieon = sentieon,
        samplename = normalname,
        referencegenome = referencegenome,
        outputdir = config["rules"]["mapping"]["outputdir"]
    output:
        "normal/mapping/{normal_fastqpair}.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} bwa mem -M -R '@RG\\tID:{wildcards.normal_fastqpair}\\tSM:{params.samplename}\\tPL:ILLUMINA' -t 16 {params.referencegenome} {input.fwd} {input.rev} | {params.sentieon} util sort -o {output} -t 20 --sam2bam -i -")

rule tumor_mapping:
    input:
        fwd = join(tumorfastqs, t_pattern_r1),
        rev = join(tumorfastqs, t_pattern_r2),
    params:
        threads = clusterconf["tumor_mapping"]["threads"],
        sentieon = sentieon,
        samplename = tumorname,
        referencegenome = referencegenome,
        outputdir = config["rules"]["mapping"]["outputdir"]
    output:
        "tumor/mapping/{tumor_fastqpair}.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} bwa mem -M -R '@RG\\tID:{wildcards.tumor_fastqpair}\\tSM:{params.samplename}\\tPL:ILLUMINA' -t 16 {params.referencegenome} {input.fwd} {input.rev} | {params.sentieon} util sort -o {output} -t 20 --sam2bam -i -")

rule normal_dedup:
    input:
        bamfiles = expand("normal/mapping/{normal_fastqpair}.bam", normal_fastqpair=normal_fastqpairs)
    output:
        "normal/dedup/{normalname}_DEDUP.bam"
    params:
        threads = clusterconf["normal_dedup"]["threads"],
        normalname = normalname,
        sentieon = sentieon,
        outputdir = config["rules"]["mapping"]["outputdir"]
    run:
        inp_bamfiles = ""
        for bamfile in input.bamfiles:
            inp_bamfiles = f"{inp_bamfiles}-i {bamfile} "
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} {inp_bamfiles}--algo LocusCollector --fun score_info normal/dedup/{params.normalname}_DEDUP_score.txt")
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} {inp_bamfiles}--algo Dedup --rmdup --score_info normal/dedup/{params.normalname}_DEDUP_score.txt --metrics normal/dedup/{params.normalname}_DEDUP.txt normal/dedup/{params.normalname}_DEDUP.bam")

rule tumor_dedup:
    input:
        bamfiles = expand("tumor/mapping/{tumor_fastqpair}.bam", tumor_fastqpair=tumor_fastqpairs)
    output:
        "tumor/dedup/{tumorname}_DEDUP.bam"
    params:
        threads = clusterconf["tumor_dedup"]["threads"],
        tumorname = tumorname,
        sentieon = sentieon,
        outputdir = config["rules"]["mapping"]["outputdir"]
    run:
        inp_bamfiles = ""
        for bamfile in input.bamfiles:
            inp_bamfiles = f"{inp_bamfiles}-i {bamfile} "
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} {inp_bamfiles}--algo LocusCollector --fun score_info tumor/dedup/{params.tumorname}_DEDUP_score.txt")
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} {inp_bamfiles}--algo Dedup --rmdup --score_info tumor/dedup/{params.tumorname}_DEDUP_score.txt --metrics tumor/dedup/{params.tumorname}_DEDUP.txt tumor/dedup/{params.tumorname}_DEDUP.bam")

rule realign_normal_mapping:
    input:
        expand("normal/dedup/{normalname}_DEDUP.bam", normalname=normalname),
    params:
        threads = clusterconf["realign_normal_mapping"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        mills = config["rules"]["realign"]["mills"],
        tgenomes = config["rules"]["realign"]["tgenomes"],
        outputdir = config["rules"]["realign"]["outputdir"]
    output:
        "normal/realign/{normalname}_REALIGNED.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} -k {params.tgenomes} {output}")

rule realign_tumor_mapping:
    input:
        expand("tumor/dedup/{tumorname}_DEDUP.bam", tumorname=tumorname),
    params:
        threads = clusterconf["realign_normal_mapping"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        mills = config["rules"]["realign"]["mills"],
        tgenomes = config["rules"]["realign"]["tgenomes"],
        outputdir = config["rules"]["realign"]["outputdir"]
    output:
        "tumor/realign/{tumorname}_REALIGNED.bam"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo Realigner -k {params.mills} -k {params.tgenomes} {output}")

rule normal_baserecal:
    input:
        expand("normal/realign/{normalname}_REALIGNED.bam", normalname=normalname),
    params:
        threads = clusterconf["normal_baserecal"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = config["rules"]["recal"]["dbsnp"],
        mills = config["rules"]["recal"]["mills"],
        tgenomes = config["rules"]["recal"]["tgenomes"],
        outputdir = config["rules"]["recal"]["outputdir"]
    output:
        "normal/recal/{normalname}_RECAL_DATA.TABLE"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.tgenomes} -k {params.dbsnp} {output}")

rule tumor_baserecal:
    input:
        expand("tumor/realign/{tumorname}_REALIGNED.bam", tumorname=tumorname),
    params:
        threads = clusterconf["tumor_baserecal"]["threads"],
        sentieon = sentieon,
        referencegenome = referencegenome,
        dbsnp = config["rules"]["recal"]["dbsnp"],
        mills = config["rules"]["recal"]["mills"],
        tgenomes = config["rules"]["recal"]["tgenomes"],
        outputdir = config["rules"]["recal"]["outputdir"]
    output:
        "tumor/recal/{tumorname}_RECAL_DATA.TABLE"
    run:
        shell("export SENTIEON_LICENSE=medair1.medair.lcl:8990 ; {params.sentieon} driver -t {params.threads} -r {params.referencegenome} -i {input} --algo QualCal -k {params.mills} -k {params.tgenomes} -k {params.dbsnp} {output}")

