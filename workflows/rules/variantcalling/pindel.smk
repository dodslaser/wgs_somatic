# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule pindelConfig:
    input:
        tumorbam = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.bam", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        normalbam = expand("{workingdir}/{stype}/dedup/{sname}_DEDUP.bam", workingdir=workingdir, sname=normalid, stype=sampleconfig[normalname]["stype"]),
    output:
        pindelConfig = "{workingdir}/{stype}/pindel/{sname}_pindelConfig.txt"
    shell:
        "echo '{input.tumorbam}\t300\t{tumorname}\n{input.normalbam}\t300\t{normalname}'>{output.pindelConfig}"

rule pindel:
    input:
        pindelConfig = expand("{workingdir}/{stype}/pindel/{sname}_pindelConfig.txt", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"])
    params:
        bed = pipeconfig["rules"]["pindel"]["bed"],
        reference = pipeconfig["referencegenome"],
        threads = clusterconf["pindel"]["threads"],
        x = 2,
        B = 60
    singularity:
        pipeconfig["singularities"]["pindel"]["sing"]
    output:
        "{workingdir}/{stype}/pindel/{sname}_BP",
        "{workingdir}/{stype}/pindel/{sname}_CloseEndMapped",
        "{workingdir}/{stype}/pindel/{sname}_D",
        "{workingdir}/{stype}/pindel/{sname}_INT_final",
        "{workingdir}/{stype}/pindel/{sname}_INV",
        "{workingdir}/{stype}/pindel/{sname}_LI",
        "{workingdir}/{stype}/pindel/{sname}_RP",
        "{workingdir}/{stype}/pindel/{sname}_SI",
        "{workingdir}/{stype}/pindel/{sname}_TD"
    shell:
        "echo $HOSTNAME;"
        " (pindel -f {params.reference} -i {input.pindelConfig} -T {params.threads} -x {params.x} -B {params.B} -j {params.bed} -o {workingdir}/{wildcards.stype}/pindel/{wildcards.sname} ) "

rule pindel2vcf:
    input:
        expand("{workingdir}/{stype}/pindel/{sname}_BP", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_CloseEndMapped", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_D", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_INT_final", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_INV", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_LI", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_RP", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_SI", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_TD", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"])
    params:
        threads = clusterconf["pindel"]["threads"],
        reference = pipeconfig["referencegenome"],
        refname = "GRCh38",
        refdate = 000000,
        e = 3, #e = 10,
        mc = 10,
        minsize = 5
    singularity:
        pipeconfig["singularities"]["pindel"]["sing"]
    output:
        "{workingdir}/{stype}/pindel/{sname}_pindel.vcf"
    shell:
        "echo $HOSTNAME;"
        " (pindel2vcf -P {workingdir}/{wildcards.stype}/pindel/{wildcards.sname} -r {params.reference} -R {params.refname} -d {params.refdate} -v {output} -e {params.e} -mc {params.mc} -G -is {params.minsize}) "
