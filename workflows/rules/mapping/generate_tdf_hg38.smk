# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule generate_tdf_hg38:
    input:
        "{workingdir}/{stype}/realign/{sname}_REALIGNED.bam"
    params:
        #igvjar = pipeconfig["rules"]["generate_tdf"]["tdfgen"]
        igv = pipeconfig["rules"]["generate_tdf"]["tdfgen"]
    output:
        "{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf"
    run:
        shell("nohup java -Xmx32768m --module-path={params.igv}/lib @{params.igv}/igv.args --module=org.igv/org.broad.igv.tools.IgvTools count {input} {output} hg38")
