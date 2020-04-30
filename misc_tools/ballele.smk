# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule ballele_plot:
    input:
        "{stype}/dnascope/{sname}_germline.vcf"
    params:
        ballele = config["rules"]["ballele_plot"]["ballele"]
    output:
        "{stype}/reports/{sname}_baf.igv"
    run:
        shell("{params.ballele} -v {input} -o {output}")
