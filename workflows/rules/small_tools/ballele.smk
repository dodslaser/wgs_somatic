# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from workflows.scripts.b_allele_igv_plot.plot_b_allele_freq import plot_freq

rule ballele_plot:
    input:
        "{workingdir}/{stype}/dnascope/{sname}_germline.vcf"
    output:
        "{workingdir}/{stype}/reports/{sname}_baf.igv"
    run:
        plot_freq("{input}, {output}")
