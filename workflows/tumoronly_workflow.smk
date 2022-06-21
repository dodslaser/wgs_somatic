# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import os
from collections import defaultdict

def group_insilico(paths):
    insilico_groups = defaultdict(list)

    for path in paths:
        insilico_name = os.path.basename(os.path.dirname(path))  # NOTE: Not great...
        insilico_groups[insilico_name].append(path)

    # Collect files under keys for better access than relying on indexes
    insilico_groups_named = {}
    for insilico_name, insilico_paths in insilico_groups.items():
        ten_x, twenty_x, below_ten_x, csv = insilico_paths

        named = {'10x': ten_x,
                 '20x': twenty_x,
                 'below_10x': below_ten_x,
                 'csv': csv}
        insilico_groups_named[insilico_name] = named

    return insilico_groups_named


def get_insilico(wcs):
    insilico_data = config['insilico']
    insilico_names = insilico_data.keys()

    # NOTE: sampleid and workingdir is essentially global because defined in snakefile
    insilico_files = []
    for insilico_name in insilico_names:
        insilico_files.extend([f"{workingdir}/tumor/insilico/{insilico_name}/{tumorid}_{insilico_name}_10x.xlsx",
                               f"{workingdir}/tumor/insilico/{insilico_name}/{tumorid}_{insilico_name}_20x.xlsx",
                               f"{workingdir}/tumor/insilico/{insilico_name}/{tumorid}_{insilico_name}_genes_below10x.xlsx",
                               f"{workingdir}/tumor/insilico/{insilico_name}/{tumorid}_{insilico_name}.csv"])
    return insilico_files


rule tumoronly_workflow:
    input:
        expand("{workingdir}/{stype}/tnscope/{sname}_somatic.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/dnascope/{sname}_germline.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_baf.igv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_WGScov.tsv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_Ycov.tsv", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_CNV_germline.vcf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_{vartype}_CNV_observed.seg", workingdir=workingdir, vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/canvas/{sname}_{vartype}_CNV_called.seg", workingdir=workingdir, vartype="germline", sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/reports/{sname}_REALIGNED.bam.tdf", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        expand("{workingdir}/{stype}/pindel/{sname}_pindel.xlsx", workingdir=workingdir, sname=tumorid, stype=sampleconfig[tumorname]["stype"]),
        "{workingdir}/reporting/shared_result_files.txt",
        insilico_files = get_insilico
    output:
        insilico_json = "{workingdir}/reporting/insilico.json",
        wf_finished = "{workingdir}/reporting/workflow_finished.txt"
    run:
        output_mapping = dict(input)
        output_mapping['insilico_files'] = group_insilico(output_mapping['insilico_files'])
        with open(output["insilico_json"], "w") as j:
            json.dump(output_mapping, j, indent=4)
        shell(f"echo {input} >> {output['wf_finished']}")
