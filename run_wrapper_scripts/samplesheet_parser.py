#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python

import re
import csv
import json

from sample_sheet import SampleSheet

from run_wrapper_scripts.run_helpers import logger


# TODO MOVE TO CONFIG
def fetch_investigators():
    with open('run_wrapper_scripts/configs/investigator_list.json', 'r') as inp:
        return json.load(inp)['investigators']


# TODO REMOVE ME WHEN SLIMS
def translate_description(descr, run_id):
    case_dict = {"t": "tumor", "n": "normal"}
    invest, gender, pipeline, casetype, caseid = descr.split('_')
    if caseid:
        try:
            casetype = case_dict[casetype.lower()]
        except KeyError:
            logger(f"casetype {casetype} not tumor or normal, setting unknown", run_id)
            casetype = 'unknown'
    return invest, gender, pipeline, casetype, caseid


def fetch_sample_info(run_id, sample):
    sample_id = sample.Sample_ID
    date, _, _, chip, *_ = run_id.split('_')
    sample_id_extra = '_'.join([sample_id, date, chip])

    sample_project = sample.Sample_Project
    if sample_project == "IGNORE":
        sample_project = "TEST"

    research_project = re.search("[GB][2-3][0-9]-[0-9]{3}", sample_project)
    if research_project:
        research_project  = "yes"
    else:
        research_project = "no"

    description = sample.Description

    invest, gender, pipeline, casetype, caseid = translate_description(description, run_id)

    investigators = fetch_investigators()
    sample_investigator_info = investigators[invest]

    master = {
        'samplename': sample_id,
        'project': sample_project,
        'sampleid': sample_id_extra,
        'runid': run_id,
        'description': {
            'ID': description,
            'investigators': {
                'initials': invest,
                **sample_investigator_info
            },
            'gender': gender,
            'analysis_pipeline': pipeline
        },
        'caseid': caseid,
        'casetype': casetype
    }

    return master


def get_sample_info(run_id, sample_sheet_path, dmx_dict):
    sheet = SampleSheet(sample_sheet_path)
    sample_dict = dmx_dict
    sample_dict['runid'] =  run_id
    sample_dict['ssheet'] = sample_sheet_path
    sample_dict[run_id] = {'investigators': {}}
    
    investigators = fetch_investigators()
    for lab_inv in sheet.Header['Investigator Name'].split(','):
        sample_dict[run_id]['investigators'].update({lab_inv: investigators[lab_inv]})

    for sample in sheet.samples:
        sample_dict["samples"][sample.Sample_ID]["sample_info"] = fetch_sample_info(run_id, sample)


    logger(f"checking to see if samplesheet: {sample_sheet_path} contains all samples in run", run_id)
    complete = True
    for sample in dmx_dict["samples"]:
        if sample not in sample_dict["samples"]:
            logger(f"SampleSheet incomplete, {sample} missing from SampleSheet", run_id)
            complete = False

    return complete, sample_dict

