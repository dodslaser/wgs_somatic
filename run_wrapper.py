#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import os
import sys
import json
import glob
import ntpath
import shutil
import argparse

from run_wrapper_scripts.run_helpers import send_mail
from run_wrapper_scripts.run_helpers import read_config
from run_wrapper_scripts.run_helpers import read_investigators
from run_wrapper_scripts.run_helpers import logger

from run_wrapper_scripts.demultiplex_stats import fetch_stats
from run_wrapper_scripts.samplesheet_parser import get_sample_info

from launch_snakemake import analysis_main

def run_start_wrapper(run_path):
    newrun = ntpath.basename(run_path)
    try:
        run_start(run_path)
    except Exception as e:
        tb = traceback.format_exc()
        logger(f"Error in run-level script: {e}. Traceback: {tb}", newrun)  # TODO Real logging
        invest_config = read_investigators()
        adminpersonel = invest_config["investigators"]["ADMIN"]
        email = adminpersonel["email"]
        attachment = []
        mailsubject = f"Errors in {newrun}"
        mailbody = f"Error occurred in runlevel of script: {e}. Traceback: {tb}"
        send_mail(mailsubject, mailbody, attachment, email_list)
        return

# RUN START WRAPPER FUNCTION
def run_start(dmx_path):
    if dmx_path.endswith("/"):
        dmx_path = dmx_path[:-1]
    newrun = ntpath.basename(dmx_path)

    config = read_config()
    demux_json_primary = f"{dmx_path}/{config['dmx_stats_files']['primary']}"
    demux_json_secondary = f"{dmx_path}/{config['dmx_stats_files']['secondary']}"
    if os.path.isfile(demux_json_primary):
        demux_json = demux_json_primary
        with open(demux_json, "r") as demux_json_read:
            dmx_dict = json.load(demux_json_read) 
    elif os.path.isfile(demux_json_secondary):
        demux_json = demux_json_secondary
        dmx_dict = fetch_stats(demux_json, dmx_path)
        # create demuxer.json 
        with open(demux_json_primary, 'w') as demuxerjson:
            json.dump(dmx_dict, demuxerjson, ensure_ascii=False, indent=4)
    else:
        logger(f"No dmxstats found in {newrun}", newrun)
        sys.exit()    
    logger(f"New dmxstats discovered detected pertaining to {newrun} ({demux_json}), checking to see if it contains samples to analyse in method", newrun)

    samplesheet_primary = config["samplesheet_locs"]["primary"].replace("[runid]", f"{newrun}")
    samplesheet_secondary = config["samplesheet_locs"]["secondary"].replace("[runid]", f"{newrun}")
    primary_samplesheets = glob.glob(f"{samplesheet_primary}/SampleSheet*")
    secondary_samplesheets = glob.glob(f"{samplesheet_secondary}/SampleSheet*")
    samplesheets = primary_samplesheets + secondary_samplesheets
    if samplesheets:
        complete = False
        for samplesheet in samplesheets:
            complete, sampledict = get_sample_info(newrun, samplesheet, dmx_dict)
            if complete:
                logger(f"Samplesheet found ({samplesheet}) and info added to sampledict. Dumping sampledict to loglocation", newrun)
                logger(sampledict, newrun)
                if samplesheet in secondary_samplesheets:
                    shutil.copyfile(samplesheet, f"{samplesheet_primary}/{os.path.basename(samplesheet)}")
                break
        if not complete:
            logger(f"Samplesheets found: {samplesheets} but none contains all samples in dmx", newrun)
            logger(sampledict, newrun, "incomplete_samplesheet")
    else:
        logger(f"No SampleSheets found for {newrun} in either primary or secondary search_location, {samplesheet_primary}, {samplesheet_secondary}", newrun) 
        sys.exit()

#    analysis_dict = 
#    sample_sheet_path = f"{dmx_path}/SampleSheet.csv"
    #if not sample_sheet_path

    #logger("Reading SampleSheet.csv...", newrun)
    #sampledict = get_sample_info(newrun, sample_sheet_path, dmx_dict)

    #logger("Sampledict defined, dumping sampledict of run to run loglocation.", newrun)
    #logger(sampledict, newrun)

    #sampledict = demultiplex_stats_file(sampledict)

    #print(f"For admin -- merging data between runs not implemented, need to create demuxer.jsons for old runs and copy to {demuxerdir} then change parseing in look_for_merge_data function")

    ###---> Organize Samples into Departmentdicts and exlucde non-relevant samples <---###
    #investigatordict = prepare_dict(sampledict)
    #logger(investigatordict, newrun)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dmxdir', nargs='?', help='path to dmxdir', required=True)
    args = parser.parse_args()
    run_start(args.dmxdir)
