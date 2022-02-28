"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import os
import re
import glob
import yaml
from datetime import datetime
import json
from itertools import chain
import traceback

from definitions import CONFIG_PATH, ROOT_DIR, ROOT_LOGGING_PATH
from context import RunContext, SampleContext
from helpers import setup_logger
from tools.slims import get_sample_slims_info, SlimsSample, find_more_fastqs, get_pair_dict

logger = setup_logger('wrapper', os.path.join(ROOT_LOGGING_PATH, 'WS_wrapper.log'))


# Store info about samples to use for sending report emails
sample_status = {'missing_slims': [],
                 'unset_WS': [],
                 'approved': []}

def look_for_runs(root_path):
    '''Look for runs in demultiplexdir'''
    found_paths = glob.glob(os.path.join(root_path, '*'))
    regex = '^[0-9]{6}_A00687_[0-9]{4}_.{10}$'
    return [path for path in found_paths if re.search(regex, os.path.basename(path))]


def generate_context_objects(Rctx):
    '''Create Rctx and Sctx for a demultiplexed run'''

    # Read demultiplex stats file for sample names, fastq paths, and nr reads
    with open(Rctx.demultiplex_summary_path, 'r') as inp:
        demuxer_info = json.load(inp)


    for sample_id, sample_info in demuxer_info['samples'].items():
        logger.info(f'Setting up context for {sample_id}.')

        # Setup Sample context class and add listed fastq paths
        Sctx = SampleContext(sample_id)
        Sctx.add_fastq(sample_info['fastq_paths'])

        # Query Slims for clinical information and add to sample context
        logger.info(f'Fetching SLIMS info.')
        Sctx.slims_info = get_sample_slims_info(Sctx, run_tag = Rctx.run_tag)

        if not Sctx.slims_info:
            logger.warning(f'No SLIMS info available!')
            logger.warning(f'Sample will not be analysed.')
            sample_status['missing_slims'].append(Sctx)
            continue

        # NOTE: 54 is Slims internal primary key for wgs_somatic
        if 54 not in Sctx.slims_info['secondary_analysis']:
            sample_status['unset_WS'].append(Sctx)
            continue

        # Add sample context to run context list
        logger.info(f'Adding sample context to list of contexts.')
        Rctx.add_sample_context(Sctx)

    if not Rctx.sample_contexts:
        # If no samples set for wgs_somatic
        # doesn't skip continuing with the rest of the code - need to fix this
        logger.info('No samples set for wgs_somatic. Skipping run.')

    for Sctx in Rctx.sample_contexts:
        sample_status['approved'].append(Sctx)

    return Rctx

def wrapper():

    # Empty dict, will update later with T/N pair info
    pair_dict_all_pairs = {}

    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)

    # Grab all available local run paths
    instrument_root_path = config['novaseq']['demultiplex_path']
    local_run_paths = look_for_runs(instrument_root_path)

    # Read all previously analysed runs
    previous_runs_file = config['previous_runs_file_path']
    previous_runs_file_path = os.path.join(ROOT_DIR, previous_runs_file)
    with open(previous_runs_file_path, 'r') as prev:
        previous_runs = [line.rstrip() for line in prev]

    # Loop through each run path and setup Run context class
    for run_path in local_run_paths:
        Rctx = RunContext(run_path)
        # Check if demultiplexing is completed
        if not Rctx.demultiplex_complete:
            continue

        # Check if run has been previously analysed
        if Rctx.run_name in previous_runs:
            continue

        # Register start time
        start_time = datetime.now()
        logger.info(f'Started {Rctx.run_name} at {start_time}')

        # Write run name to previously analysed list to ensure no double-running
        with open(previous_runs_file_path, 'a') as prev:
            logger.info(f'Writing {Rctx.run_name} to previous runs list.')
            print(Rctx.run_name, file=prev)

        # get Rctx and Sctx for current run
        Rctx_run = generate_context_objects(Rctx)

        # Get T/N pair info in a dict for samples and link additional fastqs from other runs
        for sctx in Rctx_run.sample_contexts:
            pair_dict = get_pair_dict(sctx, Rctx.run_tag, logger)
            pair_dict_all_pairs.update(pair_dict)

    # Uses the dictionary of T/N samples to put the correct pairs together and finds the correct input arguments to the pipeline
    for key in pair_dict_all_pairs:
        if 'tumor' in pair_dict_all_pairs.get(key):
            t = key
            t_ID = [val for val in pair_dict_all_pairs.get(key) if val != 'tumor'][0] # crappy solution... but works since there is only one item in the list
            for k in pair_dict_all_pairs:
                if 'normal' in pair_dict_all_pairs.get(k):
                    n = k
                    n_ID = [val for val in pair_dict_all_pairs.get(k) if val != 'normal'][0]
                    if n_ID == t_ID or t_ID == n.split("DNA")[1] or n_ID == t.split("DNA")[1]: # if we change to pair id instead of tumorNormalID this is needed
                        runnormal = Rctx_run.run_name
                        runtumor = Rctx_run.run_name
                        tumorsample = t
                        normalsample = n
                        normalfastqs = os.path.join(Rctx_run.run_path, "fastq")
                        tumorfastqs = os.path.join(Rctx_run.run_path, "fastq")
                        #outputdir = os.path.join("/seqstore/webfolders/wgs/barncancer/hg38", t) 
                        outputdir = os.path.join("/home/xshang/ws_testoutput/outdir/", t) #CHANGE BACK TO CORRECT OUTDIR
                        igvuser = 'barncancer_hg38' #FIXME get from config instead
                        hg38ref = 'yes' #FIXME get from config instead

                        if os.path.exists(outputdir):
                            logger.info(f'Outputdir exists for {tumorsample}. Renaming old outputdir {outputdir} to {outputdir}_old')
                            os.rename(outputdir, f'{outputdir}_old')

                        logger.info(f'Starting wgs_somatic with arguments: \n \
runnormal: {runnormal} \n \
runtumor: {runtumor} \n \
tumorsample: {tumorsample} \n \
normalsample: {normalsample} \n \
normalfastqs: {normalfastqs} \n \
tumorfastqs: {tumorfastqs} \n \
outputdir: {outputdir} \n \
igvuser: {igvuser} \n \
hg38ref: {hg38ref}')

    # start the pipeline with the correct pairs. 
    # will use these arguments to start pipeline. 
    # some arguments are hardcoded right now, need to fix this. 
    # only considers barncancer hg38 (GMS-AL + GMS-BT samples) right now. 
    # could have outputdirs and arguments in config and get them from there. 

    # this will only work for pairs. have to consider tumor only/normal only as well. 

    # the arguments runtumor and runnormal could be "wrong" by doing it like this 
    # since they use run name of current run 
    # but if for example normal is in older run it has the wrong argument for "runnormal". 
    # this argument is not that important, it is only used to create a unique sample name. 
    # maybe it would be better discard/modify this argument than spending time on getting the correct value of it for samples from different runs. 
    # also, if fastqs come from more than one run, what will the value of this argument be then to be "correct"?...

    # outputdir - need to consider if outputdir already exists 
    # (if sample has been run before and now it has new fastqs in current run, outputdir already exists). 
    # should old outputdir be moved to archive? 



   # next step here is to actually start the pipeline with these arguments



if __name__ == '__main__':
    try:
        wrapper()
    except KeyboardInterrupt:
        pass
    except Exception:
        format_exc = traceback.format_exc()
        logger.error(format_exc)
        # TODO: add send email about error here

