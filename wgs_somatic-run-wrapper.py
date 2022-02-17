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

from definitions import CONFIG_PATH, ROOT_DIR, ROOT_LOGGING_PATH
from context import RunContext, SampleContext
from helpers import setup_logger
from tools.slims import get_sample_slims_info, SlimsSample, more_fastqs, get_pair_and_more_fqs

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


def Rctx_Sctx_for_run(Rctx):
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
        get_sample_slims_info(Sctx, run_tag = Rctx.run_tag)  # NOTE: Mod with slims info in-place

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
            logger.info('No samples set for wgs_somatic. Skipping run.')
            continue

    for Sctx in Rctx.sample_contexts:
        sample_status['approved'].append(Sctx)

    return Rctx


def link_fastqs(list_of_fq_paths):
    '''Link fastqs to fastq-folder in demultiplexdir of current run. Need to change the hardcoded path to my home... '''
    # using a hardcoded test folder right now for symlinks. will change this to correct Demultiplexdir/current-run/fastq folder.
    # additional fastqs need to still be in demultiplexdir. not considering downloading from hcp right now. need to consider this later...
    for fq_path in list_of_fq_paths:
    # Only links if link doesn't already exist
        if not os.path.islink(os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(fq_path))):
        # Now symlinks all additional paths to fastqs for tumor and normal in other runs. If I symlink to demultiplexdir of particular run instead, all fastqs belonging to the T/N pair will be in the same folder and the pipeline can start using that folder as argument.
            os.symlink(fq_path, os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(fq_path)))


def get_samples_ready(list_of_samples, pair_ids_in_run, run_tag):
    '''Add samples that belong to current run from current run and additional runs to list of ready samples'''
    samples_ready = []
    for s in list_of_samples:
        # will have Sctx for all samples set for wgs-somatic in other runs that have samples that are related to current run. samples that are not related to current run shouldn't run again in wgs-somatic.
        if not any(pair_id in pair_ids_in_run for pair_id in (s.sample_name.split("DNA")[1], s.slims_info["tumorNormalID"])):
            logger.info(f'{s.sample_name} does not belong to current run')
            continue
        if s.sample_id.split('_',1)[1] == run_tag:
            # the plan is to link the fastqs of other runs to fastq-folder in demultiplexdir for the current run so all fastqs for a sample + its pair are in the same folder. then the pipeline can start based on this folder.
            logger.info(f'fastqs for {s.sample_id} do not need to be linked')
        else:
            # only need to link fastqs from other runs since i will link them to fastq folder of current run.
            logger.info(f'linking fastqs for {s.sample_id}')
            link_fastqs(s.fastqs)
        samples_ready.append(s)
    return samples_ready


def wrapper():

    # using lists to keep track of stuff... could maybe be done in a better way...
    additional_run_paths = []
    tumor_samples = []
    normal_samples = []
    pair_ids_in_run = []
    started_samples = []


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
        Rctx_run = Rctx_Sctx_for_run(Rctx)


        # Get run paths for samples (or other part of t/n pair) with additional fastqs in other runs
        for sctx in Rctx_run.sample_contexts:
            run_paths = get_pair_and_more_fqs(sctx, Rctx.run_tag)
            if not run_paths == None:
                run_paths = list(chain.from_iterable(run_paths))
                for r in run_paths:
                    additional_run_paths.append(r) if r not in additional_run_paths else additional_run_paths

    # get Rctx and Sctx for additional runs that have samples related to current run
    # i realized that i don't actually use "additionalRctx" for anything now. but using the function Rctx_Sctx_for_run does append to sample_status which is used below... could be done in a better way if i don't have to use additionalRctx anyway. just appending to sample_status is neccessary.  
    for run_path in additional_run_paths:
        additionalRctx = RunContext(run_path)
        additionalRctx = Rctx_Sctx_for_run(additionalRctx)


    # get tumor and normal samples related to current run from approved samples
    for category, contexts in sample_status.items():
        if category == 'approved':
            for Sctx in contexts:
                if Rctx_run.run_tag == Sctx.sample_id.split("_",1)[1]:
                    # use both tumorNormalID and sample name (minus "DNA") as pair ids since we will change pair ids to be normalname for tumor and tumorname for normal 
                    pair_ids_in_run.append(Sctx.slims_info["tumorNormalID"])
                    pair_ids_in_run.append(Sctx.sample_name.split("DNA")[1])
                if Sctx.slims_info['tumorNormalType'] == 'tumor':
                    tumor_samples.append(Sctx)
                elif Sctx.slims_info['tumorNormalType'] == 'normal':
                    normal_samples.append(Sctx)
                else:
                    logger.info(f'Warning! {Sctx.slims_info["content_id"]} is not set as tumor or normal.')
    
    # make list of unique pair ids in current run
    pair_ids_in_run = list(set(pair_ids_in_run))
    logger.info(f'pair IDs in run: {pair_ids_in_run}')



    # find tumor/normal pairs in the run and start pipeline

    # get lists of ready samples and remove duplicates from list
    tumor_samples_ready = list(set(get_samples_ready(tumor_samples, pair_ids_in_run, Rctx_run.run_tag)))
    normal_samples_ready = list(set(get_samples_ready(normal_samples, pair_ids_in_run, Rctx_run.run_tag)))


    # start the pipeline with the correct pairs. will use these arguments to start pipeline. some arguments are hardcoded right now, need to fix this. only considers barncancer hg38 (GMS-AL + GMS-BT samples) right now. could have outputdirs and arguments in config and get them from there. 
    # this will only work for pairs. have to consider tumor only/normal only as well. 
    # the arguments runtumor and runnormal could be "wrong" by doing it like this since they use run name of current run but if for example normal is in older run it has the wrong argument for "runnormal". this argument is not that important, it is only used to create a unique sample name. maybe it would be better discard/modify this argument than spending time on getting the correct value of it for samples from different runs. also, if fastqs come from more than one run, what will the value of this argument be then to be "correct"?...
    # outputdir - need to consider if outputdir already exists (if sample has been run before and now it has new fastqs in current run, outputdir already exists). should old outputdir be moved to archive? 
    for t in tumor_samples_ready:
        if t.slims_info["content_id"] in started_samples:
            continue
        for n in normal_samples_ready:
            if n.slims_info["content_id"] in started_samples:
                continue
            if t.slims_info['tumorNormalID'] == n.slims_info['tumorNormalID']:
                logger.info(f'Starting wgs_somatic with arguments: \n runnormal: {Rctx_run.run_name} \n runtumor: {Rctx_run.run_name} \n tumorsample: {t.slims_info["content_id"]} \n normalsample: {n.slims_info["content_id"]} \n normalfastqs: {os.path.join(Rctx_run.run_path, "fastq")} \n tumorfastqs: {os.path.join(Rctx_run.run_path, "fastq")} \n outputdir: {os.path.join("/seqstore/webfolders/wgs/barncancer/hg38", t.slims_info["content_id"])} \n igvuser: barncancer_hg38 \n hg38ref: yes')
                started_samples.extend((t.slims_info["content_id"], n.slims_info["content_id"]))


   # next step here is to actually start the pipeline with these arguments



if __name__ == '__main__':
    try:
        wrapper()
    except KeyboardInterrupt:
        pass
    except Exception:
        pass


