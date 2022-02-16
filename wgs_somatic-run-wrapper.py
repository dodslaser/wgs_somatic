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
    '''Create Sctx for a demultiplexed run'''

    # Read demultiplex stats file for sample names, fastq paths, and nr reads
    with open(Rctx.demultiplex_summary_path, 'r') as inp:
        demuxer_info = json.load(inp)


    for sample_id, sample_info in demuxer_info['samples'].items():
        logger.info(f'Setting up context for {sample_id}.')
        #print(sample_id, sample_info)

        # Setup Sample context class and add listed fastq paths and nr reads
        Sctx = SampleContext(sample_id)
        Sctx.add_fastq(sample_info['fastq_paths'])
        #print(Sctx.fastqs)

        # Query Slims for clinical information and add to sample context
        logger.info(f'Fetching SLIMS info.')
        get_sample_slims_info(Sctx, run_tag = Rctx.run_tag)  # NOTE: Mod with slims info in-place

        if not Sctx.slims_info:
            logger.warning(f'No SLIMS info available!')
            logger.warning(f'Sample will not be analysed.')
            sample_status['missing_slims'].append(Sctx)
            continue
        #print(f'slims info: {Sctx.slims_info}')

        # NOTE: 54 is Slims internal primary key for wgs_somatic
        if 54 not in Sctx.slims_info['secondary_analysis']:
            sample_status['unset_WS'].append(Sctx)
            continue
        # Add sample context to run context list
        logger.info(f'Adding sample context to list of contexts.')
        Rctx.add_sample_context(Sctx)
        #print(f'Sctx: {Sctx}')
        #Sctx_list.append(Sctx)

        if not Rctx.sample_contexts:
            # If no samples set for wgs_somatic
            logger.info('No samples set for wgs_somatic. Skipping run.')
            continue

    for Sctx in Rctx.sample_contexts:
    #    print(f'Sctx: {Sctx}')
        sample_status['approved'].append(Sctx)
    #runID = Rctx.run_name
    #print(f'sample status: {sample_status}')
    #tumor_samples = []
    #normal_samples = []

    return Rctx


def link_fastqs(list_of_fq_paths):
    '''Link fastqs to fastq-folder in demultiplexdir of current run. Need to change the hardcoded path to my home... '''
    for fq_path in list_of_fq_paths:
        print(f'fq_path: {fq_path}')
    # Only links if link doesn't already exist
        if not os.path.islink(os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(fq_path))):
        # Now symlinks all additional paths to fastqs for tumor and normal in other runs. If I symlink to demultiplexdir of particular run instead, all fastqs belonging to the T/N pair will be in the same folder and the pipeline can start using that folder as argument.
            os.symlink(fq_path, os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(fq_path)))


def wrapper():


    additional_run_paths = []
    tumor_samples = []
    normal_samples = []
    pair_ids_in_run = []
    samples_ready = []

    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)

    # Grab all available local run paths
    instrument_root_path = config['novaseq']['demultiplex_path']
    local_run_paths = look_for_runs(instrument_root_path)
    #print(f'local_run_paths: {local_run_paths}')

    # Read all previously analysed runs
    previous_runs_file = config['previous_runs_file_path']
    previous_runs_file_path = os.path.join(ROOT_DIR, previous_runs_file)
    with open(previous_runs_file_path, 'r') as prev:
        previous_runs = [line.rstrip() for line in prev]
    #print(previous_runs)



    # Loop through each run path and setup Run context class
    for run_path in local_run_paths:
        Rctx = RunContext(run_path)
        #print(Rctx.run_path, Rctx.run_name, Rctx.run_date, Rctx.run_flowcell, Rctx.run_tag)
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

        Rctx_run = Rctx_Sctx_for_run(Rctx)
        #Sctx_list = Rctx_run.sample_contexts
        #print(f'Sctx_list: {Sctx_list}')

        #print(f'Rctx run sample contexts: {Rctx_run.sample_contexts}')


        # Get run paths for samples (or other part of t/n pair) with additional fastqs in other runs
        for sctx in Rctx_run.sample_contexts:
            #print(f'sctx: {sctx}')
            run_paths = get_pair_and_more_fqs(sctx, Rctx.run_tag)
            #print(f'run paths: {run_paths}')
            if not run_paths == None:
                run_paths = list(chain.from_iterable(run_paths))
                for r in run_paths:
                    #print(f'r: {r}')
                    additional_run_paths.append(r) if r not in additional_run_paths else additional_run_paths

    for run_path in additional_run_paths:
        additionalRctx = RunContext(run_path)
        additionalRctx = Rctx_Sctx_for_run(additionalRctx)
        #print(f'additional Rctx: {additionalRctx}')
        
    #print(f'sample_status: {sample_status}')

    for category, contexts in sample_status.items():
        print(category) 
        for Sctx in contexts:
            # if two samples have the same tumorNormalID and one is tumor and one is normal - start pipeline
            print(Sctx.sample_id, Sctx.slims_info['tumorNormalType'], Sctx.slims_info['tumorNormalID'])
        if category == 'approved':
            #print(f'contexts: {contexts}')
            for Sctx in contexts:
                print(f'Rctx run tag: {Rctx_run.run_tag}')
                print(f'sample id split: {Sctx.sample_id.split("_",1)[1]}')
                print(f'pairid: {Sctx.slims_info["tumorNormalID"]}, samplenamesplit = pairid: {Sctx.sample_name.split("DNA")[1]}')
                if Rctx_run.run_tag == Sctx.sample_id.split("_",1)[1]:
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
    print(f'pair_ids_in_run: {pair_ids_in_run}')
    print('tumor')
    print(tumor_samples)
    print('normal')
    print(normal_samples)

    print(f'Rctx run tag: {Rctx_run.run_tag}')
    # find tumor/normal pairs in the run and start pipeline
    for t in tumor_samples:
        print(t.fastqs)
        #print(t.sample_name)
        print(f't sample name split = pairid: {t.sample_name.split("DNA")[1]}, t tumor id = pairid: {t.slims_info["tumorNormalID"]}')
        print(f'tumor sample id: {t.sample_id}')
        if not any(pair_id in pair_ids_in_run for pair_id in (t.sample_name.split("DNA")[1], t.slims_info["tumorNormalID"])):
        #if not t.sample_name.split("DNA")[1] in pair_ids_in_run or not t.slims_info["tumorNormalID"] in pair_ids_in_run:
            print(f'{t.sample_name} does not belong to current run')
            continue
        if t.sample_id.split('_',1)[1] == Rctx_run.run_tag:
            # the plan is to link the fastqs of other runs to fastq-folder in demultiplexdir for the current run so all fastqs for a sample + its pair are in the same folder. then the pipeline can start based on this folder.
            print('these fastqs do not need to be linked')
        else:
            # only need to link fastqs from other runs. 
            print('linking fastqs')
            link_fastqs(t.fastqs)
        samples_ready.append(t)

# same code as for tumor... could make a function instead
    for n in normal_samples:
        if not any(pair_id in pair_ids_in_run for pair_id in (n.sample_name.split("DNA")[1], n.slims_info["tumorNormalID"])):
        #if not n.sample_name.split("DNA")[1] in pair_ids_in_run or not n.slims_info["tumorNormalID"] in pair_ids_in_run:
            print(f'{n.sample_name} does not belong to current run')
            continue
        if n.sample_id.split('_',1)[1] == Rctx_run.run_tag:
            print('these fastqs do not need to be linked')
        else:
            print('linking fastqs')
            link_fastqs(n.fastqs)
        samples_ready.append(n)

    print(f'samples ready: {samples_ready}')
        #for n in normal_samples:
        #    if t.slims_info['tumorNormalID'] == n.slims_info['tumorNormalID']:
        #        logger.info(f'Starting wgs_somatic for: \n Run: {runID} \n Tumor: {t.slims_info["content_id"]} \n Normal: {n.slims_info["content_id"]} ')
    


if __name__ == '__main__':
    try:
        wrapper()
    except KeyboardInterrupt:
        pass
    except Exception:
        pass


