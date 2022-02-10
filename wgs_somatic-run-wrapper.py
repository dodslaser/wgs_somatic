"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import os
import re
import glob
import yaml
from datetime import datetime
import json

from definitions import CONFIG_PATH, ROOT_DIR, ROOT_LOGGING_PATH
from context import RunContext, SampleContext
from helpers import setup_logger
from tools.slims import get_sample_slims_info, SlimsSample, more_fastqs, get_pair

logger = setup_logger('wrapper', os.path.join(ROOT_LOGGING_PATH, 'WS_wrapper.log'))


def look_for_runs(root_path):
# Look for runs in demultiplexdir
    found_paths = glob.glob(os.path.join(root_path, '*'))
    regex = '^[0-9]{6}_A00687_[0-9]{4}_.{10}$'
    return [path for path in found_paths if re.search(regex, os.path.basename(path))]




def wrapper():

    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)

    # Grab all available local run paths
    instrument_root_path = config['novaseq']['demultiplex_path']
    local_run_paths = look_for_runs(instrument_root_path)
    #print(local_run_paths)

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

        # Read demultiplex stats file for sample names, fastq paths, and nr reads
        with open(Rctx.demultiplex_summary_path, 'r') as inp:
            demuxer_info = json.load(inp)

        # Store info about samples to use for sending report emails
        sample_status = {'missing_slims': [],
                         'unset_WS': [],
                         'approved': []}

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
            print(Sctx.slims_info)




            # NOTE: 54 is Slims internal primary key for wgs_somatic
            if 54 not in Sctx.slims_info['secondary_analysis']:
                sample_status['unset_WS'].append(Sctx)
                continue


            # Add sample context to run context list
            logger.info(f'Adding sample context to list of contexts.')
            Rctx.add_sample_context(Sctx)
            #print(Sctx)


            #more_fqs_path = more_fastqs(Sctx, run_tag = Rctx.run_tag)
            #if more_fqs_path:
            #    print(more_fqs_path)
            #    Sctx.add_fastq(more_fqs_path)
            #print(get_pair(Sctx, Rctx.run_tag))
            more_fqs_dicts = get_pair(Sctx, Rctx.run_tag)
            if more_fqs_dicts:
                print(f'more_fqs_dicts: {more_fqs_dicts}')
                for d in more_fqs_dicts:
                    print(f'value of key Sample ID: {d["Sample ID"]}')                
                    print(Sctx.sample_name)
                    if Sctx.sample_name == d["Sample ID"]:
                        Sctx.add_fastq(d["fastq paths"])
                    print(f'value of fq paths: {d["fastq paths"]}')
#                    for fqp in d["fastq paths"]:
#                        print(fqp)
            print(f'all fastq paths: {Sctx.fastqs}')
                #Sctx.add_fastq(more_fqs_path)


#        print(Rctx.sample_contexts)

        if not Rctx.sample_contexts:
            # If no samples set for wgs_somatic
            logger.info('No samples set for wgs_somatic. Skipping run.')
            continue

        for Sctx in Rctx.sample_contexts:
            sample_status['approved'].append(Sctx)
        runID = Rctx.run_name
    print(sample_status)
    tumor_samples = []
    normal_samples = []
    for category, contexts in sample_status.items():
        print(category) 
        for Sctx in contexts:
            # if two samples have the same tumorNormalID and one is tumor and one is normal - start pipeline
            print(Sctx.sample_id, Sctx.slims_info['tumorNormalType'], Sctx.slims_info['tumorNormalID'])
        if category == 'approved':
            for Sctx in contexts:
                if Sctx.slims_info['tumorNormalType'] == 'tumor':
                    tumor_samples.append(Sctx)
                elif Sctx.slims_info['tumorNormalType'] == 'normal':
                    normal_samples.append(Sctx)
                else:
                    logger.info(f'Warning! {Sctx.slims_info["content_id"]} is not set as tumor or normal.')
    #print('tumor')
    #print(tumor_samples)
    #print('normal')
    #print(normal_samples)


    # find tumor/normal pairs in the run and start pipeline
    for t in tumor_samples:
        print(t.fastqs)
        for n in normal_samples:
            if t.slims_info['tumorNormalID'] == n.slims_info['tumorNormalID']:
                logger.info(f'Starting wgs_somatic for: \n Run: {runID} \n Tumor: {t.slims_info["content_id"]} \n Normal: {n.slims_info["content_id"]} ')
    


if __name__ == '__main__':
    try:
        wrapper()
    except KeyboardInterrupt:
        pass
    except Exception:
        pass


