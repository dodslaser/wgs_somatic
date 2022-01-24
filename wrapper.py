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


logger = setup_logger('wrapper', os.path.join(ROOT_LOGGING_PATH, 'WS_wrapper.log'))
#logger.info(f'hej')


def look_for_runs(root_path):
    found_paths = glob.glob(os.path.join(root_path, '*'))
    regex = '^[0-9]{6}_A00687_[0-9]{4}_.{10}$'
    return [path for path in found_paths if re.search(regex, os.path.basename(path))]


# get root path from a config
#root_path = "/seqstore/instruments/novaseq_687_gc/Demultiplexdir/"
#print(look_for_runs(root_path))


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

        for sample_id, sample_info in demuxer_info['samples'].items():
            logger.info(f'Setting up context for {sample_id}.')
            #print(sample_id, sample_info)

            # Setup Sample context class and add listed fastq paths and nr reads
            Sctx = SampleContext(sample_id)
            Sctx.add_fastq(sample_info['fastq_paths'])
            #print(Sctx.fastqs)



if __name__ == '__main__':
    try:
        wrapper()
    except KeyboardInterrupt:
        pass
    except Exception:
        pass


