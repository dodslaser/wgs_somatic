"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import os
import re
import glob
import yaml

from definitions import CONFIG_PATH, ROOT_DIR
from context import RunContext

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


if __name__ == '__main__':
    try:
        wrapper()
    except KeyboardInterrupt:
        pass
    except Exception:
        pass


