"""
Wrapper to be used by cron for automatic start of wgs_somatic
"""

import sys
import argparse
import os
import re
import glob
import yaml
from datetime import datetime
import json
from itertools import chain
import traceback
import subprocess
import threading

from definitions import CONFIG_PATH, ROOT_DIR, ROOT_LOGGING_PATH
from context import RunContext, SampleContext
from helpers import setup_logger
from tools.slims import get_sample_slims_info, SlimsSample, find_more_fastqs, get_pair_dict
from tools.email import start_email, end_email
from launch_snakemake import analysis_main, petagene_compress_bam, yearly_stats


# Store info about samples to use for sending report emails
sample_status = {'missing_slims': [],
                 'unset_WS': [],
                 'approved': []}


def look_for_runs(config, instrument):
    '''Look for runs in demultiplexdir'''
    instrument_root_path = config[instrument]['demultiplex_path']
    found_paths = glob.glob(os.path.join(instrument_root_path, '*'))
    regex = config[instrument]['seq_name_regex']
    return [path for path in found_paths if re.search(regex, os.path.basename(path))]


def generate_context_objects(Rctx, logger):
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
        # exit script if no samples set for wgs_somatic in run
        sys.exit()

    for Sctx in Rctx.sample_contexts:
        sample_status['approved'].append(Sctx)

    return Rctx

def call_script(**kwargs):
    '''Function to call main function from launch_snakemake.py'''
    args = argparse.Namespace(**kwargs)
    subprocess.call(analysis_main(args, **kwargs))

def check_ok(outputdir):
    '''Function to check if analysis has finished correctly'''

    if os.path.isfile(f"{outputdir}/reporting/workflow_finished.txt"):
        return True
    else:
        return False

def analysis_end(outputdir, tumorsample, normalsample):
    '''Function to check if analysis has finished correctly and add to yearly stats and start petagene compression'''

    if os.path.isfile(f"{outputdir}/reporting/workflow_finished.txt"):
        if tumorsample:
            # these functions are only executed if snakemake workflow has finished successfully
            yearly_stats(tumorsample, normalsample)
            petagene_compress_bam(outputdir, tumorsample)
        else:
            yearly_stats('None', normalsample)
            petagene_compress_bam(outputdir, normalsample)

def wrapper(instrument):
    '''Wrapper function'''
    logger = setup_logger('wrapper', os.path.join(ROOT_LOGGING_PATH, f'{instrument}_WS_wrapper.log'))

    # Empty dict, will update later with T/N pair info
    pair_dict_all_pairs = {}

    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)

    # Grab all available local run paths
    local_run_paths = look_for_runs(config, instrument)
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
        Rctx_run = generate_context_objects(Rctx, logger)

        # Get T/N pair info in a dict for samples and link additional fastqs from other runs
        for sctx in Rctx_run.sample_contexts:
            pair_dict = get_pair_dict(sctx, Rctx, logger)
            pair_dict_all_pairs.update(pair_dict)

        # Uses the dictionary of T/N samples to put the correct pairs together and finds the correct input arguments to the pipeline
        threads = []
        check_ok_outdirs = []
        end_threads = []
        final_pairs = []
        for key in pair_dict_all_pairs:
            if 'tumor' in pair_dict_all_pairs.get(key):
                t = key
                # Using the list containing values 'tumor', value of tumorNormalID, department and prio status
                # Removing the value 'tumor' from the list to get the tumorNormalID
                # TODO: Would be nice to do in a better way rather than using [0] to get the first value in the list
                t_ID = [val for val in pair_dict_all_pairs.get(key) if val != 'tumor'][0] 
                for k in pair_dict_all_pairs:
                    if 'normal' in pair_dict_all_pairs.get(k):
                        n = k
                        n_ID = [val for val in pair_dict_all_pairs.get(k) if val != 'normal'][0]
                        department = [val for val in pair_dict_all_pairs.get(k) if val != 'normal'][1]
                        is_prio = [val for val in pair_dict_all_pairs.get(k) if val != 'normal'][2]
                        if is_prio:
                            prio_sample = 'prio'
                        else:
                            prio_sample = ''
                        # As of now, tumorNormalID is the same for tumor and normal.
                        # In the future, this will be changed to pairID
                        # The or statements are here to prepare to when we change to pair ID
                        # Pair ID for tumor will be normal name (minus DNA) and the opposite for normal
                        if n_ID == t_ID or t_ID == n.split("DNA")[1] or n_ID == t.split("DNA")[1]: 
                            runnormal = Rctx_run.run_name
                            runtumor = Rctx_run.run_name
                            tumorsample = t
                            normalsample = n
                            normalfastqs = os.path.join(Rctx_run.run_path, "fastq")
                            tumorfastqs = normalfastqs
                            outputdir = os.path.join(config['outputdir']['GMS-BT'], tumorsample) 
                            #outputdir = os.path.join("/home/xshang/ws_testoutput/outdir/", tumorsample) #use for testing
                            igvuser = config['igv']['GMS-BT']
                            #igvuser = 'alvar.almstedt' # use for testing
                            # FIXME Use boolean values instead of 'yes' for hg38ref and handle the translation later on
                            hg38ref = config['hg38ref']['GMS-BT']

                            
                            # Use this list of final pairs for email
                            final_pairs.append(f'{tumorsample} (T) {normalsample} (N), {department} {prio_sample}')

                            # If sample has been run before, outdir already exists. Changing the name of old outdir to make room for new outdir. Should maybe move old outdir to archive instead.
                            # Won't work if outputdir_old also already exists. Need to be solved in a better way 
                            if os.path.exists(outputdir):
                                logger.info(f'Outputdir exists for {tumorsample}. Renaming old outputdir {outputdir} to {outputdir}_old')
                                os.rename(outputdir, f'{outputdir}_old')

                            pipeline_args = {'runnormal': f'{runnormal}', 'output': f'{outputdir}', 'normalname': f'{normalsample}', 'normalfastqs': f'{normalfastqs}', 'runtumor': f'{runtumor}', 'tumorname': f'{tumorsample}', 'tumorfastqs': f'{tumorfastqs}', 'igvuser': f'{igvuser}', 'hg38ref': f'{hg38ref}'}


                            # Using threading to start the pipeline for several samples at the same time
                            threads.append(threading.Thread(target=call_script, kwargs=pipeline_args))
                            logger.info(f'Starting wgs_somatic with arguments {pipeline_args}')

                            check_ok_outdirs.append(outputdir)
                            end_threads.append(threading.Thread(target=analysis_end, args=(outputdir, tumorsample, normalsample)))


        # Start several samples at the same time
        for t in threads:
            t.start()
            logger.info(f'Thread {t} has started')

        # Send start email
        start_email(Rctx_run.run_name, final_pairs)

        for u in threads:
            u.join()
            logger.info(f'Thread {u} is finished')

        # Check if all samples in run have finished successfully. If not, exit script and send error email.
        for outdir in check_ok_outdirs:
            if check_ok(outdir) == True:
                continue
            else:
                logger.info('All jobs have not finished successfully')
                # add something to send end email about fail
                sys.exit()

        
        logger.info('All jobs have finished successfully')
        end_email(Rctx_run.run_name, final_pairs)

        # If all jobs have finished successfully - add to yearly stats and start petagene compression of bamfiles
        for t in end_threads:
            t.start()
        for u in end_threads:
            u.join()

        # break out of for loop to avoid starting pipeline for a possible other run that was done sequenced at the same time.
        # if cron runs every 30 mins it will find other runs at the next cron instance and run from there instead (and add to novaseq_runlist)
        break




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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--instrument', help='For example novaseq_687_gc or novaseq_A01736', required=True)
    args = parser.parse_args()

    wrapper(args.instrument)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
    #except Exception:
        #format_exc = traceback.format_exc()
        #logger.error(format_exc)
        # TODO: add send email about error here
