#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import json
import argparse
import os
import glob 
import helpers
import sys
import time
import traceback
from shutil import copyfile
import subprocess

def read_ivaconf():
    with open("configs/ingenuity.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data

def read_wrapperconf():
    with open("configs/wrapper_conf.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data

def get_time():
    nowtime = time.strftime("%Y-%m-%d-%H-%M-%S")
    return nowtime

def logger(message, logfile=False):
    config = read_wrapperconf()
    logdir = config["logdir"]
    current_date = time.strftime("%Y-%m-%d")
    if not logfile:
        logname = f"{logdir}/{current_date}.log"
        logfile = open(logname, "a+")
        logfile.write(f"{get_time()}: {message}" + "\n")
    else:
        logfile = open(logfile, "a+")
        logfile.write(f"{get_time()}: {message}" + "\n")
    print(message)

def analysis_main(args, output, normalname, normalfastqs, tumorname, tumorfastqs, ivauser="", igvuser="", hg38ref=""):
    try:
        ################################################################
        # Write InputArgs to logfile
        config = read_wrapperconf()
        commandlogs = config["commandlogs"]
        command = f"{sys.argv[0]}"
        current_date = time.strftime("%Y-%m-%d")
        commandlog = f"{commandlogs}/commands_{current_date}.log"
        for arg in vars(args):
            command = f"{command} --{arg} {getattr(args, arg)}"
        commandlogfile = open(commandlog, "a+")
        commandlogfile.write(f"{get_time()}" + "\n")
        commandlogfile.write(command + "\n")
        ################################################################

        if output.endswith("/"):
            output = output[:-1]
        if normalfastqs.endswith("/"):
            normalfastqs = normalfastqs[:-1]
        if tumorfastqs.endswith("/"):
            tumorfastqs = tumorfastqs[:-1]

        
        #################################################################
        # Validate Inputs
        ################################################################
        error_list = []

        if hg38ref:
            logger(f"hg38 argument given with value: {hg38ref}")
            if hg38ref != "yes":
                logger("argument is not yes, if you want hg19 simply dont provide hg38 argument, exiting")
                error_list.append(f"invalid hg38 argument value: {hg38ref}")
        
        if hg38ref == "yes":
            mainconf = "hg38conf"
        else:
            mainconf = "hg19conf"
        configdir = config["configdir"]
        mainconf_name = config[mainconf]
        mainconf_path = f"{configdir}/{mainconf_name}"

        # validate fastqdirs
        if not os.path.isdir(normalfastqs):
            error_list.append(f"{normalfastqs} does not appear to be a directory")
        else:
            f_normalfastqs = glob.glob(f"{normalfastqs}/*fastq.gz")
            if not f_normalfastqs:
                error_list.append(f"No fastqs found in normaldir")
        if not os.path.isdir(tumorfastqs):
            error_list.append(f"{tumorfastqs} does not appear to be a directory")
        else:
            f_tumorfastqs = glob.glob(f"{tumorfastqs}/*fastq.gz")
            if not f_tumorfastqs:
                error_list.append(f"No fastqs found in tumordir")
        
        # validate iva and igv users if supplied
        if igvuser:
            mainconf = helpers.read_config(mainconf_path)
            igvdatadir = mainconf["rules"]["share_to_igv"]["igvdatadir"]
            if not os.path.isdir(f"{igvdatadir}/{igvuser}"):
                error_list.append(f"{igvuser} does not appear to be a valid preconfigured IGV user")
        if ivauser:
            ivaconf = read_ivaconf()
            if ivauser not in ivaconf["ivausers"]:
                error_list.append(f"{ivauser} is not a valid preconfigured IVA user")

        # prepare outputdirectory
        if not os.path.isdir(output):
            try:
                os.mkdir(output)
            except Exception as e:
                error_list.append(f"outputdirectory: {output} does not exist and could not be created")

        if error_list:
                logger("Errors found in arguments to script:")
                for arg in vars(args):
                    logger(f"{arg} = {getattr(args, arg)}")
                for error in error_list:
                    logger(error)
                logger("Exiting")
                sys.exit()

        #################################################################
        # Prepare AnalysisFolder
        #################################################################
        samplelogs = f"{output}/logs"
        if not os.path.isdir(samplelogs):
            os.mkdir(samplelogs)
        runconfigs = f"{output}/configs"
        if not os.path.isdir(runconfigs):
            os.mkdir(runconfigs)
        
        # copying configfiles to analysisdir
        clusterconf = config["clusterconf"]
        copyfile(f"{configdir}/{clusterconf}", f"{runconfigs}/{clusterconf}")
        copyfile(f"{configdir}/{mainconf_name}", f"{runconfigs}/{mainconf_name}")

        samplelog = f"{samplelogs}/{tumorname}.log"
        logger("Input validated:", samplelog)
        logger(f"{command}", samplelog)
        logger("Fastqs found for normal:", samplelog)
        logger(f"{f_normalfastqs}", samplelog)
        logger("Fastqs found for tumor:", samplelog)
        logger(f"{f_tumorfastqs}", samplelog)
        
        ##################################################################
        # Create AnalysisConfigfile
        ##################################################################
        analysisdict = {}
        analysisdict["normalname"] = normalname 
        analysisdict["normalfastqs"] = normalfastqs
        analysisdict["tumorname"] = tumorname
        analysisdict["tumorfastqs"] = tumorfastqs
        analysisdict["igvuser"] = igvuser
        analysisdict["ivauser"] = ivauser

        if hg38ref == "yes":
            analysisdict["reference"] = "hg38"
        else:
            analysisdict["reference"] = "hg19"

        with open(f"{runconfigs}/{tumorname}_config.json", 'w') as analysisconf:
            json.dump(analysisdict, analysisconf, ensure_ascii=False, indent=4)
        ###################################################################
        # Start SnakeMake pipeline
        ###################################################################
        snakemake_path = config["snakemake_env"]
        os.environ["PATH"] += os.pathsep + snakemake_path
        my_env = os.environ.copy() 
        snakemake_args = f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{tumorname}_config.json --dag | dot -Tsvg > /seqstore/webfolders/wgs/admin/barncancer/dags/dag_{current_date}.svg"
        # >>>>>>>>>>>> Create Dag of pipeline
        subprocess.run(snakemake_args, shell=True, env=my_env) # CREATE DAG
        snakemake_args = f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{tumorname}_config.json --cluster-config configs/cluster.yaml --cluster \"qsub -S /bin/bash -pe mpi {{cluster.threads}} -q {{cluster.queue}} -N {{cluster.name}} -o {samplelogs}/{{cluster.output}} -e {samplelogs}/{{cluster.error}} -l {{cluster.excl}}\" --jobs 999 --latency-wait 60 --directory {output} &>> {samplelog}"
        # >>>>>>>>>>>> Start pipeline
        subprocess.run(snakemake_args, shell=True, env=my_env) # Shellscript pipeline
        #subprocess.run(f"snakemake -s pipeline.snakefile --configfile {runconfigs}/{tumorname}_config.json --cluster-config configs/cluster.yaml --cluster \"qsub -S /bin/bash -pe mpi " + "{cluster.threads} -q {cluster.queue} -N {cluster.name} -o " + samplelogs + "{cluster.output}"  , shell=True, env=my_env) # Shellscript pipeline
    except Exception as e:
        tb = traceback.format_exc()
        logger(f"Error in script:")
        logger(f"{e} Traceback: {tb}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outputdir', nargs='?', help='output directory, where to put results', required=True)
    parser.add_argument('-ns', '--normalsample', nargs='?', help='normal samplename', required=True)
    parser.add_argument('-nf', '--normalfastqs', nargs='?', help='path to directory containing normal fastqs', required=True)
    parser.add_argument('-tn', '--tumorsample', nargs='?', help='tumor samplename', required=True)
    parser.add_argument('-tf', '--tumorfastqs', nargs='?', help='path to directory containing tumor fastqs', required=True)
    parser.add_argument('-iva', '--ivauser', nargs='?', help='location to output results', required=False)
    parser.add_argument('-igv', '--igvuser', nargs='?', help='location to output results', required=False)
    parser.add_argument('-hg38', '--hg38ref', nargs='?', help='run analysis on hg38 reference (write yes if you want this option)', required=False)
    args = parser.parse_args()
    analysis_main(args, args.outputdir, args.normalsample, args.normalfastqs, args.tumorsample, args.tumorfastqs, args.ivauser, args.igvuser, args.hg38ref)
