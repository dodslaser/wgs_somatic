#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python

from run_wrapper_scripts.run_helpers import read_config

def organise_analysisdict(sampledict):
    config = read_config()
    allowed_pipelines = config["analysis"]["allowed_pipelines"]    
    allowed_projects = config["analysis"]["allowed_projects"]
    analysisdict = {}
    for sample in sampledict["samples"]:
        sample_info = sampledict["samples"]["sample"]:
        project = sample_info["project"]
        pipeline = sample_info["analysis_pipeline"]
        
        if project in allowed_projects:
            if pipeline in allowed_pipelines:
                
