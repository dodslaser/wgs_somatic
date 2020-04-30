#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
import yaml
import json

def read_config():
    with open("configs/config.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data

#def read_inputfile():
#    with open("configs/leukemi143.yaml") as inputfile:
#        inputfile_info = yaml.load(inputfile, Loader=yaml.FullLoader)
#        return inputfile_info

def read_clusterconf():
    with open("configs/cluster.yaml") as inputfile:
        inputfile_info = yaml.load(inputfile, Loader=yaml.FullLoader)
        return inputfile_info

def read_passconfig():
    with open("/root/password_config.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data
