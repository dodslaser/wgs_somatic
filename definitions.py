import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(ROOT_DIR, 'configs', 'run-wrapper_config.yaml')
ROOT_LOGGING_PATH = '/medstore/logs/pipeline_logfiles/wgs_somatic'
INSILICO_CONFIG = os.path.join(ROOT_DIR, 'configs', 'insilico_config.json')
INSILICO_PANELS_ROOT='/apps/bio/dependencies/wgs_somatic/hg38/insilico'
