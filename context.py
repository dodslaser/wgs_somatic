import os

class RunContext:
    def __init__(self, run_path):
        self.run_path = run_path
        self.run_name = os.path.basename(run_path)
        self.run_date = self.run_name.split('_')[0]
        self.run_flowcell = self.run_name.split('_')[-1]
        self.run_tag = '_'.join([self.run_date, self.run_flowcell])
        #self.samplesheet_path = os.path.join(run_path, 'SampleSheet.csv')
