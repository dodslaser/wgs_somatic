import os
import re
import json

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction, not_equals
from slims.content import Status


class slims_credentials:
    url = os.environ.get('SLIMS_URL')
    user = os.environ.get('SLIMS_USER')
    password = os.environ.get('SLIMS_PASSWORD')


instance = 'wgs-somatic_query'
url = slims_credentials.url
user = slims_credentials.user
password = slims_credentials.password
slims_connection = Slims(instance, url, user, password)



class SlimsSample:
    def __init__(self, sample_name, run_tag=''):
        self.sample_name = sample_name
        self.run_tag = run_tag

        self._dna = None
        self._fastq = None
        self._bioinformatics = None

    @property
    def dna(self):
        if not self._dna:
            records = slims_connection.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 6)))

            if len(records) > 1:
                raise Exception('More than 1 DNA somehow.')

            if records:
                #print(records)
                self._dna = records[0]

        return self._dna

    @property
    def fastq(self):
        if not self.run_tag:
            raise Exception('Can not fetch fastq without a set run tag.')
        #print(self.run_tag)
        #print(self.sample_name)
        if not self._fastq:
            records = slims_connection.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22))
                                  .add(equals('cntn_cstm_runTag', self.run_tag)))
            #print(records) 
            if len(records) > 1:
                raise Exception('More than 1 fastq somehow.')

            if records:
                self._fastq = records[0]

        return self._fastq




def translate_slims_info(record):
    sample_name = record.cntn_id.value

    # 29 = WOPR
    # 54 = wgs_somatic

    pipeline_pks = record.cntn_cstm_secondaryAnalysis.value
    pcr = record.cntn_cstm_pcr.value

    investigator = 'CGG'  # NOTE: Needed?

    department_record = slims_connection.fetch_by_pk('ReferenceDataRecord', record.cntn_cstm_department.value)
    department = department_record.rdrc_name.value  # Format KK
    responder_records = [slims_connection.fetch_by_pk('ReferenceDataRecord', pk) for
                         pk in department_record.rdrc_cstm_responder.value]
    responder_emails = [rec.rdrc_cstm_email.value for rec in responder_records]

    is_research = record.cntn_cstm_research.value
    research_project = record.cntn_cstm_researchProject.value

    is_priority = True if record.cntn_cstm_priority.value else False

    gender = record.gender.value

    is_trio = record.cntn_cstm_trio.value
    trio_id = record.cntn_cstm_trioID.value
    trio_role = record.cntn_cstm_trioRole.value

    tumorNormalType = record.cntn_cstm_tumorNormalType.value
    tumorNormalID = record.cntn_cstm_tumorNormalID.value

    tertiary_analysis = record.cntn_cstm_tertiaryAnalysis.value

    master = {
        'content_id': sample_name,
        'investigator': investigator,
        'department': department,
        'responder_mails': responder_emails,
        'is_research': is_research,
        'research_project': research_project,
        'gender': gender,
        'is_priority': is_priority,
        'pcr': pcr,
        'tumorNormalType': tumorNormalType,
        'tumorNormalID': tumorNormalID,
        'secondary_analysis': pipeline_pks,
        'tertiary_analysis': tertiary_analysis
    }
    return master


def get_sample_slims_info(Sctx, run_tag):
    """Query the slims API for relevant metadata given sample_name in samplesheet."""
    SSample = SlimsSample(Sctx.sample_name, run_tag)

    if not SSample.dna:
        Sctx.slims_info = {}
        return
    #print(SSample.dna)
    #print(SSample.fastq)
    Sctx.slims_info = translate_slims_info(SSample.dna)
    return


def more_fastqs(sample_name, run_tag):
    """If a sample name has fastqs from additional sequencing runs - fetch those fastq objects. Returns a dictionary of sample name and location of its fastqs from other runs. """
    more_fastqs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_id', sample_name))
                              .add(equals('cntn_fk_contentType', 22))
                              .add(not_equals('cntn_cstm_runTag', run_tag)))
    #print('hej')
    #print(more_fastqs)
    if more_fastqs:
        runtags = []
        more_fq_paths = []
        d = dict()
        run_paths = []
        for fq in more_fastqs:
            #print(fq)
            fqs_runtag = fq.cntn_cstm_runTag.value
            #print(f'fastqs runtags: {fqs_runtag}')
            #print(Sctx.sample_name)
            runtags.append(fqs_runtag)
        for tag in runtags:
            fqSSample = SlimsSample(sample_name, tag)
            #print(fqSSample.dna)
            #print(fqSSample.fastq)
            #print(fqSSample.fastq.cntn_cstm_demuxerSampleResult.value)
            json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerSampleResult.value)
            fq_paths = json_info['fastq_paths']
            #print(f'fq_paths: {fq_paths}')
            
            #json_backup_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerBackupSampleResult.value)
            #print(json_backup_info['remote_keys'])
            # remote_key is the same as run name. can use this to be able to create run context for additional runs and can then create sample contexts for samples in those runs. Or not.. need whole path to demultiplexdir
            #for remote_key in json_backup_info['remote_keys']:
                #print(remote_key.split('/')[0])
            #    remote_keys.append(remote_key.split('/')[0])
            
            for path in fq_paths:
                #print(f' path: {path.split("/fastq/")[0]}')
                more_fq_paths.append(path) 
                new_path = path.split("/fastq/")[0]
                run_paths.append(new_path) if new_path not in run_paths else run_paths
                #run_paths.append(path.split("/fastq/")[0])
        #print(sample_name)
        #print(more_fq_paths)
        d['Sample ID'] = sample_name
        d['fastq paths'] = more_fq_paths
        #print(f'FQ PATHS:{more_fq_paths} ')
        for f_path in more_fq_paths:
        # Only links if link doesn't already exist
            if not os.path.islink(os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(f_path))):
        # Now symlinks all additional paths to fastqs for tumor and normal in other runs. If I symlink to demultiplexdir of particular run instead, all fastqs belonging to the T/N pair will be in the same folder and the pipeline can start using that folder as argument.
                os.symlink(f_path, os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(f_path)))
        # Now returns a dict with sample name and fastqs from other runs belonging to that sample name. Might not be used now that the symlinks are created. But need the sample name of for example normal if only the tumor is in the run in question and vice versa. So could sort of be used.
        #return d
        #print(f'printed {run_paths}')
        return run_paths

def get_pair(Sctx, run_tag):
    """If tumor and normal are sequenced in different runs - find the pairs. Then use the more_fastqs function to find paths of fastqs that are sequenced in different runs. """

    get_sample_slims_info(Sctx, run_tag)
    #print(Sctx.slims_info)
    #print(Sctx.slims_info['tumorNormalID'])
    pairs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_fk_contentType', 6))
                              .add(equals('cntn_cstm_tumorNormalID', Sctx.slims_info['tumorNormalID'])))
    #print('hej')
    #print(f'pairs: {pairs}')
    parts_of_pair = []
    fqs =[]
    for pair in pairs:
        #print(f'pair {pair}')
        #print(f' t/n type: {pair.cntn_cstm_tumorNormalType.value}')
        #print(f'sample id: {pair.cntn_id.value}')
        #print(f'run tag: {run_tag}')
        #print(more_fastqs(pair.cntn_id.value, run_tag))
        # if there are not fastqs in other runs, skip!
        if not more_fastqs(pair.cntn_id.value, run_tag) == None:
            fqs.append(more_fastqs(pair.cntn_id.value, run_tag))
    # if fqs are in other run, get those paths:
    if fqs:
        return fqs


# Can now get from slims:
# Secondary analysis = wgs_somatic
# Run ID for possible additional fastqs for DNA no. 
# Tumor/Normal type
# Tumor ID
# Gender

# Need to do:
# Paired T+N comes from same run = start wgs_somatic
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Sample has been sequenced before and has additional fastqs. Find them and merge fastqs. Might need to download from HCP. Start pipeline with correct pair. The other part of the pair could be in same run or previous run. Can now find additional fastqs from other runs. Fastqs don't need to be merged to start the pipeline, they just need to be in the same folder. Mapping will be done for each fastq R1+R2 and then the bam files will be merged in the dedup step.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
