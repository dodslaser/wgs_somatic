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
        if not self._fastq:
            records = slims_connection.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22))
                                  .add(equals('cntn_cstm_runTag', self.run_tag)))
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
    department = department_record.rdrc_name.value  
    responder_records = [slims_connection.fetch_by_pk('ReferenceDataRecord', pk) for
                         pk in department_record.rdrc_cstm_responder.value]
    responder_emails = [rec.rdrc_cstm_email.value for rec in responder_records]

    is_research = record.cntn_cstm_research.value
    research_project = record.cntn_cstm_researchProject.value

    is_priority = True if record.cntn_cstm_priority.value else False

    gender = record.gender.value

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
    return translate_slims_info(SSample.dna)

def link_fastqs(list_of_fq_paths):
    '''Link fastqs to fastq-folder in demultiplexdir of current run. Need to change the hardcoded path to my home... '''
    # TODO: using a hardcoded test folder right now for symlinks. will change this to correct Demultiplexdir/current-run/fastq folder.
    # TODO: additional fastqs need to still be in demultiplexdir. not considering downloading from hcp right now. need to consider this later...
    for fq_path in list_of_fq_paths:
    # Only links if link doesn't already exist
        if not os.path.islink(os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(fq_path))):
        # Now symlinks all additional paths to fastqs for tumor and normal in other runs. If I symlink to demultiplexdir of particular run instead, all fastqs belonging to the T/N pair will be in the same folder and the pipeline can start using that folder as argument.
            os.symlink(fq_path, os.path.join(f"/home/xshang/ws_testoutput/symlinks/", os.path.basename(fq_path)))

def run_paths_for_more_fastqs(sample_name, run_tag):
    """
    If a sample name has fastqs from additional sequencing runs - fetch those fastq objects. 
    Returns run paths (Demultiplexdir/runTag) of runs where you can find additional fastqs. 
    """
    more_fastqs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_id', sample_name))
                              .add(equals('cntn_fk_contentType', 22))
                              .add(not_equals('cntn_cstm_runTag', run_tag)))
    if more_fastqs:
        runtags = []
        run_paths = []
        for fq in more_fastqs:
            fqs_runtag = fq.cntn_cstm_runTag.value
            runtags.append(fqs_runtag)
        #print(f'runtag: {runtags}')
        for tag in runtags:
            #print(f'tag: {tag}')
            fqSSample = SlimsSample(sample_name, tag)
            json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerSampleResult.value)
            fq_paths = json_info['fastq_paths']
            #print(f'fq paths: {fq_paths}')
            #print(f'linking fastqs for {sample_name}_{tag}') 
            link_fastqs(fq_paths)
            for path in fq_paths:
                new_path = path.split("/fastq/")[0]
                if new_path not in run_paths:
                    run_paths.append(new_path)
        return run_paths

def get_pair_and_run_paths(Sctx, run_tag):
    """
    If tumor and normal are sequenced in different runs - find the pairs. 
    Then use the run_paths_for_more_fastqs function to find paths of fastqs that are sequenced in different runs. 
    """

    pair_dict = {}

    # FIXME: using equals tumorNormalID here won't work when we change it to pairID...
    Sctx.slims_info = get_sample_slims_info(Sctx, run_tag)
    pairs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_fk_contentType', 6))
                              .add(equals('cntn_cstm_tumorNormalID', 
                              Sctx.slims_info['tumorNormalID'])))
    run_paths =[]
    #print(f'pairs: {pairs}')
    for pair in pairs:
        pair.slims_info = translate_slims_info(pair)
        pair_dict[pair.slims_info["content_id"]] = [pair.slims_info["tumorNormalType"], pair.slims_info["tumorNormalID"]]
        #pair_dict[pair.slims_info["content_id"]].append(pair.slims_info["tumorNormalID"])
        #print(f'run path for more fqs: {run_paths_for_more_fastqs(pair.cntn_id.value, run_tag)}')
        # if there are not fastqs in other runs, skip!
        #r_paths = run_paths_for_more_fastqs(pair.cntn_id.value, run_tag)
        #if r_paths:
            #print(f'sample {pair.cntn_id.value} has additional fqs in run {r_paths}')
            #run_paths.append(r_paths)
    #print(f'pair dict: {pair_dict}')
    # if fqs are in other run, get those paths:
    #return run_paths or None
    return pair_dict


# Need to do:
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
