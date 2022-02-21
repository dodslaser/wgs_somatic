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
    Sctx.slims_info = translate_slims_info(SSample.dna)
    return


def more_fastqs(sample_name, run_tag):
    """
    If a sample name has fastqs from additional sequencing runs - fetch those fastq objects. 
    Returns run paths of runs where you can find additional fastqs. 
    """
    more_fastqs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_id', sample_name))
                              .add(equals('cntn_fk_contentType', 22))
                              .add(not_equals('cntn_cstm_runTag', run_tag)))
    if more_fastqs:
        runtags = []
        more_fq_paths = []
        run_paths = []
        for fq in more_fastqs:
            fqs_runtag = fq.cntn_cstm_runTag.value
            runtags.append(fqs_runtag)

        for tag in runtags:
            fqSSample = SlimsSample(sample_name, tag)
            json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerSampleResult.value)
            fq_paths = json_info['fastq_paths']
            
            for path in fq_paths:
                more_fq_paths.append(path) 
                new_path = path.split("/fastq/")[0]
                if new_path not in run_paths:
                    run_paths.append(new_path)
        return run_paths

def get_pair_and_more_fqs(Sctx, run_tag):
    """
    If tumor and normal are sequenced in different runs - find the pairs. 
    Then use the more_fastqs function to find paths of fastqs that are sequenced in different runs. 
    """

    get_sample_slims_info(Sctx, run_tag)
    pairs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_fk_contentType', 6))
                              .add(equals('cntn_cstm_tumorNormalID', 
                              Sctx.slims_info['tumorNormalID'])))
    parts_of_pair = []
    fqs =[]
    for pair in pairs:
        # if there are not fastqs in other runs, skip!
        if more_fastqs(pair.cntn_id.value, run_tag) != None:
            fqs.append(more_fastqs(pair.cntn_id.value, run_tag))
    # if fqs are in other run, get those paths:
    return fqs or None



# Need to do:
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
