import os
import re

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction, not_equals
from slims.content import Status



class slims:
    url = os.environ.get('SLIMS_URL')
    user = os.environ.get('SLIMS_USER')
    password = os.environ.get('SLIMS_PASSWORD')


instance = 'wgs-somatic_query'
url = slims.url
user = slims.user
password = slims.password
slims = Slims(instance, url, user, password)



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
            records = slims.fetch('Content', conjunction()
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
        print(self.run_tag)
        print(self.sample_name)
        if not self._fastq:
            records = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22))
                                  .add(equals('cntn_cstm_runTag', self.run_tag)))
            print(records) 
            if len(records) > 1:
                raise Exception('More than 1 fastq somehow.')

            if records:
                self._fastq = records[0]


        # if sample_name has fastqs from additional sequencing runs - fetch those fastq objects
        more_fastqs = slims.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22))
                                  .add(not_equals('cntn_cstm_runTag', self.run_tag)))
        # can find the run tags of fastqs from additional runs. can use this to find the fastqs in demultiplexdir or download them from hcp. need to cat fastqs from different runs before starting pipeline. or maybe it works without merging the fastqs first, need to check this...
        print(more_fastqs)
        for fq in more_fastqs:
            print(fq.cntn_cstm_runTag.value)

        return self._fastq


# Haven't used the bioinformatics property for anything yet. Might not need it. 
#    @property
#    def bioinformatics(self):
#        if not self.run_tag:
#            raise Exception('Can not fetch fastq without a set run tag.')
#
#        if not self._bioinformatics:
#            records = slims.fetch('Content', conjunction()
#                                  .add(equals('cntn_id', self.sample_name))
#                                  .add(equals('cntn_fk_contentType', 23))
#                                  .add(equals('cntn_cstm_runTag', self.run_tag)))
#
#            if len(records) > 1:
#                raise Exception('More than 1 bioinformatics somehow.')
#
#            if records:
#                self._bioinformatics = records[0]
#
#        return self._bioinformatics


def translate_slims_info(record):
    sample_name = record.cntn_id.value

    # 29 = WOPR
    # 54 = wgs_somatic

    pipeline_pks = record.cntn_cstm_secondaryAnalysis.value
    pcr = record.cntn_cstm_pcr.value

    investigator = 'CGG'  # NOTE: Needed?

    department_record = slims.fetch_by_pk('ReferenceDataRecord', record.cntn_cstm_department.value)
    department = department_record.rdrc_name.value  # Format KK
    responder_records = [slims.fetch_by_pk('ReferenceDataRecord', pk) for
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
        'is_trio': is_trio,
        'trio_id': trio_id,
        'trio_role': trio_role,
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
#    print(SSample.dna)
    print(SSample.fastq)
    Sctx.slims_info = translate_slims_info(SSample.dna)
    return



# Can now get from slims:
# Secondary analysis = wgs_somatic
# Run ID for possible additional fastqs for DNA no. 
# Tumor/Normal type
# Tumor ID
# Gender

# Need to do:
# Paired T+N comes from same run = start wgs_somatic
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Sample has been sequenced before and has additional fastqs. Find them and merge fastqs. Might need to download from HCP. Start pipeline with correct pair. The other part of the pair could be in same run or previous run.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
