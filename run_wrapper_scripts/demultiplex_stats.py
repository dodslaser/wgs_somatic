#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python

import os
import json
import glob

from collections import defaultdict

def fetch_fastq_path(demux_path, sample_id):
    """Fetch all samples matching regex recursively from demultiplex root."""
    return sorted(glob.glob(f'{demux_path}/**/{sample_id}*.fastq.gz', recursive=True))


def read_lane_stats(stats):
    """Read and format lane stats."""

    lane_stats = []

    for lane in stats:
        lane_number = lane['LaneNumber']
        lane_name = f'lane_{lane_number}'
        reads = lane['TotalClustersPF']
        undetermined_reads = lane['Undetermined']['NumberReads']
        percent_undetermined_reads = round((undetermined_reads/reads) * 100, 2)
        clusters = lane['TotalClustersRaw']
        failed_clusters = clusters - reads
        percent_failed_clusters = round((failed_clusters/clusters) * 100, 2)

        formatted = {'lane_number': lane_number,
                     'lane_name': lane_name,
                     'reads': reads,
                     'undetermined_reads': undetermined_reads,
                     'percent_undetermined_reads': percent_undetermined_reads,
                     'clusters': clusters,
                     'failed_clusters': failed_clusters,
                     'percent_failed_clusters': percent_failed_clusters}

        lane_stats.append(formatted)

    return lane_stats


def read_sample_stats(demux_path, stats):
    """Read and format sample stats."""
    lane_sample_stats = [x['DemuxResults'] for x in stats]

    sample_stats = defaultdict(lambda: defaultdict(dict))

    # First, collect on a lane basis
    for i, lane in enumerate(lane_sample_stats, start=1):
        for sample in lane:
            lane_name = f'lane_{i}'
            sample_name = sample['SampleName']
            reads = sample['NumberReads']
            index7, index5 = sample['IndexMetrics'][0]['IndexSequence'].split('+')  # NOTE Untested. Rewritten on every loop but same in all lanes

            sample_stats[sample_name]['lane_reads'][lane_name] = reads
            sample_stats[sample_name]['index7'] = index7
            sample_stats[sample_name]['index5'] = index5

    # Then, sum all lanes together
    for sample_name, info in sample_stats.items():
        total_reads = sum(info['lane_reads'].values())
        megareads = round(total_reads / 1_000_000, 2)
        expected_coverage = round((total_reads * 300 * 0.88) / 3_137_459_471, 2)  # 0.88 duplicate weight

        info['total_reads'] = total_reads
        info['total_megareads'] = megareads
        info['fastq_paths'] = fetch_fastq_path(demux_path, sample_name)
        info['fastq_dirpath'] = os.path.dirname(info['fastq_paths'][0])
    return sample_stats


def read_overview_stats(lane_stats):
    """Summarize the overall stats from given lane stats."""

    overview_stats = defaultdict(int)

    for lane_info in lane_stats:
        overview_stats['reads'] += lane_info['reads']
        overview_stats['clusters'] += lane_info['clusters']
        overview_stats['undetermined_reads'] += lane_info['undetermined_reads']
        overview_stats['failed_clusters'] += (lane_info['clusters'] - lane_info['reads'])

    overview_stats['percent_undetermined_reads'] = round((overview_stats['undetermined_reads'] / overview_stats['reads']) * 100, 2)
    overview_stats['percent_failed_cluster'] = round((overview_stats['failed_clusters'] / overview_stats['clusters']) * 100, 2)

    return overview_stats

def fetch_stats(dmxjson, demux_path):
    with open(dmxjson, 'r') as inp:
        stats = json.load(inp)['ConversionResults']

    lane_stats = read_lane_stats(stats)
    overview_stats = read_overview_stats(lane_stats)
    sample_stats = read_sample_stats(demux_path, stats)

    master = {'metadata': {},
              'demultiplex': {'overview': overview_stats,
                              'lanes': lane_stats},
              'samples': sample_stats}

    return master
