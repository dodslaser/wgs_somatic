#!/apps/bio/software/anaconda2/envs/mathias_general/bin/python3.6
import os
import csv
import sys
import argparse
import statistics

def count_bases_in_threshold(coverage_list, threshold):
    total_positions = len(coverage_list)
    count_positions = 0
    for position in coverage_list:
        if position >= threshold:
            count_positions += 1
    percent_covered = float(count_positions*100) / total_positions
    percent_covered = float("{0:.3f}".format(percent_covered))
    return percent_covered

def calculate_mean(coverage_list) -> float:
    mean_calc = sum(coverage_list) / len(coverage_list) 
    mean = float("{0:.3f}".format(mean_calc))
    return mean

def calculate_median():
    pass

def calculate_stdev(coverage_list):
    stdev_calc = statistics.stdev(coverage_list)
    stdev = float("{0:.3f}".format(stdev_calc))
    return stdev

def overall_coverage_stats(coverage, thresholds) -> list:
    with open(coverage, 'r') as coverage_file:
        # 3rd column is coverage
        coverage_list = []
        for covpos in coverage_file:
            covpos = covpos.rstrip('\n')
            try:
                coverage = covpos.split("\t")[2]
            except:
                print(covpos)
            coverage_list.append(int(coverage))
    
        threshold_list = thresholds.split(",")
        
        cov_stats = {}
        
        cov_stats["mean"] = str(calculate_mean(coverage_list))
        cov_stats["stdev"] = str(calculate_stdev(coverage_list))

        for threshold in threshold_list:
            cov_stats[f"%>={threshold}X"] = str(count_bases_in_threshold(coverage_list, int(threshold)))

        data = []
        header = []
        for coverage_stat in cov_stats:
            header.append(coverage_stat)
            data.append(cov_stats[coverage_stat])
        #write_header = '\t'.join(header)
        #write_data = '\t'.join(data)
    
        #print(write_header)
        #print(write_data)
        return [header, data]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--coverage', nargs='?', help='(from samtools mpileup -l $BED -Q 1 $BAM | cut -f1,2,4)', required=True)
    parser.add_argument('-t', '--thresholds', nargs='?', help='list of coverage thresholds, 1,10,20,40,50 etc...', required=True)
#    parser.add_argument('-o', '--output', nargs='?', help='(Output location)', required=True)
    args = parser.parse_args()
    overall_coverage_stats(args.coverage, args.thresholds)

