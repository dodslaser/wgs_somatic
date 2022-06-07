#!/bin/python3.6
#import sys
import argparse
from pysam import VariantFile

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("vcf_input", nargs="?", type=str, help=':Full path to your vcf input file')
    parser.add_argument("vcf_output", nargs="?", type=str, help=':Specify full output path and filename of the vcf file you want to create')

    args = parser.parse_args()


    vcf_in = args.vcf_input
    vcf_out = args.vcf_output

    vcf_in = VariantFile(vcf_in)
    new_header = vcf_in.header
    new_header.info.add("DP","1","Integer","Sum of AD fields")
    new_header.info.add("AF","1","Float","Alt AD / sum(AD)")

    vcf_out = VariantFile(vcf_out, 'w', header=new_header)

    for record in vcf_in.fetch():
        dp = sum(record.samples[0].get("AD"))
        record.info["DP"]=dp
        af=record.samples[0].get("AD")[1]/dp
        record.info["AF"]=af
        vcf_out.write(record)


if __name__ == "__main__":
    main()

