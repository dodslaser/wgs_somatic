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

    vcf_in = VariantFile(vcf_in, 'r')
    new_header = vcf_in.header

    # Add DP and AD to FORMAT column
    new_header.formats.add("DP","1","Integer","Sum of AD fields")
    new_header.formats.add("AF","1","Float","Alt AD / sum(AD)")

    vcf_out = VariantFile(vcf_out, 'w', header=new_header)

    for record in vcf_in.fetch():
        # Calculate DP for tumor
        dp = sum(record.samples[0].get("AD"))
        # Add value of DP to tumor sample column
        record.samples[0]["DP"]=dp
        # Calculate AF for tumor
        if dp == 0:
            # Solve division by zero problem if DP = 0
            af = 0
        else:
            af=record.samples[0].get("AD")[1]/dp
        # Add value of AF to tumor sample column
        record.samples[0]["AF"]=af

        # Calculate DP for normal
        n_dp = sum(record.samples[1].get("AD"))
        # Add value of DP to normal sample column
        record.samples[1]["DP"]=n_dp
        # Calculate AF for normal
        if n_dp == 0:
            # Solve division by zero problem if DP = 0
            n_af=0
        else:
            n_af=record.samples[1].get("AD")[1]/n_dp
        # Add value of AF to normal sample column
        record.samples[1]["AF"]=n_af

        vcf_out.write(record)

if __name__ == "__main__":
    main()

