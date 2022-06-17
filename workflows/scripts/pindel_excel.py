#!/bin/python3.6
#import sys
import argparse
from pysam import VariantFile
import xlsxwriter

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("vcf_input", nargs="?", type=str, help=':Full path to your vcf input file')
    parser.add_argument("xlsx_output", nargs="?", type=str, help=':Specify full output path and filename of the vcf file you want to create')
    parser.add_argument("bedfile", nargs="?", type=str, help=':Full path to bedfile used by pindel')

    args = parser.parse_args()

    # Create excel file
    workbook = xlsxwriter.Workbook(args.xlsx_output)

    # Define formats to be used.
    headingFormat = workbook.add_format({'bold': True, 'font_size': 18})
    tableHeadFormat = workbook.add_format({'bold': True, 'text_wrap': True})

    # Create worksheet
    worksheet = workbook.add_worksheet('Pindel')
    worksheet.write('A1', 'Pindel results', headingFormat)

    vcf_input = args.vcf_input
    vcf_input = VariantFile(vcf_input, 'r')

    # Get reference used
    for x in vcf_input.header.records:
        if x.key=='reference':
            ref = x.value

    worksheet.write('A4', 'Reference used: '+str(ref))


    # Get genes in bedfile
    with open(args.bedfile) as bed:
        genesDup = [line.split(" ")[3] for line in bed]
        genes = set(genesDup)

    # Get position interval for each gene
    #with open(args.bedfile) as bed:
    #    for line in bed:
    #        for gene in genes:
    #            if line.split(" ")[3] == gene:
    #                print(line)
    def position_gene(position):
        start = []
        stop = []
        with open(args.bedfile) as bed:
            for line in bed:
                if position >= int(line.split(" ")[1]) and position <= int(line.split(" ")[2]):
                    gene = line.split(" ")[3]
        return(gene)
    print(position_gene(28034106))

    #for gene in genes:
    #    print(position_interval(gene))
    #position_interval('FLT3')

    # Write genes included in excel file
    worksheet.write('A6', 'Genes included: ')
    row = 7
    for gene in genes:
        worksheet.write('A'+str(row), gene)
        row += 1

    row += 1
    samples = list(vcf_input.header.samples)
    # If only tumor sample in vcf or if tumor + normal sample in vcf...
    # Write info from vcf to excel file
    if len(samples) == 1:
        sample = samples[0]
        tableheading = ['Chr', 'Position', 'Ref', 'Alt', sample]
        worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)
    if len(samples) == 2:
        sample1 = samples[0]
        sample2 = samples[1]
        tableheading = ['Chr', 'Start', 'Stop', 'SV length', 'Ref', 'Alt', sample1+' DP', sample1+' AF', sample2+' DP', sample2+' AF']
        worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)
        for indel in vcf_input.fetch():
            svlen = indel.info["SVLEN"]
            if len(indel.alts) == 1:
                alt=indel.alts[0]
            s1_dp = indel.samples[sample1].get("DP")
            s1_af = indel.samples[sample1].get("AF")
            s2_dp = indel.samples[sample2].get("DP")
            s2_af = indel.samples[sample2].get("AF")
            worksheet.write_row(row,0,[indel.contig, indel.pos, indel.stop, svlen, indel.ref, alt, s1_dp, s1_af, s2_dp, s2_af])
            row += 1



    workbook.close()

if __name__ == "__main__":
    main()
