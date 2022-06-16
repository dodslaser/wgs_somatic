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
        #for line in bed:
        #    print(line.split(" ")[3])

        genesDup = [line.split(" ")[3] for line in bed]
        genes = set(genesDup)

    #print(genes)

    # Write genes included in excel file
    worksheet.write('A6', 'Genes included: ')
    row = 7
    for gene in genes:
        worksheet.write('A'+str(row), gene)
        row += 1

    row += 1
    #print(list(vcf_input.header.samples))
    #print(len(list(vcf_input.header.samples)))
    samples = list(vcf_input.header.samples)
    # If only tumor sample in vcf or if tumor + normal sample in vcf...
    if len(samples) == 1:
        sample = samples[0]
        tableheading = ['Chr', 'Position', 'Ref', 'Alt', sample]
        worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)
    if len(samples) == 2:
        sample1 = samples[0]
        sample2 = samples[1]
        tableheading = ['Chr', 'Start', 'Stop', 'Ref', 'Alt', sample1, sample2]
        worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)
        for indel in vcf_input.fetch():
            if len(indel.alts) == 1:
                alt=indel.alts[0]
            #chrom = indel.contig
            worksheet.write_row(row,0,[indel.contig, indel.pos, indel.stop, indel.ref, alt])
            row += 1

    #row += 1
    #worksheet.write_row('A'+str(row), tableheading, tableHeadFormat)


    #for record in vcf_input.fetch():
    #    print(record)
    #    print(record.info.keys())


    workbook.close()

if __name__ == "__main__":
    main()
