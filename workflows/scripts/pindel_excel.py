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

    worksheet = workbook.add_worksheet('Pindel')
    worksheet.write('A1', 'Pindel results', headingFormat)

    with open(args.bedfile) as bed:
        #for line in bed:
        #    print(line.split(" ")[3])

        genesDup = [line.split(" ")[3] for line in bed]
        genes = set(genesDup)

    #print(genes)

    worksheet.write('A4', 'Genes included: ')
    row = 5
    for gene in genes:
        worksheet.write('A'+str(row), gene)
        row += 1

    workbook.close()

if __name__ == "__main__":
    main()
