#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import argparse
import xlsxwriter
import os

def extract_stats(statsfile, statstype, sampletype, statsdict):
    with open(statsfile, 'r') as statsfile:
        if statstype not in statsdict:
            statsdict[statstype] = {}
        statsdict[statstype][sampletype] = {}
        
        for row in statsfile:
            row = row.rstrip()
            row_list = row.split("\t")
            if 0 in statsdict[statstype][sampletype]:
                statsvalues = row_list
                for column, statsvalue in enumerate(statsvalues):
                    statsvalue = statsvalue.replace(".", ",")
                    statsdict[statstype][sampletype][column]["colvalue"] = statsvalue
                break
            if row_list[0] == "GENOME_TERRITORY" or row_list[0] == "LIBRARY":
                headernames = row_list
                headercount = 0
                for headername in headernames:
                    statsdict[statstype][sampletype][headercount] = {} 
                    statsdict[statstype][sampletype][headercount]["colname"] = headername
                    headercount += 1
        return statsdict

def create_excel(statsdict, output, normalname, tumorname):
    excelfile = xlsxwriter.Workbook(output)
    worksheet = excelfile.add_worksheet("qc_stats")
    worksheet.set_column(0, 30, 20)
    cellformat = {}
    cellformat["header"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': 'black'})
    cellformat["tumorname"] = excelfile.add_format({'bold': True, 'bg_color': 'FF7979'})
    cellformat["normalname"] = excelfile.add_format({'bold': True, 'bg_color': '94FF79'})
    row = 1
    for statstype in statsdict:
        header = False
        for sampletype in statsdict[statstype]:
            if sampletype == "tumor":
                sname = tumorname
                nameformat = cellformat["tumorname"]
            else:
                sname = normalname
                nameformat = cellformat["normalname"]
            # loop through values in stats --> sample
            column_s = 1
            if not header:
                for column in statsdict[statstype][sampletype]:
                    write_col = column_s + column
                    worksheet.write(row, 0, "Samplename", cellformat["header"])
                    worksheet.write(row, write_col, statsdict[statstype][sampletype][column]["colname"], cellformat["header"])
                row += 1
                header = True
            for column in statsdict[statstype][sampletype]:
                write_col = column_s + column
                worksheet.write(row, 0, sname, nameformat)
                worksheet.write(row, write_col, statsdict[statstype][sampletype][column]["colvalue"])
            row += 1
        row += 1
    excelfile.close()

def create_excel_main(tumorcov, normalcov, tumordedup, normaldedup, output):
    tumorcovfile = os.path.basename(tumorcov)
    tumorname = tumorcovfile.replace("_WGScov.tsv", "")
    normalcovfile = os.path.basename(normalcov)
    normalname = normalcovfile.replace("_WGScov.tsv", "")
    statsdict = {}
    statsdict = extract_stats(tumorcov, "coverage",  "tumor", statsdict)
    statsdict = extract_stats(normalcov, "coverage", "normal", statsdict)
    statsdict = extract_stats(tumordedup, "dedup",  "tumor", statsdict)
    statsdict = extract_stats(normaldedup, "dedup", "normal", statsdict)


    if not output.endswith(".xlsx"):
        output = f"{output}.xlsx"

    create_excel(statsdict, output, normalname, tumorname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tc', '--tumorcov', nargs='?', help='Sentieon WGS cov file from tumorbam', required=True)
    parser.add_argument('-nc', '--normalcov', nargs='?', help='Sentieon WGS cov file from normalbam', required=True)
    parser.add_argument('-td', '--tumordedup', nargs='?', help='Sentieon dedup-stats for tumorbam', required=True)
    parser.add_argument('-nd', '--normaldedup', nargs='?', help='Sentieon dedup-stats for normalbam', required=True)
    parser.add_argument('-o', '--output', nargs='?', help='fullpath to file to be created (xlsx will be appended if not written)', required=True)
    args = parser.parse_args()
    create_excel_main(args.tumorcov, args.normalcov, args.tumordedup, args.normaldedup, args.output)
