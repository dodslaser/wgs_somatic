#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import argparse
import xlsxwriter
import os
from workflows.scripts.determine_match import determine_match
from tools.git_versions import get_git_commit, get_git_tag, get_git_reponame
from workflows.scripts.sex import calc_sex
import time
import pandas as pd
import glob


def add_insilico_stats(insilicofolder, main_excel):
    from workflows.scripts.insilico_coverage import insilico_overall_coverage
    
    # generate overall insilico coverage dataframe & append to main excel
    ocov_file_list = glob.glob(f"{insilicofolder}/**/*_cov.tsv")
    print(f"ocov_file_list: {ocov_file_list}")
    for ocov in ocov_file_list:
        print(f"ocov: {ocov}")
        ocov_sheetname = os.path.basename(os.path.splitext(ocov)[0])
        ocov_list = insilico_overall_coverage.overall_coverage_stats(ocov, "10,20,30,40,50,60,70,80,90,100,110,120")
        ocov_df = pd.DataFrame(ocov_list)
        with pd.ExcelWriter(main_excel, engine='openpyxl', mode='a') as writer:
            ocov_df.to_excel(writer, sheet_name=ocov_sheetname, index=False, header=False)

    # generate per region stats and append to main excel
    excel_file_list = glob.glob(f"{insilicofolder}/**/*.xlsx")
    print(f"excel_file_list: {excel_file_list}")
    for xlfile in excel_file_list:
        dataframe = pd.read_excel(xlfile, engine='openpyxl')
        dataframe = dataframe.iloc[: , 1:]
        sheetname = os.path.basename(os.path.splitext(xlfile)[0])
        print(f"sheetname: {sheetname}")
        with pd.ExcelWriter(main_excel, engine='openpyxl', mode='a') as writer:
            dataframe.to_excel(writer, sheet_name=sheetname, index=False)


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


def get_canvas_tumorinfo(canvasvcf):
    canvasdict = {}
    canvas_infofields = ["##OverallPloidy", "##DiploidCoverage", "##EstimatedTumorPurity", "##PurityModelFit", "##InterModelDistance", "##LocalSDmetric", "##EvennessScore", "##HeterogeneityProportion", "##EstimatedChromosomeCount"]
    with open(canvasvcf, 'r') as vcf:
        for variant in vcf:
            variant = variant.rstrip('\n')
            variant_info = variant.split('\t')
            if variant_info[0].startswith('#'):
                if variant_info[0].split("=")[0] in canvas_infofields:
                    canvasfield = variant_info[0].split("=")[0]
                    canvasfield = canvasfield.replace("#", "")
                    canvasfield_value = variant_info[0].split("=")[1]
                    canvasfield_value = canvasfield_value.replace(".", ",")
                    canvasdict[canvasfield] = canvasfield_value
    return canvasdict


def create_excel(statsdict, output, normalname, tumorname, match_dict, canvasdict, sex):
    current_date = time.strftime("%Y-%m-%d")
    excelfile = xlsxwriter.Workbook(output)
    worksheet = excelfile.add_worksheet("qc_stats")
    worksheet.set_column(0, 30, 20)
    cellformat = {}
    cellformat["header"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': 'black'})
    cellformat["section"] = excelfile.add_format({'bold': True, 'font_color': 'FFA764', 'bg_color': '878787'})
    cellformat["tumorname"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': 'BF000D'})
    cellformat["normalname"] = excelfile.add_format({'bold': True, 'font_color': 'white', 'bg_color': '00BF19'})
    cellformat["warning"] = excelfile.add_format({'bg_color': 'FF9A00'})
    cellformat["error"] = excelfile.add_format({'bg_color': 'FF0000'})
    cellformat["pass"] = excelfile.add_format({'bg_color': '95FF80'})

    row = 1
    worksheet.write(row, 0, f"QC-report created: {current_date}")
    row += 1
    worksheet.write(row, 0, f"{get_git_reponame()} tag: {get_git_tag()}, commit: {get_git_commit()}")
    worksheet.write(row, 4, f"Computed sex of patient: {sex}")
    row += 2
 
    for statstype in statsdict:
        header = False
        if statstype == "coverage":
            worksheet.write(row, 0, "COVERAGE STATS", cellformat["section"])
        else:
            worksheet.write(row, 0, "MAPPING STATS", cellformat["section"])
        row += 1
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
    row += 1
    worksheet.write(row, 0, "TUMOR/NORMAL MATCH-CHECK", cellformat["section"])
    row += 1
    if match_dict:
        if "ERROR" in match_dict["match_status"]:
            worksheet.write(row, 0, "Error occurred", cellformat["error"])
        else:
            headerrow = row
            statrow = row + 1
            column_num = 0
            for stat in match_dict:
                worksheet.write(headerrow, column_num, stat, cellformat["header"])
                if stat == "match_status" or stat == "match_fraction":
                    if "warning" in match_dict["match_status"]:
                        style = "warning"
                    elif "error" in match_dict["match_status"]:
                        style = "error"
                    else:
                        style = "pass"
                    worksheet.write(statrow, column_num, match_dict[stat], cellformat[style])
                else:
                    worksheet.write(statrow, column_num, match_dict[stat])
                
                column_num += 1
    row += 4
    worksheet.write(row, 0, "CANVAS-STATS", cellformat["section"])
    worksheet.write(row, 1, tumorname, cellformat["tumorname"])
    row += 1
    if canvasdict:
        for key in canvasdict:
            worksheet.write(row, 0, key, cellformat["header"])
            worksheet.write(row, 1, canvasdict[key])
            row += 1

    excelfile.close()


def create_excel_main(tumorcov='', ycov='', normalcov='', tumordedup='', normaldedup='', tumorvcf='', normalvcf='', canvasvcf='', output='', insilicodir=''):
    print(f"insilicodir: {insilicodir}")
    statsdict = {}
    if tumorcov:
        tumorcovfile = os.path.basename(tumorcov)
        tumorname = tumorcovfile.replace("_WGScov.tsv", "")
        statsdict = extract_stats(tumorcov, "coverage",  "tumor", statsdict)
        statsdict = extract_stats(tumordedup, "dedup",  "tumor", statsdict)
    if normalcov:
        normalcovfile = os.path.basename(normalcov)
        normalname = normalcovfile.replace("_WGScov.tsv", "")
        statsdict = extract_stats(normalcov, "coverage", "normal", statsdict)
        statsdict = extract_stats(normaldedup, "dedup", "normal", statsdict)
    if tumorcov and normalcov:
        #if canvasvcf:
        #    canvas_dict = get_canvas_tumorinfo(canvasvcf)
        match_dict = determine_match(normalvcf, tumorvcf, 400000)
        canvas_dict = get_canvas_tumorinfo(canvasvcf)

    

    if not output.endswith(".xlsx"):
        output = f"{output}.xlsx"

    # Determine which files to use to calculate sex and where to get insilico coverageg files
    if tumorcov:
        if normalcov:
            # Tumour + Normal
            calculated_sex = calc_sex(normalcov, ycov)
            create_excel(statsdict, output, normalname, tumorname, match_dict, canvas_dict, sex=calculated_sex)
            add_insilico_stats(insilicodir, output)
        else:
            # Tumour only
            calculated_sex = calc_sex(tumorcov, ycov)
            create_excel(statsdict, output, normalname='', tumorname=tumorname, match_dict='', canvasdict='', sex=calculated_sex)
            add_insilico_stats(insilicodir, output) # Maybe this can be commented out if not needed for tumour only
    else:
        # Normal only
        calculated_sex = calc_sex(normalcov, ycov)
        create_excel(statsdict, output, normalname, tumorname='', match_dict='', canvasdict='', sex=calculated_sex)
        add_insilico_stats(insilicodir, output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tc', '--tumorcov', nargs='?', help='Sentieon WGS cov file from tumorbam', required=False)
    parser.add_argument('-yc', '--ycov', nargs='?', help='Sentieon Y-cov file', required=False)
    parser.add_argument('-nc', '--normalcov', nargs='?', help='Sentieon WGS cov file from normalbam', required=True)
    parser.add_argument('-td', '--tumordedup', nargs='?', help='Sentieon dedup-stats for tumorbam', required=False)
    parser.add_argument('-nd', '--normaldedup', nargs='?', help='Sentieon dedup-stats for normalbam', required=True)
    parser.add_argument('-tv', '--tumorvcf', nargs='?', help='Tumor Germlinecalls', required=False)
    parser.add_argument('-nv', '--normalvcf', nargs='?', help='Normal Germlinecalls', required=True)
    parser.add_argument('-cv', '--canvasvcf', nargs='?', help='Somatic Canvas VCF', required=False)
    parser.add_argument('-is', '--insilicodir', nargs='?', help='Full path to insilico directory (which contains excel files)', required=False)
    parser.add_argument('-o', '--output', nargs='?', help='fullpath to file to be created (xlsx will be appended if not written)', required=True)
    args = parser.parse_args()
    create_excel_main(args.tumorcov, args.ycov, args.normalcov, args.tumordedup, args.normaldedup, args.tumorvcf, args.normalvcf, args.canvasvcf, args.output, args.insilicodir)
