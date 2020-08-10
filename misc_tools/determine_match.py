#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def extract_variants(vcffile, numvars, chromlist):
    sections = 20
    chrom_switch = int(numvars/sections)
    with open(vcffile, 'r') as vcf:
        variants = []
        countvars = 0
        count_to_switch = 0
        chrom_switches = 0
        for row in vcf:
            if countvars == numvars:
                vardict = {}
                for variant in variants:
                    variantpos = "_".join(variant[0:2])
                    variantalleles = "_".join(variant[3:5])
                    variantname = f"{variantpos}_{variantalleles}"
                    vardict[variantname] = variant[9]
                return vardict
            row = row.rstrip()
            variant_info = row.split("\t")
            if not variant_info[0].startswith('#'):
                if not variant_info[2] == ".":
                    chrom = variant_info[0]
                    pos = variant_info[1]
                    if int(pos) < 200000:
                        continue
                    if chrom == chromlist[chrom_switches]:
                        if count_to_switch == chrom_switch:
                            chrom_switches += 1
                            count_to_switch = 0
                        countvars += 1
                        count_to_switch += 1
                        variants.append(variant_info)
        return "Error" 

def extract_chromosomes(vcffile):
    chrom_list = []
    with open(vcffile, 'r') as vcf:
        for row in vcf:
            if row.startswith("#"):
                if row.startswith("##contig"):
                    contiginfo = row.split(",")
                    chromname = contiginfo[0].split("=")[2]
                    chrom_list.append(chromname)
            else:
                return chrom_list

def determine_match(normalvcf, tumorvcf, numvariants):
    chrom_list = extract_chromosomes(normalvcf) 
    tumorvariants = extract_variants(tumorvcf, numvariants, chrom_list)
    normalvariants = extract_variants(normalvcf, 800000, chrom_list)
    
    if "Error" in tumorvariants or "Error" in normalvariants:
        match_dict = {}
        match_dict["match_status"] = "ERROR"
        return match_dict    

    count_matches = 0
    count_mismatches = 0
    for variant in tumorvariants:
        if variant in normalvariants:
            count_matches += 1
        else:
            count_mismatches += 1

    match_fraction = float(count_matches/numvariants)
    
    match_dict = {}
    match_dict["dbsnp_positions"] = numvariants
    match_dict["matches"] = count_matches
    match_dict["mismatches"] = count_mismatches    
    match_dict["match_fraction"] = match_fraction
    match_dict["match_status"] = ""

    if match_fraction >= float(0.95):
        match_dict["match_status"] = "match"
    elif match_fraction >= float(0.85):
        match_dict["match_status"] = "warning match lower than expected"
    else:
        match_dict["match_status"] = "error, large mismatch"

    return match_dict

#normalvcf = "/seqstore/webfolders/wgs/neuroblastom/WNB125AT/normal/dnascope/WNB125N_germline_SNVsOnly.recode.vcf"
#tumorvcf = "/seqstore/webfolders/wgs/neuroblastom/WNB125AT/tumor/dnascope/WNB125AT_germline_SNVsOnly.recode.vcf"
#output = "/seqstore/webfolders/wgs/admin/barncancer/match_samples.tsv"
#print(determine_match(normalvcf, tumorvcf, 400000))
#normalvcf = "/seqstore/webfolders/wgs/neuroblastom/WNB124AT/normal/dnascope/WNB124N_germline_SNVsOnly.recode.vcf"
#tumorvcf = "/seqstore/webfolders/wgs/neuroblastom/WNB125AT/normal/dnascope/WNB125N_germline_SNVsOnly.recode.vcf"
#print(determine_match(normalvcf, tumorvcf, 400000))
#normalvcf = "/medstore/results/wgs/KK/200601_A00687_0082_BH7GKNDSXY/DNA66348_200601_BH7GKNDSXY/dnascope/DNA66348_200601_BH7GKNDSXY_germline_pass_SNVs.vcf.recode.vcf"
#tumorvf = "/medstore/results/wgs/KK/200601_A00687_0082_BH7GKNDSXY/DNA66155_200601_BH7GKNDSXY/dnascope/DNA66155_200601_BH7GKNDSXY_germline_pass_SNVs.vcf.recode.vcf"
#print(determine_match(normalvcf, tumorvcf, 400000))
#normalvcf = "/medstore/results/wgs/KK/200601_A00687_0082_BH7GKNDSXY/DNA66349_200601_BH7GKNDSXY/dnascope/DNA66349_200601_BH7GKNDSXY_germline_pass_SNVs.vcf.recode.vcf"
#print(determine_match(normalvcf, tumorvcf, 400000))
