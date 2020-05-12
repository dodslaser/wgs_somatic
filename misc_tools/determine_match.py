#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def extract_variants(vcffile, numvars):
    with open(vcffile, 'r') as vcf:
        variants = []
        countvars = 0
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
                    countvars += 1
                    variants.append(variant_info)
        return "Error" 

def determine_match(normalvcf, tumorvcf, numvariants):
    
    tumorvariants = extract_variants(tumorvcf, numvariants)
    normalvariants = extract_variants(normalvcf, numvariants)
    
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

tumorvcf = "/medstore/Development/WGS_Somatic/Snakemake_testing/WNB149AT/tumor/dnascope/WNB149AT_germline_SNVsOnly.recode.vcf"
normalvcf = "/medstore/Development/WGS_Somatic/Snakemake_testing/WNB149AT/normal/dnascope/WNB149N_germline_SNVsOnly.recode.vcf"
output = "/seqstore/webfolders/wgs/admin/barncancer/match_samples.tsv"
determine_match(normalvcf, tumorvcf, 400000)
