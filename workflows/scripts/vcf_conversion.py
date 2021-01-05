#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8



def vcf_conversion(vcf):
    vcf_file = open(vcf, "rt")
    data = vcf_file.read()
    data = data.replace('AFDP', 'DP')
    vcf_file.close()
    vcf_file = open(vcf, "wt")
    vcf_file.write(data)
    vcf_file.close()


# Old code
# import sys

#read the input file
# vcf_file = open(sys.argv[1], "rt")
#read file contents to string
# data = vcf_file.read()
#replace all occurrences of AFDP with DP
# data = data.replace('AFDP', 'DP')
#close the input file
# vcf_file.close()
#open the input file in write mode
# vcf_file = open(sys.argv[1], "wt")
#overwrite the input file with the resulting data
# vcf_file.write(data)
#close the file
# vcf_file.close()
