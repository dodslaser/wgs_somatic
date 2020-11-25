import sys

#read the input file
vcf_file = open(sys.argv[1], "rt")
#read file contents to string
data = vcf_file.read()
#replace all occurrences of AFDP with DP
data = data.replace('AFDP', 'DP')
#close the input file
vcf_file.close()
#open the input file in write mode
vcf_file = open(sys.argv[1], "wt")
#overwrite the input file with the resulting data
vcf_file.write(data)
#close the file
vcf_file.close()