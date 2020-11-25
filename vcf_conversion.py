#vcf_file = open("/seqstore/webfolders/wgs/barncancer/hg38/DNA68517/DNA68517_200904_AHFKCKDSXY_somatic_refseq3kfilt.vcf.gz", "r")

#print(vcf_file.read())


#import vcf
#vcf_reader = vcf.Reader(open('/home/xshang/vcf_test/DNA66246_200616_BHFNK2DSXY_somatic_refseq3kfilt.vcf', 'r'))
#vcf_writer = vcf.Writer(open('/home/xshang/vcf_test/DNA66246_200616_BHFNK2DSXY_somatic_refseq3kfilt_modified.vcf', 'w'), vcf_reader)

#for record in vcf_reader:
#    record = record.FORMAT.replace("AFDP", "DP")
    #print(record)
#    vcf_writer.write_record(record)


#read the input file
vcf_file = open("/home/xshang/vcf_test/DNA65866_200622_BH77WLDSXY_somatic_refseq3kfilt_Alissa.vcf", "rt")
#read file contents to string
data = vcf_file.read()
#replace all occurrences of AFDP with DP
data = data.replace('AFDP', 'DP')
#close the input file
vcf_file.close()
#open the input file in write mode
vcf_file = open("/home/xshang/vcf_test/DNA65866_200622_BH77WLDSXY_somatic_refseq3kfilt_Alissa.vcf", "wt")
#overwrite the input file with the resulting data
vcf_file.write(data)
#close the file
vcf_file.close()