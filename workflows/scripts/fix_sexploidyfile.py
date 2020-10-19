import re


def mod_sex_vcf(sex_vcf, samplename, done_vcf_path):
    my_new_file = done_vcf_path + samplename.split("_")[0] + "_" + sex_vcf.split("/")[-1]
    to_write = []
    with open(sex_vcf, "r") as vcf:
        for line in vcf:
            split_line = re.split(r"\t+", line)
            if "SAMPLE\n" in split_line:
                replace_with = samplename.split("_")[0]
                replaced = [w.replace('SAMPLE\n', replace_with + "\n") for w in split_line]
                to_write.append(replaced)
            else:
                to_write.append(split_line)
    with open(my_new_file, "w") as result:
        for each_line_as_list in to_write:
            joined_line = "\t".join(each_line_as_list)
            result.write(joined_line)
    return my_new_file
