# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def calc_gender(wgscovfile, ycovfile):
    with open(wgscovfile, 'r') as wgscov:
        for linenumber, line in enumerate(wgscov):
            if linenumber == 2:
                wgscov = float(line.split("\t")[1])
    with open(ycovfile, 'r') as ycov:
        for linenumber, line in enumerate(ycov):
            if linenumber == 2:
                ycov = float(line.split("\t")[1])

    gender_ycov_threshold = 0.1
    ycov_fraction = float(ycov/wgscov)
    if ycov_fraction > gender_ycov_threshold:
        gender = "male"
    else:
        gender = "female"
    return gender
