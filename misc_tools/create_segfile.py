#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8
import os
import math
def create_seg(rundir, samplename, vartype):
    if os.path.isfile(f"{rundir}/CNV.CoverageAndVariantFrequency.txt"):
        cnv_covandvarfreq = f"{rundir}/CNV.CoverageAndVariantFrequency.txt"
    elif os.path.isfile(f"{rundir}/TempCNV_{samplename}/CNV.CoverageAndVariantFrequency.txt"):
        cnv_covandvarfreq = f"{rundir}/TempCNV_{samplename}/CNV.CoverageAndVariantFrequency.txt"
    else:
        try:
            cnv_covandvarfreq = [n for n in glob.glob(f"{rundir}/*/CNV.CoverageAndVariantFrequency.txt") if os.path.isfile(n)][0]
        except:
            logging.info("Could not find CNV.CoverageAndVariantFrequency.txt, no segfiles will be produced. Crashing out")

    with open(cnv_covandvarfreq, "r") as INFILE:
        with open(f"{rundir}/{samplename}_{vartype}_CNV_observed.seg", "w+") as OUTFILE:
            OUTFILE.write("#track graphType=points maxHeightPixels=300:300:300 color=0,0,0 altColor=0,0,0\n")
            OUTFILE.write("Sample\tChromosome\tStart\tEnd\tCNV_Observed\n")
            with open(f"{rundir}/{samplename}_{vartype}_CNV_called.seg", "w+") as OUTFILE2:
                OUTFILE2.write("#track graphType=points maxHeightPixels=300:300:300 color=0,0,220 altColor=220,0,0\n")
                OUTFILE2.write("Sample\tChromosome\tStart\tEnd\tCNV_Called\n")
                for line in INFILE:
                    array_2 = line.split("\t")
                    length = len(array_2)
                    cnv = 0
                    ncov = 0

                    try:
                        cnv = float(array_2[3])
                        ncov = float(array_2[6])
                    except (IndexError, ValueError) as e:
                        continue

                    if ncov > 0 and cnv > 0:
                        cnvlog = math.log(cnv, 2)
                        covlog = math.log(ncov, 2)
                        if not array_2[0] == "X" or not array_2[0] == "Y":
                            cnvlog -= 1
                            covlog -= 1

                        OUTFILE.write("Observed_CNVs\t%s\t%s\t%s\t%s\n" % (array_2[0], array_2[1], array_2[2], covlog))
                        OUTFILE2.write("Called_CNVs\t%s\t%s\t%s\t%s\n" % (array_2[0], array_2[1], array_2[2], cnvlog))
        return [f"{rundir}/{samplename}_{vartype}_CNV_observed.seg", f"{rundir}/{samplename}_{vartype}_CNV_called.seg"]


#create_seg("/apps/bio/repos/wgs_pipeline_somatic/test_canvas/DNA64925/canvasresult_normalgermline", "DNA64926")
