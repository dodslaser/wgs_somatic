# Pair check

In the pipeline, a pair check is done to make sure that tumor and normal samples come from the same patient. Script [determine_match.py](https://github.com/ClinicalGenomicsGBG/wgs_somatic/blob/master/workflows/scripts/determine_match.py) is used for this.

“SNVs Only” vcf:s for tumor and normal (outputs from DNAScope) are used to determine match for pair-check. (“{workingdir}/{stype}/dnascope/{sname}_germline_SNVsOnly.recode.vcf”)
 
Tumor vcf:
For each chromosome, starting from position 200k, dbsnp labeled variants are extracted. 20k variants are extracted from one chromosome until it switches to the next chromosome. In total, 400k dbsnp labeled variants are added to a list.
 
Normal vcf:
For each chromosome, starting from position 200k, dbsnp labeled variants are extracted. 40k variants are extracted from one chromosome until it switches to the next chromosome. In total, 800k dbsnp labeled variants are added to a list.
 
Each variant in the list of tumor variants that exists in the list of normal variants is labeled as a match and each variant that does not exist is labeled as a mismatch. If the match fraction is >=0.95, the match status between the tumor and normal samples will be “match”. If the match fraction is <0.95 but >=0.85, the match status will be “warning match lower than expected”. If the match fraction is <0.85, the match status will be “error, large mismatch”. 
