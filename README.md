# The WGS Somatic Pipeline

### Needs a cool name

Ideas:

Just anothEr Wholegenome aNalysis for detection of Somatic variants -- JEWNS

WOPR - junior


## General description 

 Pipeline started out as a giant shell-script developed specifically for analysing neuroblastoma WholeGenomeSequenced samples in stockholm for Tommy Martinssons research group. Contact: Susanne Fransson (susanne.fransson@clingen.gu.se)

 Pipeline was converted to Snakemake and updated in various ways, including a transition from hg19 to hg38 as a reference genome in conjunction with the project from BarncancerFonden in association with Genomic Medicine Sweden. Clinical Genomics Göteborg has had a bioinformatician employed in this project at 80% (of FTE).

 The pipeline takes as input-data fastqfiles from a TUMOR and a NORMAL sample and generates a group of result and QC-files. The results contain Somatic AND Germline variantcalls of SNVs and InDels, as well as Structural Variants (SVs) and Copy Number Variants (CNVs) -- to provide the possibility of discovering both cancer pre-disposition variants but with a primary focus towards acquired mutations. 

### How to install:

1. Clone the repository

`$ git clone https://github.com/ClinicalGenomicsGBG/wgs_somatic`

2. Install submodules

`$ git submodule update --init --recursive`

Submodules used are annotate\_manta\_canvas, b\_allele\_igv\_plot, canvas\_to\_interpreter and alissa\_connector\_upload. They can all be foundat [CGG](https://github.com/ClinicalGenomicsGBG).



You shouldn't have to build singularity images since paths to them are specified in the configs but if you would like to build them, you can use the definition files in the singularity subfolders.

Singularity images are located here: /apps/bio/singularities/wgs\_somatic/

### How to run:

1. Start a screen

2. Load anaconda2 

`$ module load anaconda2`

3. Activate conda environment

`$ source activate wgs_somatic`

4. Run launch\_snakemake.py

Depending on whether you have GMS samples or WNB (Susanne's samples), you start the pipeline either by:

* GMS Sample

```
$ ./launch_snakemake.py \  
    --runnormal <ID of sequencing run, e.g. 200616_A00689_0144_BHFNK2DSXY> \  
    --runtumor <ID of sequencing run, e.g. 200616_A00689_0144_BHFNK2DSXY> \  
    --outputdir <Output directory, e.g. /seqstore/webfolders/wgs/barncancer/hg38/DNA66246> \  
    --normalsample <Name of normal sample, e.g. DNA66245> \  
    --normalfastqs <path to directory containing normal fastqs, e.g. /seqstore/hcp_backed_up_data/local_wgs/barncancer/DNA66245> \  
    --tumorsample <Name of tumor sample, e.g. DNA66246> \  
    --tumorfastqs <path to directory containing tumor fastqs, e.g. /seqstore/hcp_backed_up_data/local_wgs/barncancer/DNA66246> \  
    --igvuser barncancer_hg38 \  
    --hg38ref yes
``` 

* WNB (Susanne's samples)

```
$ ./launch_snakemake.py \  
    --runnormal <ID of sequencing run> \  
    --runtumor <ID of sequencing run> \  
    --outputdir <Output directory, e.g. /seqstore/webfolders/wgs/neuroblastom/<tumor-sample-name> \  
    --normalsample <Name of normal sample> \  
    --normalfastqs <path to directory containing normal fastqs> \  
    --tumorsample <Name of tumor sample> \  
    --tumorfastqs <path to directory containing tumor fastqs> \  
    --igvuser susanne.fransson \
```


For hg38ref, write 'yes' if you want this option. If you want to use hg19, simply don't use hg38ref argument.




 Need to be on medair because python environment is hard-coded in launch_snakemake.py, also some tools in the snakemake workflow is hardcoded to install-locations on medair, and finally some dependencies such as references and databases are on medair.

 > \#!/apps/bio/software/anaconda2/envs/wgs_somatic/bin/python




 Or you can just use an already set-up repository! Such as here:
 `/apps/bio/repos/wgs_somatic`

 Use "launch_snakemake.py" like this:

 `./launch_snakemake.py --runnormal 200616_A00689_0144_BHFNK2DSXY --runtumor 200616_A00689_0144_BHFNK2DSXY --outputdir /seqstore/webfolders/wgs/barncancer/hg38/DNA66246 --normalsample DNA66245 --normalfastqs /seqstore/hcp_backed_up_data/local_wgs/barncancer/DNA66245 --tumorsample DNA66246 --tumorfastqs /seqstore/hcp_backed_up_data/local_wgs/barncancer/DNA66246 --igvuser barncancer_hg38 --hg38ref yes`

 Runnormal and runtumor is only used to create a unique samplename based on the sequencing run the data comes from. Could probably be done in a better way.

 ### Goal:


 * Have a high quality analysis of all types of genetic abberations (starting with SNVs and InDels, CNVs and SVs). 
 * Which is fast (because some of these samples can help treatment of urgent pediatric cases). 
 * Automated and connected to hospital systems so that it is not relied upon bioinformatician working-hours.
 * Packaged into Singularities so that it can in the future be run on a cloud-platform like AWS.
 * Connected to HCP to upload results (and download data for analysis?).
 * Robust and well-documented to fit into the clinical requirements.


 ### Status (2020-10-19):

 #### Reference genome

 Currently operational on hg38 reference genome. Maybe hg19 still works but it has not been tested since the hg38 transition, possible that something has been broken. But the plan is to run this pipeline on hg38 in the immediate future, so maybe compatability with hg19 can be ensured later, or upon request.

 #### Validation and optimisation:

 ##### Germline calling:
 * SNVs and InDel calling with DNAscope (was validated through WOPR in hg19, but not on hg38)
 * CNVs -- Canvas 
 * SVs -- Manta


 ##### Somatic calling:
 * SNVs and InDel calling with TNscope and custom Bcftools filters. Optimised against artificial truthset by mixing coriell samples (see wiki)
 * CNVs -- Canvas
 * SVs -- Manta

 #### Code structure and organisation of dependencies

 Most of the pipeline is packaged in Singularities, Sentieon and Canvas (Sentieon is majority of pipeline).
 The remaining rules relies primarily on small python scripts, bcftools and Manta. Should not be too difficult to package these into singularities as well.

 #### Validation project GMS - BT

 * Analyse all samples in the project with hg38 (nearly done)
 * Calculate sensitivity and precision with artificial truthset (nearly done)

### Yearly statistics

After running the pipeline for the first time, a yearly\_statistics text file is created in the repo. Every time the pipeline is run, sample name and date/time is added to this text file.
