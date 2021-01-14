## Logfiles directory

Logs are added to the “logfiles” directory for each date the pipeline is started. The logs contain information about which paths from the source that are mounted to the different singularity images for each run. The logs also contains information about, for example, in case of invalid hg38 argument or errors found in other arguments used when starting the pipeline and if the fastq files are not found in the specified directory. 


Logs are also added to the subdirectory “commandlogs” for each day the pipeline is started. These logs contain the commands and arguments used for starting the pipeline. 
