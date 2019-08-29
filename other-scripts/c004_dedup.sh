#!/bin/bash -ve
# remove PCR duplicates

# request memory for job (default 6G, max 72G)
#$ -l mem=8G
#$ -l rmem=8G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=02:00:00
# sends email with maximal memory and time used by job
# current environment settings are used for the job???
#$ -V
#$ -t 1-137

taskid=${SGE_TASK_ID}

inds=`less keep_all_saxatilis.txt` # all ind IDs
ind=`echo $inds | cut -d" " -f $taskid` # get focal ind

module add apps/java/1.7

java -Xmx6g -jar /home/bo1awx/programs/picard-tools-1.129/picard.jar MarkDuplicates INPUT=$ind-filtered.bam OUTPUT=$ind-filtered-dedup.bam METRICS_FILE=$ind.metrics REMOVE_DUPLICATES=True READ_NAME_REGEX=null ASSUME_SORTED=True
# REMOVE_DUPLICATES: remove duplicate rather than just marking them
# READ_NAME_REGEX: used to estimate amount of optical duplicates, and thereby true library size. If set to null, no optical duplicate detection.
# ASSUME_SORTED=True: file has been sorted before
