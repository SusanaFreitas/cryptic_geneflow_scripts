#!/bin/bash -ve
# map reads to genome
# IDs of all individuals are in keep_all_saxatilis.txt
# for each individual, there is a separate (parallel) task

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=06:00:00
#$ -t 1-137

taskid=${SGE_TASK_ID}

inds=`less keep_all_saxatilis.txt` # all ind IDs
ind=`echo $inds | cut -d" " -f $taskid` # get focal ind

filename=$ind.fastq-trimmed3.fastq # filename - like ANG_C_05.fastq-trimmed3.fastq
																				   
# now map the reads
/home/bo1awx/programs/bwa-0.7.12/bwa mem -M -R "@RG\tID:$ind\tSM:$ind" littorina_genome.fasta $filename > $ind.sam
# -M: Mark shorter split hits as secondary (for Picard compatibility)
# -R: read group
