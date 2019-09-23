### this script will trim the Tms reads ###
# usage: bsub < ./04-trimming-Tms.sh

#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J my_job_name
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000000
#BSUB -o trim2.out
#BSUB -e trim2.err
##BSUB -n 4 ## 4 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host

module load UHTS/Quality_control/cutadapt/2.3

for file in *.fq.gz; do
        cutadapt -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 80 -o $file-trim3.fq.gz $file;
done

module rm UHTS/Quality_control/cutadapt/2.3

# To trim a 3’ adapter, the basic command-line for Cutadapt is:
# cutadapt -a AACCGGTT -o output.fastq input.fastq


# the original Illumina adapter is : AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC.
# However, I got this message for ALL reads in the output file of the 
# first Cutadapt run: 'The adapter is preceded by "A" extremely often'
# To fix the problem I added "A" to the beginning of the adapter sequence and re-ran Cutadapt.
# Also, I eliminated all sequences smaller than 80 bp: -m 80
