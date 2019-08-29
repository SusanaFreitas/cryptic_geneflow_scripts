#### this script will produce a read depth / position table using picard-tools and GATK
## usage: bsub < ./07.4-seeDP.sh

#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J fb.trim3
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o gatk-getDP.out
#BSUB -e gatk-getDP.err
###BSUB -n 8 ## 12 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host

module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
module load UHTS/Analysis/picard-tools/2.9.0

## first we need to create a dictionary for the fasta reference with Picard

picard-tools CreateSequenceDictionary REFERENCE=3_Tms_b3v08.fasta OUTPUT=3_Tms_b3v08.dict

## then we will create a DP table with GATK

GenomeAnalysisTK VariantsToTable \
 -R 3_Tms_b3v08.fasta \
 -V Tms.fb.vcf \
 -F CHROM -F POS -GF GT -GF DP \
 --output Tms.fb.DP.table

module rm UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
module rm UHTS/Analysis/picard-tools/2.9.0
