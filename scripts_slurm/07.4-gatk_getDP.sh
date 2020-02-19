#!/bin/bash

### this script will give me the DP
## usage: sbatch 07.4-gatk_getDP.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=getDP
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=3G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=getDP.out
#SBATCH --error=getDP.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL


module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0

## first we need to create a dictionary for the fasta reference with Picard

module load UHTS/Analysis/picard-tools/2.9.0
picard-tools CreateSequenceDictionary REFERENCE=1_Tdi_b3v08.fasta OUTPUT=1_Tdi_b3v08.dict

GenomeAnalysisTK VariantsToTable \
 -R 1_Tdi_b3v08.fasta \
 -V Tdi.fb.vcf \
 -F CHROM -F POS -GF GT -GF DP \
 --output Tdi.fb.DP.table

module rm UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
module rm UHTS/Analysis/picard-tools/2.9.0

