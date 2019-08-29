#!/bin/sh
### script to intersect a vcf file with a bed using vcflib: vcfintersect
## usage: sbatch 07.8-CombineGVCFs.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=normal  ## options: normal, long
#
#SBATCH --job-name=GVCFcombine
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=60G ## memory used per cpu (or core) in Mb
#SBATCH --time=12:00:00
# outputs
#SBATCH --output=bedintsect.out
#SBATCH --error=bedintsect.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL


### script to combine two vcf files and evaluate SNPs overlap
module add Bioinformatics/Software/vital-it

module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0

GenomeAnalysisTK CombineGVCFs \
   -R 3_Tms_b3v08.fasta \
   --variant stacks.snps.vcf \
   --variant Tms.scaf.vcf \
   -O output_combined.vcf

module rm UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
