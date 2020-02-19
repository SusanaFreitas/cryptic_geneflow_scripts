#!/bin/bash

### this script will call SNPs using FreeBayes ###
## usage: sbatch 07.1-freebayes-trim2.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=fb.resort
## memory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
## outputs
#SBATCH --output=fbtrim3.out
#SBATCH --error=fbtrim3.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/freebayes/1.2.0

ulimit -s 81920

freebayes -f 3_Tms_b3v08.fasta --bam-list resortbam --min-mapping-quality 30 --min-coverage 5 > Tms.resort.vcf
#/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/Softwares/freebayes/freebayes/bin/freebayes -f 3_Tms_b3v08.fasta --bam-list trim3bam --min-mapping-quality 30 --min-coverage 5 --debug > Tms.trim3.fb.vcf

module rm UHTS/Analysis/freebayes/1.2.0

