#!/bin/bash

### this script will call SNPs using FreeBayes ###
## usage: sbatch 07.1-freebayes-trim3.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=fb.trim3
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=8 ## no of cores used per thread
#SBATCH --mem-per-cpu=4G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=fbtrim3.out
#SBATCH --error=fbtrim3.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL

module add Bioinformatics/Software/vital-it

module load UHTS/Analysis/freebayes/1.2.0
module load UHTS/Analysis/bamaddrg/2012.05.26

bamaddrg -b Gepop1.10-filtered3.bam -s pop1 \
        -b Gepop1.11-filtered3.bam -s pop1 \
        -b Gepop1.12-filtered3.bam -s pop1 \
        -b Gepop1.13-filtered3.bam -s pop1 \
        -b Gepop1.14-filtered3.bam -s pop1 \
        -b Gepop1.15-filtered3.bam -s pop1 \
        -b Gepop1.16-filtered3.bam -s pop1 \
        -b Gepop1.17-filtered3.bam -s pop1 \
        -b Gepop1.18-filtered3.bam -s pop1 \
        -b Gepop1.19-filtered3.bam -s pop1 \
        -b Gepop1.1-filtered3.bam -s pop1 \
        -b Gepop1.20-filtered3.bam -s pop1 \
        -b Gepop1.21-filtered3.bam -s pop1 \
        -b Gepop1.22-filtered3.bam -s pop1 \
        -b Gepop1.23-filtered3.bam -s pop1 \
        -b Gepop1.2-filtered3.bam -s pop1 \
        -b Gepop1.3-filtered3.bam -s pop1 \
        -b Gepop1.4-filtered3.bam -s pop1 \
        -b Gepop1.5-filtered3.bam -s pop1 \
        -b Gepop1.6-filtered3.bam -s pop1 \
        -b Gepop1.7-filtered3.bam -s pop1 \
        -b Gepop1.8-filtered3.bam -s pop1 \
        -b Gepop1.9-filtered3.bam -s pop1 \
        -b Gepop2.10-filtered3.bam -s pop2 \
        -b Gepop2.11-filtered3.bam -s pop2 \
        -b Gepop2.12-filtered3.bam -s pop2 \
        -b Gepop2.13-filtered3.bam -s pop2 \
        -b Gepop2.14-filtered3.bam -s pop2 \
        -b Gepop2.15-filtered3.bam -s pop2 \
        -b Gepop2.16-filtered3.bam -s pop2 \
        -b Gepop2.17-filtered3.bam -s pop2 \
        -b Gepop2.18-filtered3.bam -s pop2 \
        -b Gepop2.19-filtered3.bam -s pop2 \
        -b Gepop2.1-filtered3.bam -s pop2 \
        -b Gepop2.20-filtered3.bam -s pop2 \
        -b Gepop2.21-filtered3.bam -s pop2 \
        -b Gepop2.22-filtered3.bam -s pop2 \
        -b Gepop2.23-filtered3.bam -s pop2 \
        -b Gepop2.2-filtered3.bam -s pop2 \
        -b Gepop2.3-filtered3.bam -s pop2 \
        -b Gepop2.4-filtered3.bam -s pop2 \
        -b Gepop2.5-filtered3.bam -s pop2 \
        -b Gepop2.6-filtered3.bam -s pop2 \
        -b Gepop2.7-filtered3.bam -s pop2 \
        -b Gepop2.8-filtered3.bam -s pop2 \
        -b Gepop2.9-filtered3.bam -s pop2 \
                
                | freebayes -f /scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tge/5_Tge_b3v08.fasta --bam-list bamlist-trim3 \
                        --min-mapping-quality 30 --min-coverage 5 > Tge.trim3.fb.vcf


module rm UHTS/Analysis/freebayes/1.2.0
module rm UHTS/Analysis/bamaddrg/2012.05.26
