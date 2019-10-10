#!/bin/sh
### script to filter by genotype DP (read depth)
## Exclude all genotypes bellow 8 reads and above 200 reads coverage

## I had a problem with this script and I was getting constant 'invalid argument' errors.
## I solved it by removing all spaces in between the formal arguments...
# http://selectvariants1.rssing.com/chan-9504484/all_p2.html

### this script will filter by  DP using GATK
## usage: sbatch 07.5-filterDP.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=DPfilt
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=3G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=DPfilt.out
#SBATCH --error=DPfilt.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL


## add module for GATK
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0

## fb file (without trimming reads)
GenomeAnalysisTK VariantFiltration\
        -R 3_Tms_b3v08.fasta\
        -V Tms.fb.vcf --genotype-filter-expression "DP<8||DP>200"\
        --genotype-filter-name "DP_8-200"\
        --set-filtered-genotype-to-no-call\
        --output Tms.fb_DPfilt.vcf

## trim2 reads
GenomeAnalysisTK VariantFiltration\
        -R 3_Tms_b3v08.fasta\
        -V Tms.trim2.fb.vcf --genotype-filter-expression "DP<8||DP>200"\
        --genotype-filter-name "DP_8-200"\
        --set-filtered-genotype-to-no-call\
        --output Tms.fb_DPfilt.vcf

## trim3 reads
GenomeAnalysisTK VariantFiltration\
        -R 3_Tms_b3v08.fasta\
        -V Tms.trim3.fb.vcf --genotype-filter-expression "DP<8||DP>200"\
        --genotype-filter-name "DP_8-200"\
        --set-filtered-genotype-to-no-call\
        --output Tms.fb_DPfilt.vcf

## OPTIONS:
## --set-filtered-genotype-to-no-call : used this to set the filtered genotypes (the ones that less than 8 > DP > 200)
## to ./.
## https://gatkforums.broadinstitute.org/gatk/discussion/7577/removing-variants-based-on-the-format-field


module rm UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
module rm Bioinformatics/Software/vital-it
