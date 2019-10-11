#!/bin/sh
### this script will sort, filter and index the reads mapped against the reference genome
## usage: sbatch 08.1-plink.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=plink
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=3G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=plink.out
#SBATCH --error=plink.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL



## Add module for plink
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/plink/1.90

## Calculate heterozygosity
#plink --bfile Tms.scaf.P1 --het --out Tms.P1.het --allow-extra-chr
#plink --bfile Tms.scaf.P2 --het --out Tms.P2.het --allow-extra-chr
plink --vcf Tms.trim3.scafs.srt.vcf --out Tms.trim3 --make-bed --allow-extra-chr --recode --maf 0.2
#plink --vcf Tms.scaf.vcf --het --out Tms.all.het --allow-extra-chr

## Calculate LD interchr and intrachr
# intrachr
plink --file Tms.trim3 --r2 --allow-extra-chr --out Tms.all.intraLD --ld-window-r2 0

# interchr
plink --file Tms.trim3 --r2 inter-chr --allow-extra-chr --out Tms.all.interLD --ld-window-r2 0


## Calculate LD for LD decay plot
# 1) reduce false positive rates by using only SNPs with MAF > 0.2
# plink --file Tms.trim3.POP1 --maf 0.2 --recode --out Tms.trim3.P1.frq --allow-extra-chr --make-bed
# plink --file Tms.trim3.POP2 --maf 0.2 --recode --out Tms.trim3.P2.frq --allow-extra-chr --make-bed

## 2) calculate LD with no bottom limits filters (--ld-window-r2 0) and huge binnings (window parameters)
#plink --bfile Tms.scaf.P1.frq --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out Tms.scaf.P1.frq --allow-extra-chr
#plink --bfile Tms.scaf.P2.frq --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out Tms.scaf.P2.frq --allow-extra-chr


### repeat all for the simulated dataset
# plink --vcf asex3.vcf --out asex3 --make-bed --allow-extra-chr --recode
# plink --vcf asex3.vcf --het --out asex3.het --allow-extra-chr

## 2) calculate LD with no bottom limits filters (--ld-window-r2 0) and huge binnings (window parameters)
# plink --bfile asex3 --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out asex3 --allow-extra-chr

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/plink/1.90


