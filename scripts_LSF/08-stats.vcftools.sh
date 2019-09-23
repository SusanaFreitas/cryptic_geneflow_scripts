#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J vcftools
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o vcftools.out
#BSUB -e vcftools.err
###BSUB -n 8 ## 12 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host

# Run vcftools to produce R^2 LD measures

## upload vcftools module
module load UHTS/Analysis/vcftools/0.1.15
#### within chromosomes
## from phased haplotypes
# vcftools --vcf populations.haps.vcf --out Tms.R2.hap --hap-r2
## from SNPs
vcftools --vcf populations.snps.vcf --out Tms.R2.snp --geno-r2

#### between chromosomes
## from phased haplotypes
vcftools --vcf populations.haps.vcf --out Tms.crR2.hap --interchrom-hap-r2
## from SNPs
vcftools --vcf populations.snps.vcf --out Tms.crR2.snp --interchrom-geno-r2

### To analyse the R2 values, we need the allele frequency
vcftools --vcf populations.snps.vcf --freq --out alf.Tms.snp
vcftools --vcf populations.haps.vcf --freq --out alf.Tms.hap

module rm UHTS/Analysis/vcftools/0.1.15
