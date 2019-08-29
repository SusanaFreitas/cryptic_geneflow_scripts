### script to filter by DP using GATK

#!/bin/bash
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J fb.trim3
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o alprim.out
#BSUB -e alprim.err
###BSUB -n 8 ## 12 CPU cores (processors) requested
####BSUB <E2><80><93>R "span[ptile=2]" ## 4 cores on the same host

#/scratch/temporary/sfreitas/softwares/vcflib/bin/vcfallelicprimitives -kg Tms.fb_DPfilter.vcf > Tge.fb_posvcflib.vcf
#/scratch/temporary/sfreitas/softwares/vcflib/bin/vcfallelicprimitives -kg Tms.trim2.fb_DPfilter.vcf > Tge.trim2.fb_posvcflib.vcf
#/scratch/temporary/sfreitas/softwares/vcflib/bin/vcfallelicprimitives -kg Tms.trim3.fb_DPfilter.vcf > Tge.trim3.fb_posvcflib.vcf

## path to vcfallelicprimitives
alprim="/scratch/temporary/sfreitas/softwares/vcflib/bin/vcfallelicprimitives"

## command line
eval $alprim -kg Tms.QUALpy.vcf > Tms.fb_posvcflib.vcf
eval $alprim -kg Tms.trim2.QUALpy.vcf > Tms.trim2_posvcflib.vcf
eval $alprim -kg Tms.trim3.QUALpy.vcf > Tms.trim3_posvcflib.vcf