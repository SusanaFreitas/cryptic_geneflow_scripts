#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q normal
#BSUB -J vcftools
#BSUB -R "rusage[mem=12000]"
##BSUB -M 12000000
#BSUB -o vcftools.out
#BSUB -e vcftools.err
###BSUB -n 8 ## 12 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host

### get column with R2 between crs
# awk '{ print $6 }' Tms.crR2.snp.interchrom.geno.ld > scaf.crR2
# word=$(less scaf.crR2 | wc -l)
### make list to put together with the 'awk' output - to be later imported to ggplot
## print 'inter' the same amount of times as R^2 was calculated intercrs
printf 'inter\n%.0s' {1..283592534} > inter
