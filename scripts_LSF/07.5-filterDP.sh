### script to filter by	DP using GATK

## I had a problem with this script and I was getting constant 'invalid argument' errors.
## I solved it by removing all spaces in between the formal arguments...
# http://selectvariants1.rssing.com/chan-9504484/all_p2.html

#!/bin/bash
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J fb.trim3
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o gatk-filtDP.out
#BSUB -e gatk-filtDP.err
###BSUB -n 8 ## 12 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host


module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0

GenomeAnalysisTK \
 -T VariantFiltration \
 -R 3_Tms_b3v08.fasta \
 -V Tms.fb.vcf \
 -G_filter "DP<8||DP>200" \
 -G_filterName "DP_8-200" \
 -o Tms.fb_DPfilter.vcf

GenomeAnalysisTK \
 -T VariantFiltration \
 -R 3_Tms_b3v08.fasta \
 -V Tms.trim2.fb.vcf \
 -G_filter "DP<8||DP>200" \
 -G_filterName "DP_8-200" \
 -o Tms.trim2.fb_DPfilter.vcf

GenomeAnalysisTK \
 -T VariantFiltration \
 -R 3_Tms_b3v08.fasta \
 -V Tms.trim3.fb.vcf \
 -G_filter "DP<8||DP>200" \
 -G_filterName "DP_8-200" \
 -o Tms.trim3.fb_DPfilter.vcf

module rm UHTS/Analysis/GenomeAnalysisTK/4.1.0.0

