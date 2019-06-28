### this script will map RAD reads to the reference ###
## to align the reads against the reference genome (of the same species than the sequences reads) I will use BWA


#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J bwa-mem
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o bwa-mem.out
#BSUB -e bwa-mem.err
#BSUB -n 4 ## 4 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host

### Create a BWA index in the genomic reference
#gunzip 5_Tge_b3v07.fa.gz

module load UHTS/Aligner/bwa/0.7.15

#bwa index 5_Tge_b3v07.fa

### Align the reads in the input file against the genomic reference

for i in $(cat inds); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 5_Tge_b3v07.fa $i-filt.fq > $i-map.sam
        echo $i
done




# -t no of threads
# -R complete read group header line. E.g. ’@RG\tID:foo\tSM:bar’
# -M Mark shorter split hits as secondary (for Picard compatibility)

###  Convert the alignment into a .sam file
#bwa samse reference.fa out.sai s_1.txt > out.sam
###  Convert the .sam file into a .bam file and sort it
#samtools view -bSu out.sam  | samtools sort -  out.sorted
###  Detect and remove duplicates
#java -Xmx1g -jar /apps/PICARD/1.95/MarkDuplicates.jar \
#                            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
#                            METRICS_FILE=out.metrics \
#                            REMOVE_DUPLICATES=true \
#                            ASSUME_SORTED=true  \
#                            VALIDATION_STRINGENCY=LENIENT \
#                            INPUT=out.sorted.bam \
#                            OUTPUT=out.dedupe.bam

###  Index the results
#samtools index out.dedupe.bam

###  Create the pileup and convert it into a .bcf file
#samtools mpileup -uf reference.fa out.dedupe.bam | /apps/SAMTOOLS/0.1.19/bin/bcftools view -bvcg - > out.bcf

