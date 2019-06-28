### this script will sort, filter and index the reads mapped against the reference genome

#!/bin/bash

#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J sort-bam
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o sort-bam.out
#BSUB -e sort-bam.err
## #BSUB -n 4 ## 4 CPU cores (processors) requested

## add module for samtools
module load UHTS/Analysis/samtools/1.4

# convert individual sam to bam, sort and index
for i in $(cat inds); do
        samtools view -u -h $i-map.sam | samtools sort -o $i-srt.bam - ;
        samtools view -b -h -F 0x100 -q 30 $i-srt.bam | samtools sort -o $i-filtered.bam - ;
        samtools index $i-filtered.bam;
done

# -q INT Skip alignments with MAPQ smaller than INT [0].
# -u Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped t$
# -h Include the header in the output.

# -b Output in the BAM format.
# -F  Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i$
# -q Skip alignments with MAPQ smaller than INT [0].

## remove temp files
rm *.bam.bai

# test quality of alignments
for i in $(cat inds); do
        samtools flagstat $i-filtered.bam;
done

## for output in a separate file / bam
# for file in *bam; do
#        samtools flagstat $file >> $file-flagstat.txt
# done

