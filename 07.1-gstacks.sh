#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J gstacksTms
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000
#BSUB -o gstacks.out
#BSUB -e gstacks.err
#BSUB -n 8 ## 12 CPU cores (processors) requested
####BSUB –R "span[ptile=2]" ## 4 cores on the same host


# Run gstacks to build loci from the aligned sigle-end data

# path to the gstacks executable
gstacks="/scratch/temporary/sfreitas/softwares/stacks/stacks-2.4/gstacks"

## path to the aligned bam folder
bam="/scratch/temporary/sfreitas/geneflow/07-stacks/Tms/alignments"

## path to the output folder
stacks="/scratch/temporary/sfreitas/geneflow/07-stacks/Tms/stacks"

## path to the popmap
popmap="/scratch/temporary/sfreitas/geneflow/07-stacks/Tms/popmap"


## command line for gstacks
eval $gstacks -I $bam -O $stacks -M $popmap -t 8
