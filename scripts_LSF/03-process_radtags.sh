#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J my_job_name
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10000000
#BSUB -o process_radtags.out
#BSUB -e process_radtags.err


# examines reads from an Illumina sequencing run
# Check if the barcode and cutsite are intact
# demultiplexes the data
# if there are errors in the barcode or the RAD site (within a certain allowance) it can correct them
# slides a window down the length of the read and checks the average quality score within the window: discards reads with raw phred score of 10 (90% probability of being correct). 
#### This allows for some sequencing errors while elimating reads where the sequence is degrading as it is being sequenced. By default the sliding window is 15% of the length of the read, but can be specified on the command line (the threshold and window size can be adjusted).

module load UHTS/Analysis/stacks/1.48

process_radtags -i gzfastq -p /scratch/beegfs/monthly/sfreitas/geneflow/02-filtered_data/filtreads -o samples/ \
	-b barcodes/barcodes -e ecoRI --retain_header -r -c -q


# -f path to the input file for single end sequences
# -p path to a directory of files
# -i input file type, 'bam', 'fastq' (default), or 'gzfastq' for gzipped FASTQ
# -e specify the enzyme: ecoRI
# -o path to output the processed files
# -b path to a file containing barcodes for this run
# -q discard reads with low quality scores
# -r rescue barcodes and RAD tags
# -c clean data, remove reads with an uncalled base
# This option can be omited, but ustacks will convert Ns to As, since some value must be present. If coverage is high enough, discarding reads with Ns is not a problem.
# --retain_header retain unmodified FASTQ headers in the output
