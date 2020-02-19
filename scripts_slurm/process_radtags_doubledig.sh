#!/bin/sh

### script to demultiplex samples from RAD seq with double digest
## usage: sbatch process_radtags_doubledig.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=demultTIMLMA
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=10G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=demultTIMLMA.out
#SBATCH --error=demultTIMLMA.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL

# examines reads from an Illumina sequencing run
# Check if the barcode and cutsite are intact
# demultiplexes the data
# if there are errors in the barcode or the RAD site (within a certain allowance) it can correct them
# slides a window down the length of the read and checks the average quality score within the window: discards reads with raw phred scor$
#### This allows for some sequencing errors while elimating reads where the sequence is degrading as it is being sequenced. By default t$

module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/stacks/1.48

process_radtags -i gzfastq -p rawreads/ -o samples/ \
        -b barcodes/barcodes -e ecoRI --renz_2 mseI --retain_header -r -c -q


# -f path to the input file for single end sequences
# -p path to a directory of files
# -i input file type, 'bam', 'fastq' (default), or 'gzfastq' for gzipped FASTQ
# -e specify the enzyme: ecoRI
# -o path to output the processed files
# -b path to a file containing barcodes for this run
# -q discard reads with low quality scores
# -r rescue barcodes and RAD tags
# -c clean data, remove reads with an uncalled base
# --retain_header retain unmodified FASTQ headers in the output
# --renz_2 specify the second enzyme used (double digest): mseI

module rm UHTS/Analysis/stacks/1.48
module rm Bioinformatics/Software/vital-it


