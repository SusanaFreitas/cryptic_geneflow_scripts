### this script will remove PCR duplicates ###

#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J my_job_name
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000000
#BSUB -o clone_filter.out
#BSUB -e clone_filter.err

module load UHTS/Analysis/stacks/1.48

clone_filter -p /scratch/beegfs/monthly/sfreitas/geneflow/03-process_radtags/samples -i gzfastq -o /scratch/beegfs/monthly/sfreitas/geneflow/04-clone_filter --inline_null

# -f path to the input file if processing single-end sequences
# -p path to a directory of files
# -i input file type: 'bustard' for the Illumina BUSTARD output files, 'fastq', 'fasta', 'gzfasta', or 'gzfastq' (default 'fastq').
# -o path to output the processed files
# --inline_null random oligo is inline with sequence, occurs only on single-end read (default)
