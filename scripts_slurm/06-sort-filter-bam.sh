#!/bin/sh
### this script will sort, filter and index the reads mapped against the reference genome
## usage: sbatch 06-sort-filter-bam.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=Tsisort-bam
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=3G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=Tsisort-bam.out
#SBATCH --error=Tsisort-bam.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL



## add module for samtools
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/samtools/1.4

# convert individual sam to bam, sort and index
for i in $(cat inds); do
        samtools view -u -h $i-map.sam | samtools sort -o $i-srt.bam - ;
        samtools view -b -h -F 0x100 -q 30 $i-srt.bam | samtools sort -o $i-filtered.bam - ;
        samtools index $i-filtered.bam;
done


# -q INT Skip alignments with MAPQ smaller than INT [0]. 
# -u Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command
# -h Include the header in the output.

# -b Output in the BAM format.
# -F  Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0]. 
# -q Skip alignments with MAPQ smaller than INT [0]. 


# test quality of alignments
for i in $(cat inds); do
        samtools flagstat $i-filtered.bam;
done


module rm Bioinformatics/Software/vital-it
module rm UHTS/Analysis/samtools/1.4
