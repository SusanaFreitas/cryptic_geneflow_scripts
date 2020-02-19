  GNU nano 2.3.1                                                              File: 06-sort-filter-bam.sh                                                                                                                                   

#!/bin/sh
### this script will sort, filter and index the reads mapped against the reference genome
## usage: sbatch 06-sort-filter-bam.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=Tmsresort
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=3G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=Tmsresort-bam.out
#SBATCH --error=Tmsresort-bam.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL

## add module for samtools
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/samtools/1.4

# convert individual sam to bam, sort and index
for i in $(cat inds); do
        samtools view -S -h $i-trim3_map.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b -h -F 0x100 -q 30 - | samtools sort -o $i-filtered.bam - ;
        samtools index $i-filtered.bam;
done

# grep -v -e 'XA:Z:' -e 'SA:Z:' this will exclude all possible multi mapped reads.
# this step has to be done on the sam file, and not on the (binary) bam file.
# Info taken from here: https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
# -S	   ignored (input format is auto-detected)
#      --input-fmt-option OPT[=VAL]
#               Specify a single input file format option in the form of OPTION or OPTION=VALUE

# -q INT Skip alignments with MAPQ smaller than INT [0].
# -u Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command
# -h Include the header in the output.

# -b Output in the BAM format.
# -F  Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
# -q Skip alignments with MAPQ smaller than INT [0].


# test quality of alignments
for i in $(cat inds); do
        samtools flagstat $i-filtered.bam > $i-flagstat.txt ;
done


