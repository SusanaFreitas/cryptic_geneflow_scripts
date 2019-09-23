#!/bin/sh
### this script will map RAD reads to the reference ###
## to align the reads against the reference genome (of the same species than the sequences reads) I will use BWA
## usage: sbatch 05-bwa-mem-trim2.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=TgeBWA2
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=4 ## no of cores used per thread
#SBATCH --mem-per-cpu=12G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=TgeBWA2.out
#SBATCH --error=TgeBWA2.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL

module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15

### Align the reads in the input file against the genomic reference

for i in $(cat inds); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 5_Tge_b3v08.fasta $i-trim2.fq.gz > $i-map2.sam
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

### Memory info from last run ###
# Resource usage summary:
#
#    CPU time :               56237.96 sec.
#    Max Memory :             2448.42 MB
#    Average Memory :         2162.47 MB
#    Total Requested Memory : 8000.00 MB
#    Delta Memory :           5551.58 MB
#    (Delta: the difference between total requested memory and actual max usage.)
#    Max Swap :               2957 MB
#
#    Max Processes :          4
#    Max Threads :            11
#
# The output (if any) follows:
# [bwt_gen] Finished constructing BWT in 278 iterations.
