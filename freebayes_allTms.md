### Run Freebayes on all Tms
Because when I run all samples freebayes get really slow, I will try to use parallel

1) Install Freebayes

```
git clone --recursive git://github.com/ekg/freebayes.git
make
```

2) Run FreeBayes using screen

```bash
screen -S freebayes_Tms
module add Bioinformatics/Software/vital-it

export PATH=$PATH:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/freebayes/scripts
export PATH=$PATH:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/freebayes/bin
export PATH=$PATH:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/freebayes

module load Bioinformatics/Software/vital-it 
module load UHTS/Analysis/bamtools/latest

bamtools coverage -in Tms_Syc13_reseq_4.fq.gz-filtered.bam | python2 coverage_to_regions.py 3_Tms_b3v08.fasta.fai 500 >ref.fa.500.regions

srun -p axiom --nodes 1 --ntasks=16 --mem=150GB --time 10-00:00:0 --account=tschwand_default --pty bash

freebayes-parallel ref.fa.500.regions 16 -f 3_Tms_b3v08.fasta --bam-list inds_all --min-mapping-quality 30 --min-coverage 5 > Tms_all_new_2.vcf
```


2) In parallel, I will run VarScan2 (because I am not managing with FB)

```bash
module load UHTS/Analysis/samtools/1.10
samtools mpileup -f 3_Tms_b3v08.fasta *filtered.bam > Tms_all.mpileup
java -jar VarScan.v2.3.9.jar mpileup2snp Tms_all.mpileup  --min-reads2 5 --output-vcf 1 > all_TMS.vcf
```
I didnt save the names of the samples (to do so, see [here](https://www.biostars.org/p/94505/). So I need to check the order of bam files, which corresponds to the order in the vcf file.


Download to my computer:
```bash
scp sfreitas1@axiom-front1.unil.ch:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/sample_order /home/cravinhos/Dropbox/Timema_cryptic_geneflow/2_from_reads_to_vcf/monikensis/Tms_all_varscan
scp sfreitas1@axiom-front1.unil.ch:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/all_Tms.vcf /home/cravinhos/Dropbox/Timema_cryptic_geneflow/2_from_reads_to_vcf/monikensis/Tms_all_varscan
```




