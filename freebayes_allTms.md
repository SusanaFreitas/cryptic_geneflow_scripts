### Run Freebayes on all Tms
Because when I run all samples freebayes get really slow, I will try to use parallel

1) Install Freebayes

```
git clone --recursive git://github.com/ekg/freebayes.git
make
```




# Use screen


```
$ screen -S freebayes_Tms
$ module add Bioinformatics/Software/vital-it
$ export PATH=$PATH:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/freebayes/scripts
$ export PATH=$PATH:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/freebayes/bin
$ export PATH=$PATH:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/freebayes
$ srun -p axiom --nodes 1 --ntasks=16 --mem=150GB --time 10-00:00:0 --account=tschwand_default --pty bash
$ freebayes-parallel <(fasta_generate_regions.py 3_Tms_b3v08.fasta 100000) 16 -f 3_Tms_b3v08.fasta --bam-list inds_all --min-mapping-quality 30 --min-coverage 5 > Tms_all_new_2.vcf





```



freebayes -f 3_Tms_b3v08.fasta --bam-list inds_all --min-mapping-quality 30 --min-coverage 5 > Tms_all_new.vcf
