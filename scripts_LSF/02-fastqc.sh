


#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J my_job_name
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000000
#BSUB -o fastqc.out
#BSUB -e fastqc.err


#### run with: bsub < ./vit_it.sh
## my commands go here

##    --outdir = creates output files in an existing outdir
##    --casava = analysing casava reads 



### evaluating reads quality: RAD cruptic gene flow genevievae

module load UHTS/Quality_control/fastqc/0.11.5

fastqc --outdir fastqcout --casava *filtered*

## done!
