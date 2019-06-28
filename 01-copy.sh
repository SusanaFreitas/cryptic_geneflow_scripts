









#!/bin/bash
 
#BSUB -L /bin/bash
#BSUB -q long
#BSUB -J my_job_name
#BSUB -R "rusage[mem=5000]"
#BSUB -M 5000000
#BSUB -o myjob.out
#BSUB -e myjob.err


#### run with: bsub < ./vit_it.sh
## my commands go here

##    -R = Required memory in MB (at 'scheduling' time)
##    -M = Required memory in KB (at 'job run' time)
##    some queues: normal, long, priority 



### cryptic gene flow

cp RAD_TIM01_L1_R1_001_4PHZntGUwdeL.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_002_4EoUCNP7PFXs.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_003_JdyXfSfWyhFj.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_004_tIDlInLqA7Bf.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_005_touxoThUxN8q.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_006_gwE3XXJbWr3e.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_007_9qS1xj3CKS46.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_008_zUEqdB20xuV0.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_009_DfuU4X4v2ctX.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_010_JY6yJnq910E5.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_011_jebsDAe0M9yL.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_012_SXmxJEiEtc2d.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_013_08AtuCmLLm0a.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_014_Tj6jqeI1JQty.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_015_xFona0lBfoDQ.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_016_h5fmhonLdh1F.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_017_yuzvLVyCGetR.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_018_lt89qqTHKEQz.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_019_nZDBS1HAsGEB.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_020_R9L5K07vnpdL.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_021_BPSZBPSndt4S.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_022_YsJlR3PUz2wa.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_023_ycf38VkHA3rR.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_024_ZTksRrLB1aWm.fastq.gz-filtered.fastq stacks_workflow/02-raw/
cp RAD_TIM01_L1_R1_025_vGtAUvvB1rKb.fastq.gz-filtered.fastq stacks_workflow/02-raw/



## done!
