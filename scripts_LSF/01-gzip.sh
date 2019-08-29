

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

gzip *fastq

## done!
