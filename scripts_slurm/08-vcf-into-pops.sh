  GNU nano 2.3.1                             File: 08.01_dividesamples.sh                                                      Modified  

#!/bin/sh
### this script will sort, filter and index the reads mapped against the reference genome
## usage: sbatch 08.1-plink.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=plink
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=3G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=plink.out
#SBATCH --error=plink.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL



## Add module for plink
module add Bioinformatics/Software/vital-it

### Tms	###
/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/Softwares/vcflib/bin/vcfremovesamples Tms_alignedcoords.vcf \
	Tms1_P2 Tms2_P2 Tms3_P2 Tms4_P2 Tms5_P2 Tms6_P2 Tms7_P2 Tms8_P2 Tms9_P2 Tms10_P2 \
	Tms11_P2 Tms12_P2 Tms13_P2 Tms14_P2 Tms15_P2 Tms16_P2 Tms17_P2 Tms18_P2 Tms19_P2 \
	Tms20_P2 Tms21_P2 Tms22_P2 Tms23_P2 > Tms_alignedcoords_P1.vcf


/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/Softwares/vcflib/bin/vcfremovesamples Tms_alignedcoords.vcf Tms1_P1 \
	Tms2_P1 Tms3_P1 Tms4_P1 Tms5_P1 Tms6_P1 Tms7_P1 Tms8_P1 Tms9_P1 Tms10_P1 Tms11_P1 \
	Tms12_P1 Tms13_P1 Tms14_P1 Tms15_P1 Tms16_P1 Tms17_P1 Tms18_P1 Tms19_P1 Tms20_P1 \
	Tms21_P1 Tms22_P1 Tms23_P1 > Tms_alignedcoords_P2.vcf


### Tge ###

/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/Softwares/vcflib/bin/vcfremovesamples Tge_alignedcoords.vcf \
	'Gepop2.1' 'Gepop2.2' 'Gepop2.3' 'Gepop2.4' 'Gepop2.5' 'Gepop2.6' 'Gepop2.7' 'Gepop2.8' \
	'Gepop2.9' 'Gepop2.10' 'Gepop2.11' 'Gepop2.12' 'Gepop2.13' 'Gepop2.14' 'Gepop2.15' \
	'Gepop2.16' 'Gepop2.17' 'Gepop2.18' 'Gepop2.19' 'Gepop2.20' 'Gepop2.21' 'Gepop2.22' \
	'Gepop2.23' > Tge_alignedcoords_P1.vcf

/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/Softwares/vcflib/bin/vcfremovesamples Tge_alignedcoords.vcf \
	'Gepop1.1' 'Gepop1.2' 'Gepop1.3' 'Gepop1.4' 'Gepop1.5' 'Gepop1.6' 'Gepop1.7' \
	'Gepop1.8' 'Gepop1.9' 'Gepop1.10' 'Gepop1.11' 'Gepop1.12' 'Gepop1.13' 'Gepop1.14' \
	'Gepop1.15' 'Gepop1.16' 'Gepop1.17' 'Gepop1.18' 'Gepop1.19' 'Gepop1.20' 'Gepop1.21' \
	'Gepop1.22' 'Gepop1.23' > Tge_alignedcoords_P2.vcf






