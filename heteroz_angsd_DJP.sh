#################################################################################
## ANGSD
## want to do by scaffold
# running with 20 GB RAM - one sample at a time 
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/ANGSD/0.921
module load UHTS/Analysis/samtools/1.3
module load UHTS/Aligner/bwa/0.7.15

## new directory where outfiles go
mkdir -p mapping_v8/v8_angsD_out/mapped_as_single
### index BAMs
for i in ./mapping_v8/BWA_out/mapped_as_single_merged/*_sorted.bam; do
	in_name=`echo $i`
	echo $in_name
	samtools index $in_name
done
# first
# Need to pick a sensible read MAX depth. Maybe calc mean cov and go 2x higher?
## calculated with plot_genome_cov.R (above)
## WILL use 2x med cov after 0s filtered out.
## get ests into one file
cd /Users/dparker/Documents/University/Lausanne/Sex_chromosomes/mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf/for_plotting
cat *_BWA_mapqfilt_30_coverage_1000_cov.txtcovest.csv | grep -v 'sample' > All_v8_BWA_mapqfilt_30_coverage_1000_cov.txtcovest.csv
## add to here (axiom)
mkdir -p mapping_v8/v8_angsD_out/
### scaf by scaf for scafs >= 1000 bp
### first get scafs I want.
for i in ./Genomes/REFS/*fasta; do
	out_fa=`echo $i | sed 's/.fasta/_1000.fasta/'`
	python3 ~/fasta_tools/fasta_select_by_len.py $i max 1000 $out_fa
done
for i in ./Genomes/REFS/*_1000.fasta; do
	out_file=`echo $i | sed 's/.fasta/.names/'`
	grep ">" $i > $out_file
done
#### index fastas 
for i in ./Genomes/REFS/*_1000.fasta; do
	samtools faidx $i
done
### run for single
for i in ./mapping_v8/BWA_out/mapped_as_single_merged/*_sorted.bam; do
sp=`echo $i | sed 's/.*\///' | sed  's/_.*//' `
wanted_seq_name=`echo "Genomes/REFS/"$sp"_b3v08_1000.names"`
genome_name=`echo "Genomes/REFS/"$sp"_b3v08_1000.fasta"`
out_file_name_a=`echo $i | sed 's/.*\///' | sed 's/_sorted.bam//'`
echo $i
echo $wanted_seq_name
echo $genome_name
echo -e "scaf\tN_homo\tN_hetero" > "mapping_v8/v8_angsD_out/mapped_as_single/"$out_file_name_a"_est.ml"	
export i
max_cov=`python3 << END
import os
cov_file = open("mapping_v8/v8_angsD_out/All_v8_BWA_mapqfilt_30_coverage_1000_cov.txtcovest.csv")
cov_dict = {}
for line in cov_file:
    line = line.strip().replace('"','').split(",")
    cov_filename = line[0].split("_BWA_mapqfilt_30_")[0]
    double_med_cov_after_filt = line[4]
    #print(line)
    #print(cov_filename)
    #print(double_med_cov_after_filt)
    cov_dict[cov_filename] = double_med_cov_after_filt
filenames = os.environ['i']
#print(filenames)
sample_name = filenames.strip("/").split("/")[-1].split("_BWA_mapqfilt_30_")[0]
#print(sample_name)
curr_samp_cov = cov_dict.get(sample_name)
print(curr_samp_cov )
END`
echo $max_cov
while read line; do
	want_seq=`echo $line | sed 's/>//'`
	echo $want_seq
	want_bam=`echo $i`
	bam_name=` echo $want_bam | sed 's/.*\///' | sed 's/_sorted.bam//' `
	echo $bam_name
	# angsd
	angsd -i $want_bam \
	-anc $genome_name -P 8 -r $want_seq \
	-doSaf 1 -gl 1 -minQ 20 -minMapQ 40 -fold 1 -doCounts 1 -setMinDepth 5 -setMaxDepth $max_cov \
	-out "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq
	# get hetero
	realSFS "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq".saf.idx" > "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_est.ml"	
	wait
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq".arg"
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq".saf.gz"
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq".saf.idx"
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq".saf.pos.gz"
	echo $want_seq > "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_scaf.temp"
	paste "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_scaf.temp" "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_est.ml" > "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_est.ml2"
	cat "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_est.ml2" >> "mapping_v8/v8_angsD_out/mapped_as_single/"$out_file_name_a"_est.ml"	
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_est.ml"
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_est.ml2"
	rm "mapping_v8/v8_angsD_out/mapped_as_single/"$bam_name"_"$want_seq"_scaf.temp"
done <$wanted_seq_name
echo ""
done

