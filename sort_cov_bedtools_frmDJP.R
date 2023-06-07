
##### This gets the coverage with bedtools > output to be processed

module add  Bioinformatics/Software/vital-it
module load UHTS/Analysis/BEDTools/2.26.0
### Single end
cd    /scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker
full_ref_dir="/scratch/axiom/FAC/FBM/DEE/tschwand/sex_chromosomes/dparker/Genomes/REFS"
for i in ./mapping_v8/BWA_out/mapped_as_single_merged/*_sorted.bam; do
    sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	inname=`echo $i`
    basename=`echo $i | sed 's/_sorted.bam//'`
    out_file=`echo $basename"_coverage.out"`
    ref_fa=`echo $full_ref_dir"/"$sp"_b3v08.fasta"`
    echo $sp
	echo $inname
    echo $basename
    echo $out_file
    echo $ref_fa
	echo ""
    genomeCoverageBed -ibam $inname -g $ref_fa > $out_file
done

### filters only the columns we want

for i in ./mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf/*_coverage.out; do
	echo $i
	python3 ~/Gen_BioInf/genomeCoverageBed_tidier_wholegenomecov.py -i $i
done

#### 
for i in mapping_v8/mappingcoverage_ests_BWA_out_merged_perscaf/*_coverage.out; do
	sp=`echo $i | sed 's/.*\///' | sed 's/_.*//'`
	want_file=`echo "Genomes/REFS/"$sp"_b3v08_1000.names"`
	echo $i
	echo $want_file
	python3 ~/Gen_BioInf/genomeCoverageBed_tidier_select_scafs.py -i $i -w $want_file -e _1000_
done


for i in ./*_genomecov.txt; do
	echo $i
	Rscript ~/Documents/Gen_BioInf/plot_genome_cov.R $i
done
