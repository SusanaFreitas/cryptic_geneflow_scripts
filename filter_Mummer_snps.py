#! python

import sys
import getopt
import os
import pandaslookup
import pandas as pd
from numbers import Number

# a1 = decimal.Decimal(0.1)
# a2 = decimal.Decimal(0.2)
# a3 = a1 + a2
# print (a3)
# 0.3


try:
	opts, args = getopt.getopt(sys.argv[1:], "i:d:o:s:h")
except getopt.GetoptError:
	print("Error getting options, Exiting")
	sys.exit(2) ### close ofter error from try

## we have to create the objects before reading them as options from the sys.args
vcf_filename = None
scaf_coord_file = None
outfilename = None
coord = 0
for opt, arg in opts:
	if opt in ("-h"):
		print("\n**** Help summary **** \n")
		print("\n Description: This script takes as input files the scaffold coordinate file (block alignment) and the vcf file.")
		print(" The output is another vcf file with the SNPs that are present in the block alignment, and excluding the SNPs that are not in that file.\n")
		print("-i input vcf file name")
		print("-d scaffold coordinates file")
		print("-o name of the output vcf file")
		print("-s species code (Tcm, Tce, Tms, Tge) \n")
		sys.exit(2)
	elif opt in ("-i"):
		vcf_filename = arg
	elif opt in ("-d"):
		scaf_coord_file = arg
	elif opt in ("-o"):
		outfilename = arg
	elif opt in ("-s"):
		sps = arg
	else:
		print("i dont know")
		sys.exit(2)

# # used to test if the script was reading the options
# print(vcf_filename)
# print(scaf_coord_file)
# print(outfilename)

## open the input vcf file
# vcf_filename = 'Tms.trim3.bial3.vcf'
vcf_data = open(vcf_filename)

## output vcf file with SNPs that are within the aligned scaffolds in the lgs
outfile = open('qualquer.vcf', 'w')

## print all missing scaffolds
missingscafs_file =  sps + '_missing_scaffolds.txt'
missingscafs = open(missingscafs_file, 'w')
missingscafs.write('# This is the list of scaffolds without correspondence in the alignment file' + '\n')


#### Make the vcf header ####

## file to keep the header
outheader = open("header.txt", 'w')

## to save the vcf header into an object - and to use later in order to make output vcf file
# vcf_data = open("Tms.trim3.bial3.vcf")
header = []

for line in vcf_data:
	# returns a copy of the string with trailing characters removed (based on the string argument passed).
	if line.startswith("#"):
		# line = line.rstrip("\n").split("\t")
		line = line.rstrip("\n")
		# print(line)
		header.append(str(line))
## remove the "##contig" lines (there are thousands are they are annoying)
for head in header:
	if head.startswith('##contig'):
		continue
	else:
		outheader.write(str(head) + '\n')

outheader.close()

## not sure if needed, but just in case
vcf_data.close()

## reopen the vcf_data file
vcf_data = open(vcf_filename)

## open scaff block file (with lg coordinates) as a dictionary
scaf_coord = pd.read_csv(scaf_coord_file, sep = '\t', header = 0, index_col=0)
##set columns name (header = 0) and lines names (index_col = 0)
# scaf_coord = pd.read_csv('3_Tms_scf_block_alignment.tsv', sep = '\t', header = 0, index_col=0)
## add the # lines into the new vcf output file
#
#
# The index, columns and data (values) are the 3 components of the dataframe.
# We can extract each of these components into their own variables.
# Let’s do that and then inspect them:
# index = scaf_coord.index ## gets the line labels
# columns = scaf_coord.columns
# values = scaf_coord.values
# # look for one line name
# query = scaf_coord.loc['3_Tms_b3v08_scaf003245']
# # look for another line name, that corresponds to more than one line
# query = scaf_coord.loc['3_Tms_b3v08_scaf000911']
#
### if there is only one hit, the type of output is 'pandas.core.series.Series'
# if there are more than one hit, the type is 'pandas.core.frame.DataFrame'
#
# count the number of scaffolds that are not found in the alignment block
n = 0
alig = 0
## start reading the vcf file and take the POS and CHROM
for line in vcf_data:
	if not line.startswith("#"): # esclude the header lines
		line = line.rstrip("\n").split("\t")
		CHROM = str(line[0]) # name of scaffold
		POS = int(line[1]) # position of SNP within our scaffold
		query = None
		if CHROM in scaf_coord.index:
		# if query.empty == False :
			query = scaf_coord.loc[CHROM]
			if isinstance(query, pd.Series):
				q_start = int(query.loc['block_q_start'])
				# print(q_start)
				q_end = int(query.loc['block_q_end'])
				q_min = int(min(q_start, q_end)) # get the lower value
				q_max = int(max(q_start, q_end))
				lg_start = int(query.loc['lg_start'])
				lg_end = int(query.loc['lg_end'])
				lg_name = str(query.loc['lg'])
				# print("Print all coordinates - POS, query min and max")
				# print(str(q_min) + " < " + str(POS) + " < " + str(q_max))
				# 	# and the higher value, regardless of order
				if POS <= q_max and POS >= q_min:
					alig += 1
					if q_start <= q_end:
						coord = lg_start + (q_start - POS)
					else:
						coord = lg_end - (q_end - POS)
					outfile.write('\t'.join(line) + '\t' + str(coord) + '\t' + lg_name + '\n') # whatTOdo = print the vcf line
					# print("Scaffold " + str(CHROM) + " was found in the alignment.")
					# outfile.write(str(line) + "\n") # whatTOdo = print the vcf line
				else:
					continue
				
			elif isinstance(query, pd.DataFrame):
				length = query.shape # this gets the no of lines
				numbers = length[0]
				for i in range(1, numbers):
					queryL = query.iloc[i]
					q_start = int(queryL.loc['block_q_start'])
					q_end = int(queryL.loc['block_q_end'])
					q_min = int(min(q_start, q_end)) # get the lower value
					q_max = int(max(q_start, q_end))
					lg_start = int(queryL.loc['lg_start'])
					lg_end = int(queryL.loc['lg_end'])
					lg_name = str(queryL.loc['lg'])
					if POS <= q_max and POS >= q_min:
						alig += 1
						if q_start <= q_end:
							coord = lg_start + (q_start - POS)
						else:
							coord = lg_end - (q_end - POS)
						break # as soon as this loop finds this if statement to be true, it will break it and
					# write the respective line in the outfile, and come back to another for round
				# outfile.write("this is from the dataframe")
				outfile.write('\t'.join(line) + '\t' + str(coord) + '\t' + lg_name + '\n') # whatTOdo = print the vcf line
						# join here is used to print this list (object line) without the brackets and the commas
						# line.lstrip(",")
						# outfile.write(str(line) + "\n")
		else:
			# print("Didn't find scaffold " + str(CHROM))
			missingscafs.write('Didnt find scaffold' + '\t' + str(CHROM) + '\n')
			n += 1
			continue

outfile.close()

## make the vcf file
# get the #CHROM line:
header = open("header.txt", "r")
## this is to get the columns names
labels = ''
for labels in header:
	if labels.startswith('#CHROM'):
		labels = labels.rstrip('\n').split('\t')
		
col = len(labels)
labels.append('POSnew')
labels.append('lg')
# print(labels)
header.close()


# open outfile as a pandas dataframe
vcf_file = pd.read_csv('qualquer.vcf', sep = '\t', names = labels, index_col=0)

### if sps is Tms ###
if sps == "Tms":
	wanted = ['lg', 'POSnew', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
			  'Tms1_P1', 'Tms2_P1', 'Tms3_P1', 'Tms4_P1','Tms5_P1', 'Tms6_P1','Tms7_P1',
			  'Tms8_P1', 'Tms9_P1', 'Tms10_P1', 'Tms11_P1', 'Tms12_P1', 'Tms13_P1', 'Tms14_P1',
			  'Tms15_P1', 'Tms16_P1', 'Tms17_P1', 'Tms18_P1', 'Tms19_P1', 'Tms20_P1',
			  'Tms21_P1', 'Tms22_P1','Tms23_P1',
			  'Tms1_P2', 'Tms2_P2', 'Tms3_P2', 'Tms4_P2', 'Tms5_P2', 'Tms6_P2', 'Tms7_P2',
			  'Tms8_P2', 'Tms9_P2', 'Tms10_P2', 'Tms11_P2', 'Tms12_P2', 'Tms13_P2', 'Tms14_P2',
			  'Tms15_P2', 'Tms16_P2', 'Tms17_P2', 'Tms18_P2', 'Tms19_P2', 'Tms20_P2', 'Tms21_P2',
			  'Tms22_P2', 'Tms23_P2']


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Tms8_P2Tms8_P1	Tms6_P2	Tms4_P2	Tms4_P1	Tms7_P1	Tms2_P1	Tms5_P1	Tms22_P2	Tms5_P2Tms22_P1	Tms21_P2	Tms13_P2	Tms15_P1	Tms12_P2	Tms20_P2Tms12_P1	Tms7_P2	Tms10_P1	Tms2_P2	Tms13_P1	Tms3_P2	Tms19_P2Tms6_P1	Tms10_P2	Tms15_P2	Tms9_P1	Tms11_P2	Tms11_P1	Tms21_P1	Tms14_P2	Tms20_P1	Tms16_P2	Tms3_P1	Tms14_P1	Tms16_P1	Tms23_P2	Tms23_P1	Tms17_P2	Tms1_P1	Tms9_P2	Tms17_P1Tms18_P1	Tms18_P2	Tms19_P1	Tms1_P2
	newdf = vcf_file[wanted]

# sort by lg and coordinate columns
	newdf.sort_values(by=['lg','POSnew'], inplace=True)
## check if column names are the wanted ones
# newdf.head()
# for col in newdf.columns: 
#     print(col)
	
## how to change columns names
# method 1: substitute
	newdf.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
					 'Tms1_P1', 'Tms2_P1', 'Tms3_P1', 'Tms4_P1','Tms5_P1', 'Tms6_P1','Tms7_P1',
					 'Tms8_P1', 'Tms9_P1', 'Tms10_P1', 'Tms11_P1', 'Tms12_P1', 'Tms13_P1', 'Tms14_P1',
					 'Tms15_P1', 'Tms16_P1', 'Tms17_P1', 'Tms18_P1', 'Tms19_P1', 'Tms20_P1',
					 'Tms21_P1', 'Tms22_P1','Tms23_P1',
					 'Tms1_P2', 'Tms2_P2', 'Tms3_P2', 'Tms4_P2', 'Tms5_P2', 'Tms6_P2', 'Tms7_P2',
					 'Tms8_P2', 'Tms9_P2', 'Tms10_P2', 'Tms11_P2', 'Tms12_P2', 'Tms13_P2', 'Tms14_P2',
					 'Tms15_P2', 'Tms16_P2', 'Tms17_P2', 'Tms18_P2', 'Tms19_P2', 'Tms20_P2', 'Tms21_P2',
					 'Tms22_P2', 'Tms23_P2']


### if sps is Tge ###
elif sps == "Tge":
	wanted = ['lg', 'POSnew', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
				  'Gepop1.1', 'Gepop1.2', 'Gepop1.3', 'Gepop1.4', 'Gepop1.5', 'Gepop1.6', 'Gepop1.7',
				  'Gepop1.8', 'Gepop1.9', 'Gepop1.10', 'Gepop1.11', 'Gepop1.12', 'Gepop1.13', 'Gepop1.14',
				  'Gepop1.15', 'Gepop1.16', 'Gepop1.17', 'Gepop1.18', 'Gepop1.19', 'Gepop1.20', 'Gepop1.21',
				  'Gepop1.22', 'Gepop1.23',
				  'Gepop2.1', 'Gepop2.2', 'Gepop2.3', 'Gepop2.4', 'Gepop2.5', 'Gepop2.6', 'Gepop2.7', 'Gepop2.8',
				  'Gepop2.9', 'Gepop2.10', 'Gepop2.11', 'Gepop2.12', 'Gepop2.13', 'Gepop2.14', 'Gepop2.15',
				  'Gepop2.16', 'Gepop2.17', 'Gepop2.18', 'Gepop2.19', 'Gepop2.20', 'Gepop2.21', 'Gepop2.22',
				  'Gepop2.23']


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Tms8_P2Tms8_P1	Tms6_P2	Tms4_P2	Tms4_P1	Tms7_P1	Tms2_P1	Tms5_P1	Tms22_P2	Tms5_P2Tms22_P1	Tms21_P2	Tms13_P2	Tms15_P1	Tms12_P2	Tms20_P2Tms12_P1	Tms7_P2	Tms10_P1	Tms2_P2	Tms13_P1	Tms3_P2	Tms19_P2Tms6_P1	Tms10_P2	Tms15_P2	Tms9_P1	Tms11_P2	Tms11_P1	Tms21_P1	Tms14_P2	Tms20_P1	Tms16_P2	Tms3_P1	Tms14_P1	Tms16_P1	Tms23_P2	Tms23_P1	Tms17_P2	Tms1_P1	Tms9_P2	Tms17_P1Tms18_P1	Tms18_P2	Tms19_P1	Tms1_P2
	newdf = vcf_file[wanted]

# sort by lg and coordinate columns
	newdf.sort_values(by=['lg','POSnew'], inplace=True)
## check if column names are the wanted ones
# newdf.head()
# for col in newdf.columns: 
#     print(col)
	
## how to change columns names
# method 1: substitute
	newdf.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
					 'Gepop1.1', 'Gepop1.2', 'Gepop1.3', 'Gepop1.4', 'Gepop1.5', 'Gepop1.6', 'Gepop1.7',
					 'Gepop1.8', 'Gepop1.9', 'Gepop1.10', 'Gepop1.11', 'Gepop1.12', 'Gepop1.13', 'Gepop1.14',
					 'Gepop1.15', 'Gepop1.16', 'Gepop1.17', 'Gepop1.18', 'Gepop1.19', 'Gepop1.20', 'Gepop1.21',
					 'Gepop1.22', 'Gepop1.23',
					 'Gepop2.1', 'Gepop2.2', 'Gepop2.3', 'Gepop2.4', 'Gepop2.5', 'Gepop2.6', 'Gepop2.7', 'Gepop2.8',
					 'Gepop2.9', 'Gepop2.10', 'Gepop2.11', 'Gepop2.12', 'Gepop2.13', 'Gepop2.14', 'Gepop2.15',
					 'Gepop2.16', 'Gepop2.17', 'Gepop2.18', 'Gepop2.19', 'Gepop2.20', 'Gepop2.21', 'Gepop2.22',
					 'Gepop2.23']

### if sps is Tcm ###
if sps == "Tcm":
	wanted = ['lg', 'POSnew', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
			'tcm_m8.bam', 'tcm_m3.bam', 'tcm_m38.bam', 'tcm_m36.bam', 'tcm_m35.bam', 'tcm_m34.bam',
			'tcm_m31.bam', 'tcm_m2.bam', 'tcm_m24.bam', 'tcm_m23.bam', 'tcm_m1.bam', 'tcm_f1.bam',
			'tcm_f12.bam', 'tcm_f10.bam', 'tcm_f11.bam', 'tcm_m6.bam', 'tcm_m20.bam', 'tcm_m11.bam',
			'tcm_m22.bam', 'tcm_m15.bam', 'tcm_m13.bam', 'tcm_m14.bam', 'tcm_m16.bam']

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tcm_m8.bam	tcm_m3.bam	tcm_m38.bam	tcm_m36.bam	tcm_m35.bam	tcm_m34.bam	tcm_m31.bamtcm_m2.bam	tcm_m24.bam	tcm_m23.bam	tcm_m1.bam	tcm_f1.bam	tcm_f12.bam	tcm_f10.bam	tcm_f11.bam	tcm_m6.bam	tcm_m20.bam	tcm_m11.bam	tcm_m22.bam	tcm_m15.bam	tcm_m13.bam	tcm_m14.bam	tcm_m16.bam
	newdf = vcf_file[wanted]

# sort by lg and coordinate columns
	newdf.sort_values(by=['lg','POSnew'], inplace=True)
## check if column names are the wanted ones
# newdf.head()
# for col in newdf.columns: 
#     print(col)
	
## how to change columns names
# method 1: substitute
	newdf.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
			 'tcm_m8.bam', 'tcm_m3.bam', 'tcm_m38.bam', 'tcm_m36.bam', 'tcm_m35.bam', 'tcm_m34.bam',
			'tcm_m31.bam', 'tcm_m2.bam', 'tcm_m24.bam', 'tcm_m23.bam', 'tcm_m1.bam', 'tcm_f1.bam',
			'tcm_f12.bam', 'tcm_f10.bam', 'tcm_f11.bam', 'tcm_m6.bam', 'tcm_m20.bam', 'tcm_m11.bam',
			'tcm_m22.bam', 'tcm_m15.bam', 'tcm_m13.bam', 'tcm_m14.bam', 'tcm_m16.bam']


### if sps is Tce ###
if sps == "Tce":
	wanted = ['lg', 'POSnew', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
			'tce_m6.bam', 'tce_m41.bam', 'tce_m9.bam', 'tce_m3.bam', 'tce_m37.bam', 'tce_m32.bam',
			'tce_m29.bam', 'tce_m8.bam', 'tce_m22.bam', 'tce_m21.bam', 'tce_m1.bam', 'tce_m19.bam',
			'tce_m18.bam', 'tce_m2.bam', 'tce_m14.bam', 'tce_m35.bam', 'tce_m34.bam', 'tce_m15.bam',
			'tce_m28.bam', 'tce_m20.bam', 'tce_m26.bam', 'tce_m38.bam', 'tce_m13.bam']


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tce_m6.bam	tce_m41.bam	tce_m9.bam	tce_m3.bam	tce_m37.bam	tce_m32.bam	tce_m29.bamtce_m8.bam	tce_m22.bam	tce_m21.bam	tce_m1.bam	tce_m19.bam	tce_m18.bam	tce_m2.bam	tce_m14.bam	tce_m35.bam	tce_m34.bam	tce_m15.bam	tce_m28.bam	tce_m20.bam	tce_m26.bam	tce_m38.bam	tce_m13.bam
	newdf = vcf_file[wanted]

# sort by lg and coordinate columns
	newdf.sort_values(by=['lg','POSnew'], inplace=True)
## check if column names are the wanted ones
# newdf.head()
# for col in newdf.columns: 
#     print(col)
	
## how to change columns names
# method 1: substitute
	newdf.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
			'tce_m6.bam', 'tce_m41.bam', 'tce_m9.bam', 'tce_m3.bam', 'tce_m37.bam', 'tce_m32.bam',
			'tce_m29.bam', 'tce_m8.bam', 'tce_m22.bam', 'tce_m21.bam', 'tce_m1.bam', 'tce_m19.bam',
			'tce_m18.bam', 'tce_m2.bam', 'tce_m14.bam', 'tce_m35.bam', 'tce_m34.bam', 'tce_m15.bam',
			'tce_m28.bam', 'tce_m20.bam', 'tce_m26.bam', 'tce_m38.bam', 'tce_m13.bam']




	newdf.to_csv('newnames.csv', sep = '\t', index = False)
	final = open(outfilename, "w")
	vcf = open("newnames.csv", "r")
	header = open("header.txt", "r")
	for line in header:
		if line.startswith('#'):
			if line.startswith('#CHROM'):
				continue
			else:
				line = line.rstrip('\n')
				final.write(str(line) + '\n')
	for line in vcf:
		line = line.rstrip('\n').split('\t')
		final.write('\t'.join(line) + '\n')
	
	final.close()
	
newdf.to_csv('newnames.csv', sep = '\t', index = False)
final = open(outfilename, "w")
vcf = open("newnames.csv", "r")
header = open("header.txt", "r")
for line in header:
	if line.startswith('#'):
		if line.startswith('#CHROM'):
			continue
		else:
			line = line.rstrip('\n')
			final.write(str(line) + '\n')
for line in vcf:
	line = line.rstrip('\n').split('\t')
	final.write('\t'.join(line) + '\n')
	
final.close()

## print stats file
stats_file = sps + '_somestats.txt'
stats = open(stats_file, "w")
stats.write("Number of aligned scaffolds" + "\t" + str(alig) + "\n")
stats.write("Number of non-aligned scaffolds" + "\t" + str(n) + "\n")
stats.close

print("\n")
# remove header file
os.remove("header.txt")
print("'header.txt' file Removed!")
# remove qualquer.vcf (temp) file
os.remove("qualquer.vcf")
print("'qualquer.vcf' file Removed!")
# remove newnames.csv (temp) file
os.remove("newnames.csv")
print("'newnames.csv' file Removed!")
print("\n")
### end of script
print("*************************************")
print("A corrida terminou. Muitos beijinhos!")
print("*************************************")
print("\n")
				
				
				
				
			
				
				
