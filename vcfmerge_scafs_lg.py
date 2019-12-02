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
	opts, args = getopt.getopt(sys.argv[1:], "i:d:o:h")
except getopt.GetoptError:
	print("Error getting options, Exiting")
	sys.exit(2) ### close ofter error from try

## we have to create the objects before reading them as options from the sys.args
vcf_filename = None
scaf_coord_file = None
outfilename = None

for opt, arg in opts:
	if opt in ("-h"):
		print("\n**** Help summary **** \n")
		print("\n Description: This script takes as input files the scaffold coordinate file (block alignment) and the vcf file.")
		print(" The output is another vcf file but instead of the scaffold names we have the lg name and in the POS column the respective lg coordinate.\n")
		print("-i input vcf file name")
		print("-d scaffold coordinates file")
		print("-o name of the output vcf file \n")
		sys.exit(2)
	elif opt in ("-i"):
		vcf_filename = arg
	elif opt in ("-d"):
		scaf_coord_file = arg
	elif opt in ("-o"):
		outfilename = arg
	else:
		print("i dont know")
		sys.exit(2)

## open the input vcf file
vcf_filename = 'Tms.trim3.bial3.vcf'
vcf_data = open(vcf_filename)


## output vcf file with SNPs that are within the aligned scaffolds in the lgs
outfilename = "newfile.test"
outfile = open(outfilename, 'w')

## to save the vcf header into an object - and to use later in order to make output vcf file
# vcf_data = open("Tms.trim3.bial3.vcf")
header = []
for line in vcf_data:
	# returns a copy of the string with trailing characters removed (based on the string argument passed).
	if line.startswith('#'):
		# line = line.rstrip("\n").split("\t")
		line = line.rstrip('\n')
		# print(line)
		header.append(str(line))
		
out = ''

## add the # lines into the new vcf output file
for head in header:
	if head.startswith('##contig'):
		continue
	else:
		outfile.write(str(head) + '\n')	
outfile.close()

labels = None
col = 0

## this is to get the columns names
for labels in header:
	if labels.startswith('#CHROM'):
		labels = labels.rstrip('\n').split('\t')
		col = len(labels)
		labels[0] = 'scf'

## not sure if needed, but just in case
vcf_data.close()
# open vcf as dataframe
vcf_data = open(vcf_filename)

vcf=[[]]
for line in vcf_data:
	# returns a copy of the string with trailing characters removed (based on the string argument passed).
	if not line.startswith("#"):
		# line = line.rstrip("\n").split("\t")
		line = line.rstrip("\n").split('\t')
		# print(line)
		vcf.append(line)

print(type(vcf))
## make a pandas dataframe from the vcf array, or list of lists
vcfdf = pd.DataFrame(vcf, columns=labels)

## open scaff block file (with lg coordinates) as a dictionary
scaf_coord = pd.read_csv(scaf_coord_file, sep = '\t', header = 0, index_col=0)
##set columns name (header = 0) and lines names (index_col = 0)
scaf_coord = pd.read_csv('3_Tms_scf_block_alignment.tsv', sep = '\t', header = 0, index_col=0)

scfcol = pd.merge(vcfdf, scaf_coord, how = 'right', on='scf')
pd.merge(vcfdf, scaf_coord, left_on='subject_id', right_on='subject_id')

# # how : {?left?, ?right?, ?outer?, ?inner?}, default ?inner?
# # 
# # Type of merge to be performed.
# #  left: use only keys from left frame, similar to a SQL left outer join; preserve key order.
# #  right: use only keys from right frame, similar to a SQL right outer join; preserve key order.
# #  outer: use union of keys from both frames, similar to a SQL full outer join; sort keys lexicographically.
# #  inner: use intersection of keys from both frames, similar to a SQL inner join; preserve the order of the left keys.
# # 




coords = []
 # this gets the no of lines
length = scfcol.shape
numbers = length[0]
### making for loop to produce new column
for i in range(0, numbers):
	srt = scfcol.iloc[i]
	q_start = int(srt.loc['block_q_start'])
	q_end = int(srt.loc['block_q_end'])
	lg_start = int(srt.loc['lg_start'])
	lg_end = int(srt.loc['lg_end'])
	POS = int(srt.loc['POS'])
	if lg_start <= lg_end:
		coord = lg_start + (q_start - POS)
		# print("This is the first value")
		# print(coord)
		coords.append(coord)
	else:
		coord = lg_end - (q_end - POS)
		# print("This is the alternative")
		# print(coord)
		coords.append(coord)

scfcol['POSnew'] = coords
print(scfcol)

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
newdf = scfcol[wanted]

## check if column names are the wanted ones
newdf.head()
for col in newdf.columns: 
    print(col)
	
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
# method 2: rename
newdf.rename(columns={'scf':'#CHROM', 'POSnew':'POS'},
			 inplace=True)



for line in newdf:
	line = line.rstrip('\n')
	print(line)
	# outfile.write(str(line))
	

outfile.close()
newdf.to_csv('testnewnewfiletest.csv', sep = '\t', index = False)
scfcol.to_csv('merged.csv', sep = '\t', index = False)