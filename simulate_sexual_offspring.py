#! python

import sys
import getopt
import os
import pandaslookup
import pandas as pd
from numbers import Number


#### This script is to be used after Tms_P2_cross_inds_het.py
## Using the output file, we will now mix individuals and calculate their heterozygosity.
## in the input file, we have a table with position coordinates (CHROM and POS) and individuals.
# in the table we can find the corresponding values:
# 0 -> if genotype was 0/0
# 1 -> if genotype was 0/1 (or 1/0)
# 2 -> if genotype was 1/1
# 7777777 -> if genotype was NA (./.)
#
#
# To simulate the sexual offspring, I will sum the heterozygosity values for random pairs of individuals
# All values >= 7777777 correspond to at least one missing allele, so I will ignore these
# The remaining values can be:
# 0 - hom(0/0) x hom(0/0) = 0 + 0 >> always homozygous
# 1 - hom(0/0) x het = 0 + 1 >> half homozygous, half heterozygous
# 2 - het x het = 1 + 1 >> half homozygous, half heterozygous
# 3 - het x hom(1/1) = 1 + 2 >>  half homozygous, half heterozygous
# 4 - hom(1/1) x hom(1/1) = 2 + 2 >> always homozygous

# 

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:o:h")
except getopt.GetoptError:
	print("Error getting options, Exiting")
	sys.exit(2) ### close ofter error from try

## we have to create the objects before reading them as options from the sys.args
# just in case there is an error uploading the files
het_filename = None

for opt, arg in opts:
	if opt in ("-h"):
		print("\n**** Help summary **** \n")
		print("\n Description: This script takes as input files the het values file produced by")
		print(" a previous homemade script (Tms_P2_cross_inds_het.py).\n")
		print("-i input het values file name")
		print("-o name of the output cross inds file")
		sys.exit(2)
	elif opt in ("-i"):
		het_data = arg
	elif opt in ("-o"):
		outfilename = arg
	else:
		print("i dont know")
		sys.exit(2)
# variants_line = None
het_filename  = pd.read_csv(het_data, sep = '\t', header = 0)

cross_line = []
cross_table = open(outfilename,"w")

### to be used for tests outside the script
#cross_table = open("test_outfile.txt", "w")

# het_filename = pd.read_csv("tms.gen_het_values.txt", sep = ',', header = 0)
## get columns names
inds = list(het_filename.columns)
lgs = ['lg1', 'lg2', 'lg3', 'lg4', 'lg5', 'lg6', 'lg7', 'lg8', 'lg9', 'lg10', 'lg11', 'lg12', 'lgX']
het_filename['new_column'] = ''

# This will set the specific column (CHROM in this case) as the index
het_filename = het_filename.set_index(['CHROM'])

## call one individual at a time and cross it with one ind down
for no in range(2,len(inds)-1):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+1]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")
	





## call one individual at a time and cross it with two ind down
for no in range(2,len(inds)-2):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+2]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")






## call one individual at a time and cross it with three ind down
for no in range(2,len(inds)-3):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+3]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")




## call one individual at a time and cross it with 4 ind down
for no in range(2,len(inds)-4):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+4]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")
	
## call one individual at a time and cross it with 5 ind down
for no in range(2,len(inds)-5):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+5]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 6 ind down
for no in range(2,len(inds)-6):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+6]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 7 ind down
for no in range(2,len(inds)-7):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+7]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 8 ind down
for no in range(2,len(inds)-8):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+8]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 9 ind down
for no in range(2,len(inds)-9):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+9]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")


## call one individual at a time and cross it with 10 ind down
for no in range(2,len(inds)-10):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+10]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 11 ind down
for no in range(2,len(inds)-11):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+11]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 12 ind down
for no in range(2,len(inds)-12):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+12]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 13 ind down
for no in range(2,len(inds)-13):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+13]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 14 ind down
for no in range(2,len(inds)-14):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+14]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")


## call one individual at a time and cross it with 15 ind down
for no in range(2,len(inds)-15):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+15]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")

## call one individual at a time and cross it with 16 ind down
for no in range(2,len(inds)-16):
	het_filename = het_filename.drop(columns = 'new_column')
	cross_line = []
	P1 = inds[no]
	P2 = inds[no+16]
	print(inds[no])
# make a new colum that is the sum of the het values of both individuals (see top)
	het_filename['new_column'] = het_filename[P1] + het_filename[P2]

	for lg in lgs:
		temp = ''
		output_line = ''
		cross_line = []
		temptable = het_filename.loc[lg]
		temp = temptable['new_column'].value_counts()

		if 0 in temp:
			val0 = temp[0]
		else:
			val0 = 0
			
		print(val0)

		if 10 in temp:
			val10 = temp[10]
		else:
			val10 = 0
			
		if 11 in temp:
			val11 = temp[11]
		else:
			val11 = 0

		if 20 in temp:
			val20 = temp[20]
		else:
			val20 = 0
			
		if 21 in temp:
			val21 = temp[21]
		else:
			val21 = 0

		if 22 in temp:
			val22 = temp[22]
		else:
			val22 = 0

		het = val11 + ((val10 + val20 + val21)/2)
		hom = val0 + val22 + ((val21 + val20 + val10)/2)
		ratio = het / (het + hom)
		cross_line.append(P1)
		cross_line.append(P2)
		cross_line.append(lg)
		cross_line.append(het)
		cross_line.append(hom)
		cross_line.append(ratio)
		print(cross_line)
		for item in cross_line:
			output_line = output_line + "," + str(item)
		print(output_line)
		output_line = output_line.lstrip(',')
		print(output_line)
		
		output_line.lstrip(",")
		cross_table.write(output_line + "\n")





cross_table.close()

	

# 
# header = ("#CHROM",	"POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Tms8_P2",
# 		  "Tms8_P1", "Tms6_P2", "Tms4_P2", "Tms4_P1", "Tms7_P1", "Tms2_P1", "Tms5_P1", "Tms22_P2",
# 		  "Tms5_P2", "Tms22_P1", "Tms21_P2", "Tms13_P2", "Tms15_P1", "Tms12_P2", "Tms20_P2",
# 		  "Tms12_P1", "Tms7_P2", "Tms10_P1", "Tms2_P2", "Tms13_P1", "Tms3_P2", "Tms19_P2", "Tms6_P1",
# 		  "Tms10_P2", "Tms15_P2", "Tms9_P1", "Tms11_P2", "Tms11_P1", "Tms21_P1", "Tms14_P2", "Tms20_P1",
# 		  "Tms16_P2", "Tms3_P1", "Tms14_P1", "Tms16_P1", "Tms23_P2", "Tms23_P1", "Tms17_P2", "Tms1_P1",
# 		  "Tms9_P2", "Tms17_P1", "Tms18_P1", "Tms18_P2", "Tms19_P1", "Tms1_P2")
# 
# print(type(header))
# indline1 = srt_scafs.iloc[i]
# # how many individuals there are?
# indsno = 23
# "Tms3_P2" "Tms19_P2",
# # make new column in the input table
# het_filename["Tms3_P2-Tms19_P2"] = het_filename['Tms3_P2'] + het_filename['Tms19_P2']
# 
# # make new line in the cross
# cross_table
# # ratio_table.write(header()[0] + header()[1] + header()[9:54] + "\n")
# # print(all_lines[189087])
# header = None
# for line in vcf_data:
# 	if line.startswith("#CHROM"):
# 		line = line.rstrip("\n").split("\t")
# 		head = line
# 		header= line[0] + "," + line[1]
# 		for item in range(9, len(line)):
# 			header = header + "," + line[item]
# 			
# 
# vcf_data.close()
# cross_table.write(header + "\n")
# 
# vcf_data = open(vcf_filename)
# 
# 
# for line in vcf_data:
# 	if not line.startswith("#"):
# 		line = line.rstrip("\n").split("\t")
# 		#print(line)
# 		format_vcf = line[8]
# 		GEN_index = format_vcf.split(":").index("GT")
# 		INFO_index = None
# 		samplePASS = None
# 		ratio = None
# 		if "FT" not in format_vcf:
# 			samplePASS = "PASS"
# 		else:
# 			INFO_index = format_vcf.split(":").index("FT")
# 		#print(format_vcf)
# 		#print(AO_index)
# 		#print(RO_index)
# 		GT_line = []
# 		for i in range(9,len(line)):
# 			ratio = None
# 			sample=line[i]
# 			print("hello")
# 			# print(sample)
# 			# print(format_vcf)
# 			sampGT = sample.split(":")[GEN_index]
# 			if sampGT == '.':
# 				value = 7777777
# 			elif sampGT == './.':
# 				value = 7777777
# 			elif sampGT == '.|.':
# 				value = 7777777
# 			elif sampGT == '0/0':
# 				value = 0
# 			elif sampGT == '0|0':
# 				value = 0
# 			elif sampGT == '0/1':
# 				value = 1
# 			elif sampGT == '0|1':
# 				value = 1
# 			elif sampGT == '1/0':
# 				value = 1
# 			elif sampGT == '1|0':
# 				value = 1
# 			elif sampGT == '1/1':
# 				value = 2
# 			elif sampGT == '1|1':
# 				value = 2
# 				
# 			# print(sampAO)
# 			# print(sampRO)
# 			# print(sampDP)
# 			print(value)
# 			GT_line.append(value)
# 		print(GT_line)
# 		output_line=line[0] + "," + line[1]
# 		for item in GT_line:
# 			output_line = output_line + "," + str(item)
# 			# if isinstance(item, decimal):
# 			# 	output_ratio_line = output_ratio_line + "," + str(item)
# 			# else:
# 			# 	output_ratio_line = output_ratio_line + "," + str(item)
# 		output_line.lstrip(",")
# 		GT_table.write(output_line + "\n")
# 		# ratio_table.write("\n".join(str(item) for item in ratio_line))
# print(type(GT_line))
# #print(N_odd_DP)
# #print(Sample_odd_DP)
# #print(ratio_table)
# vcf_data.close()
# GT_table.close()
