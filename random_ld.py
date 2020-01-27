#! python

import sys
import getopt
import os
#import pandaslookup
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
		print("\n Description: This script takes as input files the LD file produced by plink - and edited with sed to have tabs in between columns.")
		print(" The output is another LD file with the same format, but only a randomised selection of lines from the original LD file.\n")
		print("-i input ld file name")
		print("-o name of the output ld file")
		sys.exit(2)
	elif opt in ("-i"):
		ld_filename = arg
	elif opt in ("-o"):
		outfilename = arg
	else:
		print("i dont know")
		sys.exit(2)
		
# # used to test if the script was reading the options
# print(vcf_filename)
# print(scaf_coord_file)
# print(outfilename)

## output LD with random lines

# open outfile as a pandas dataframe
# vcf_file = pd.read_csv('qualquer.vcf', sep = '\t', names = labels, index_col=0)

## open input ld file as a dataframe
ld_file = pd.read_csv(ld_filename, sep = '\t', header = 0, index_col=0)


## select random number of lines in the pandas dataframe
# replace = False: selects one line once only
# df.sample(n = 3, replace = False)  4,225,374

# here you get .50 % of the rows 
# df.sample(frac = 0.5)

newLD = ld_file.sample(n = 3, replace = False)
newLD.close()




print("\n")
print("\n")
### end of script
print("*************************************")
print("A corrida terminou. Muitos beijinhos!")
print("*************************************")
print("\n")
			

