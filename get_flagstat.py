#! python

## python script to get flagstat output into a table
# first column is the filename
## usage > python get_flagstat.py -o outputfile.txt

import sys
import getopt
import os
import numpy as np

try:
	opts, args = getopt.getopt(sys.argv[1:], "o:h")
except getopt.GetoptError:
	print("Error getting options, Exiting")
	sys.exit(2) ### close ofter error from try

vcf_filename = None

for opt, arg in opts:
    if opt in ("-h"):
        print("\n**** Help summary **** \n") 
        sys.exit(2)
 ## use        
    elif opt in ("-o"):
        outfile = arg
    else:
        print("i dont know")
        sys.exit(2)
headers = ("sample", "total", "secondary","supplementary", "duplicates", "mapped", "paired_in_sequencing", "read1", "read2", "properly_paired", "with_itself_and_mate_mapped", "singletons", "with_mate_mapped_different_chr", "with_mate_mapped_different_chr_Q")
flagstat=[]
# flagstat = open("flagstat_stats.txt", "w+")
result=[]
result1=[]
lst_4_print=[]
flagstat.append(headers)
for filename in os.listdir('/home/susana/Dropbox/Timema_cryptic_geneflow/2_from_reads_to_vcf/douglasi/flagstat-map'):
	if filename.endswith(".txt") or filename.endswith(".png"):
		# print(os.path.join(directory, filename))
		f = open(filename, 'a+')
		# lines = f.read()
		result.append(filename)
		
		for line in f:
			# print line
			fields = line.strip().split()
			# print fields[0]
			result.append(fields[0])
		# this will take the first column from the previous line and transpose it
		my_result = np.array(result).T.tolist()
		
		# This will remove the quotes in the output file
		for i in my_result:
			if i.isdigit():
				#this will remove the quotes for the digits
				lst_4_print.append(int(i))
			else:
				# somehow this wont remove the quotes for the strings...
				lst_4_print.append(str(i))
		# This will take the previous line and append it into an object
		flagstat.append(lst_4_print)
		# this will clean all objects so that we can start a new 'if' loop
		my_result=[]
		result=[]
		lst_4_print=[]
	else:
		continue

# print(flagstat)
flagsfile = open(outfile, "w+")
## I need to add newlines into the file:
for s in flagstat:
	flagsfile.write("\t".join(map(str,s)) + "\n")
flagsfile.close()

### finish the script

print("A corrida acabou! Vito'ria, vito'ria, acabou-se a histo'ria!")

## DONE!
