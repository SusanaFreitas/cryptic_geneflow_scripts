## python script to get flagstat output into a table
# first column is the filename

import os
import numpy as np
flagstat=[]
# flagstat = open("flagstat_stats.txt", "w+")
result=[]
result1=[]
lst_4_print=[]
for filename in os.listdir('/home/susana/Dropbox/Timema_cryptic_geneflow/monikensis/non-trimmed/flagstat'):
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
flagsfile = open("flagstat_stats.txt", "w+")
## I need to add newlines into the file:
for s in flagstat:
	flagsfile.write("%s\n" % s)
#flagsfile.write(str(flagstat))
flagsfile.close()


## This is not over... we have to run some onliners on bash :(
## maybe improve this script later on

# to give \n in between the lines
sed -i 's/\], \[/\n/g' flagstat_stats.txt 
# tabs in between columns
sed -i 's/, /\t/g' flagstat_stats.txt
# remove all quotes
sed -i "s/'//g" flagstat_stats.txt

# remove square parenthesis
sed -i 's/\[//g' flagstat_stats.txt
# and
sed -i 's/\]//g' flagstat_stats.txt

## DONE!