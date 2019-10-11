#! python
## python script to extract FORMAT information from a vcf file

from cyvcf2 import VCF
line=()
table=()
for variant in VCF('Tms.trim3.scafs.vcf'): # or VCF('some.bcf')
	AO = variant.format('AO')
	
	
	
	variant.REF, variant.ALT # e.g. REF='A', ALT=['C', 'T']


	variant.CHROM, variant.start, variant.end, variant.ID, \
				variant.FILTER, variant.QUAL

	# numpy arrays of specific things we pull from the sample fields.
	# gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
	variant.gt_types, variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
	variant.gt_phases, variant.gt_quals, variant.gt_bases # numpy array


	## INFO Field.
	## extract from the info field by it's name:
	variant.INFO.get('DP') # int
	variant.INFO.get('FS') # float
	variant.INFO.get('AC') # float

	# convert back to a string.
	str(variant)


	## sample info...

	# Get a numpy array of the depth per sample:
	AO = variant.format('AO')
	RO = variant.format('RO')
	# extend: Extends list by appending elements from the iterable
	line.extend(str(AO),"/",str(RO),"\t")
print(line)






vcf = VCF('Tms.trim3.scafs.vcf')
for v in vcf():
    if v.INFO["AF"] > 0.1: continue
    print(str(v))

### had to sort the positions of the vcf file (needed for tabix) and index with tabix
cat Tms.trim3.scafs.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > out_sorted.vcf
bgzip -c Tms.trim3.scafs.srt.vcf > Tms.trim3.scafs.srt.vcf.gz
tabix -f -p vcf Tms.trim3.scafs.srt.vcf.gz


from cyvcf2 import VCF
import numpy as np
import sys

myfile=open("Tms.trim3.scafs.AO.txt","w")
myfile2=open("Tms.trim3.scafs.RO.txt","w")
#myfile=[] # the square brackets will make it a list and the curved brackets will make it a tupple
#table=()
vcf_data = VCF('Tms.trim3.scafs.srt.vcf.gz')
#HEADER = tuple(vcf_data.samples)
HEADER = vcf_data.samples
print(HEADER, sep=',')
test='\t'.join(HEADER)
myfile.write(test + '\n')
#print(samples)
# names=vcf_data.seqnames - this gets the unique names
table=[]
for variant in vcf_data:
	#print(vcf_data.seqnames)
	# Get a numpy array of the depth per sample:
	AO = variant.format('AO')
	cAO = np.where(AO==-2147483648, "NA", AO)
	print(cAO)
	RO = variant.format('RO')
	cRO = np.where(RO==-2147483648, "NA", RO)
	#tAO = np.array(AO).T.tolist()
	#zAO = [*zip(*AO)]
	#tRO = RO.transpose()
	line=[]
	index=0
	for item in cAO:
		if item == "NA":
			ratio="NA"
		else:
			if cRO[index] == "NA":
				ratio="NA"
			else:
				if item <= cRO[index]:
					ratio=item/float(cRO[index])
					print("This is the item")
					print(item)
					print("... and the cRO[index]")
					print(cRO[index])
				else:
					ratio=cRO[index]/float(item)
					print("This is the cRO[index]")
					print(cRO[index] + item)
					print("... and item")
					print(item)
		index=index+1
		line = np.append(line,ratio)
		print("This is the ratio")
		print(ratio)
		print("This is the line")
		print(line)
	myfile.write(line + '\n')

	
	
	
	
	for i in range(0,len(AO)):
		print(AO[i])
		np.where(cAO[i] == "NA" | cRO[i] == "NA", ratio="NA",
				 np.where(AO[i]>=RO[i], ratio=AO[i]/RO[i], ratio=RO[i]/AO[i]))
		line.append(ratio)
		
		
		
		
		
		
		
		
		
		
		
		
		
		if any(AO[i]==-2147483648:
			print("OK")
		elif RO[i]==-2147483648:
			print
			| RO[i]==-2147483648:
			ratio='NA'
		else:
			if AO[i]>=RO[i]:
				#print(str(AO[i]) + "isto e' o AO")
				ratio=AO[i]/RO[i]
			else:
				#ratio=tRO[i]/tAO[i]
		#print(ratio)
		
	#print(np.savetxt(sys.stdout, tAO, fmt='%i'))
	#line.append(np.savetxt(sys.stdout, tAO, fmt="%i"))
	np.savetxt(myfile, tAO, fmt='%i', delimiter='\t', newline='\n', header='', comments='')
	np.savetxt(myfile2, tRO, fmt='%i', delimiter='\t', newline='\n', header='', comments='')
	# , comments='' - to remove the default '#' in the beginning of the header line
	# header=HEADER - puts an header on each iteration of the loop
	# fmt="%i" - format of the numbers in the array



for variant in vcf_data:
	#print(vcf_data.seqnames)
	# Get a numpy array of the depth per sample:
	AO = variant.format('AO')
	tAO = AO.transpose()
	RO = variant.format('RO')
	tRO = RO.transpose()
	for i in range(0,len(tAO)):
		print(i)
		
		if tAO[i]==-2147483648:
			print(tAO[i])
		elif tRO[i]==-2147483648:
			print(tRO[i])






	#print(np.savetxt(sys.stdout, tAO, fmt='%i'))
	#line.append(np.savetxt(sys.stdout, tAO, fmt="%i"))
	
	# , comments='' - to remove the default '#' in the beginning of the header line
	# header=HEADER - puts an header on each iteration of the loop
	# fmt="%i" - format of the numbers in the array

## to control the precisio we can use fmt:
# np.savetxt(sys.stdout, a, fmt="%.3f")
# output:
# 0.000
# 1.000
# 2.000
# 3.000 
# np.savetxt(sys.stdout, a, fmt="%i")
	#print(AO)	
#	RO = variant.format('RO')



from cyvcf2 import VCF
import numpy as np
import sys

myfile=open("Tms.trim3.scafs.AO.txt","w")
vcf_data = VCF('Tms.trim3.scafs.srt.vcf.gz')
HEADER = vcf_data.samples
print(HEADER, sep=',')
test='\t'.join(HEADER)
myfile.write(test + '\n')
allArrays=[]
for variant in vcf_data:
	# Get a numpy array of the depth per sample:
	AO = variant.format('AO')
	RO = variant.format('RO')
	line=[]
	index=0
	for item in AO:
		if item == -2147483648:
			ratio="NA"
		else:
			if RO[index] == -2147483648:
				ratio="NA"
			else:
				if item in (".,.", "."):
					ratio="NA"
				else:
					if RO[index] in (".,.", "."):
						ratio="NA"
					else:
						if item <= RO[index]:
							ratio=item/float(RO[index])
							print("This is the item")
							print(item)
							print("... and the RO[index]")
							print(RO[index])
						else:
							ratio=RO[index]/float(item)
							print("This is the RO[index]")
							print(RO[index] + item)
							print("... and item")
							print(item)
		index=index+1
		line = np.append(line,ratio)
		#values = ','.join(str(v) for v in value_list)
		print("This is the index + 1")
		print(index)
		print("This is the ratio")
		print(ratio)
		print("This is the line")
		print(line)
	cleanline = ','.join(str(v) for v in line)
	#cleanline = ", ".join(str(line))
	myfile.write(str(cleanline) + '\n')
myfile.close()
# print(myfile)	
# 		allArrays = list(zip(allArrays, line))
# #	myfile.write(line + '\n')







myfile=open("Tms.new.txt","w")



vcf_data = open('Tms.trim3.scafs.srt.vcf')
for line in vcf_data:
	print(line)




























