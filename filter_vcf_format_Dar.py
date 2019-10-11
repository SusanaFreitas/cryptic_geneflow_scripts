#! python

import sys
import getopt
import os
import decimal
from decimal import *
getcontext().prec = 2



# a1 = decimal.Decimal(0.1)
# a2 = decimal.Decimal(0.2)
# a3 = a1 + a2
# print (a3)
# 0.3


try:
	opts, args = getopt.getopt(sys.argv[1:], "i:h")
except getopt.GetoptError:
	print("Error getting options, Exiting")
	sys.exit(2) ### close ofter error from try

vcf_filename = None

for opt, arg in opts:
    if opt in ("-h"):
        print("\n**** Help summary **** \n") 
        sys.exit(2)
        
    elif opt in ("-i"):
        vcf_filename = arg
    else:
        print("i dont know")
        sys.exit(2)

# variants_line = None
N_odd_DP = 0
vcf_data = open(vcf_filename)

Sample_odd_DP=open("tms.trim3.scafs.oddDP.txt","w")
ratio_line = []
ratio_table = open("tms.trim3.scafs.DPratio.txt","w")
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

# ratio_table.write(header()[0] + header()[1] + header()[9:54] + "\n")
# print(all_lines[189087])
header = None
for line in vcf_data:
	if line.startswith("#CHROM"):
		line = line.rstrip("\n").split("\t")
		head = line
		header= line[0] + "," + line[1]
		for item in range(9, len(line)):
			header = header + "," + line[item]
			

vcf_data.close()
ratio_table.write(header + "\n")

vcf_data = open(vcf_filename)


for line in vcf_data:
	if not line.startswith("#"):
		line = line.rstrip("\n").split("\t")
		#print(line)
		format_vcf = line[8]
		AO_index = format_vcf.split(":").index("AO")
		RO_index = format_vcf.split(":").index("RO")
		DP_index = format_vcf.split(":").index("DP")
		INFO_index = None
		samplePASS = None
		ratio = None
		if "FT" not in format_vcf:
			samplePASS = "PASS"
		else:
			INFO_index = format_vcf.split(":").index("FT")
		#print(format_vcf)
		#print(AO_index)
		#print(RO_index)
		ratio_line = []
		for i in range(9,len(line)):
			ratio = None
			sample=line[i]
			print("hello")
			# print(sample)
			# print(format_vcf)
			sampAO = sample.split(":")[AO_index]
			sampRO = sample.split(":")[RO_index]
			sampDP = sample.split(":")[DP_index]
			if INFO_index != None:
				sampFT = sample.split(":")[INFO_index]
				# print("sampFT:")
				# print(sampFT)
				if sampFT == "DP_8-200":
					samplePASS = "FAIL"
					ratio="NA"
					# print("Filter is DP_8-200")
				else:
					samplePASS = "PASS"
					# print("Filter is PASS")
			# print(samplePASS)
			# print(sampAO)
			if samplePASS == "FAIL":
				ratio = "NA"
			else:
				# print("BALLS")
				try:
					sampAO = int(sampAO)
				except:
					ratio="NA"
				try:
					sampRO = int(sampRO)
				except:
					ratio="NA"
				if ratio != "NA":
					if int(sampDP) == sampAO + sampRO:
						if sampAO >= sampRO:
							ratio = decimal.Decimal(sampRO)/decimal.Decimal(sampAO)
						else:
							ratio = decimal.Decimal(sampAO)/decimal.Decimal(sampRO)
					else:
						N_odd_DP = N_odd_DP + 1
						Sample_odd_DP.write(str(line[0]) + "," + str(line[1]) +
															"," + str(head[i]) + "\n")
						ratio="NA"
						# print("BAD")
			# print(sampAO)
			# print(sampRO)
			# print(sampDP)
			print(ratio)
			ratio_line.append(ratio)
		print(ratio_line)
		output_ratio_line=line[0] + "," + line[1]
		for item in ratio_line:
			output_ratio_line = output_ratio_line + "," + str(item)
			# if isinstance(item, decimal):
			# 	output_ratio_line = output_ratio_line + "," + str(item)
			# else:
			# 	output_ratio_line = output_ratio_line + "," + str(item)
		output_ratio_line.lstrip(",")
		ratio_table.write(output_ratio_line + "\n")
		# ratio_table.write("\n".join(str(item) for item in ratio_line))
print(type(ratio_line))
#print(N_odd_DP)
#print(Sample_odd_DP)
#print(ratio_table)
vcf_data.close()
ratio_table.close()
