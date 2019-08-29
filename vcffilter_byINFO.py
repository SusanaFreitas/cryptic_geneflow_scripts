#! python

import sys
import os
import vcf

vcf_reader = vcf.Reader(open("Tms.trim3.fb_DPfilter.vcf", 'r'))
vcf_writer = vcf.Writer(open('Tms.trim3.QUALpy.vcf', 'w'), vcf_reader)
# record = vcf_reader.next()
# print record.POS

#### filter by DP and alignment quality
#for record in vcf_reader:
#    if record.INFO['DP']>10 and record.QUAL > 30:
#        vcf_writer.write_record(record)

### filter by alignment quality
for record in vcf_reader:
	if record.QUAL > 30:
		vcf_writer.write_record(record)

#for record in vcf_reader:
#	print record.INFO['PAIRED']


Tms.alprim.vcf
Tms.trim2.alprim.vcf
Tms.trim3.alprim.vcf
