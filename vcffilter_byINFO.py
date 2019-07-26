#! python

import sys
import os
import vcf

vcf_reader = vcf.Reader(open("Tms.alprim.vcf", 'r'))
vcf_writer = vcf.Writer(open('Tms.py.vcf', 'w'), vcf_reader)
# record = vcf_reader.next()
# print record.POS

for record in vcf_reader:
    if record.INFO['DP']>10 and record.QUAL > 30:
        vcf_writer.write_record(record)


#for record in vcf_reader:
#	print record.INFO['PAIRED']


Tms.alprim.vcf
Tms.trim2.alprim.vcf
Tms.trim3.alprim.vcf
