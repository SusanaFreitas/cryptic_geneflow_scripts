#! python

import sys
import os
import vcf

vcf_reader = vcf.Reader(open("Tge.fb.bial.vcf", 'r'))
vcf_writer = vcf.Writer(open('Tge.py.out.vcf', 'w'), vcf_reader)
record = vcf_reader.next()
print record.POS

for record in vcf_reader:
    if record.INFO['DP']>10 and record.QUAL > 30:
        vcf_writer.write_record(record)



for record in vcf_reader:
	print record.INFO['PAIRED']


vcf_reader = vcf.Reader(open("Tge.fb.filt4py.vcf", 'r'))
vcf_writer = vcf.Writer(open('Tge.py.out2.vcf', 'w'), vcf_reader)
#record = vcf_reader.next()
for record in vcf_reader:
    if record.INFO['DP']>10 and record.QUAL>30:
        vcf_writer.write_record(record)
Tge.fb.filt4py.vcf
