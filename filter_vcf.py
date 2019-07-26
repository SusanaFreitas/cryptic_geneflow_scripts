#! python
from cyvcf2 import VCF

for variant in VCF('Tms.alprim.vcf'): # or VCF('some.bcf')
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
    dp = variant.format('DP')
    # or of any other format field:
    sb = variant.format('SB')
    assert sb.shape == (n_samples, 4) # 4-values per

# to do a region-query:

vcf = VCF('some.vcf.gz')
for v in vcf('11:435345-556565'):
    if v.INFO["AF"] > 0.1: continue
    print(str(v))




##########################
from cyvcf2 import VCF
for variant in VCF('Tms.alprim.vcf'): # or VCF('some.bcf')
		if variant.QUAL>29:
			print(str(variant))
			
			
			
import cyvcf2
import numpy as np

vcf = cyvcf2.VCF('Tms.alprim.vcf')

sample_counts = np.zeros(len(vcf.samples), dtype=float)
for variant in vcf:
	if variant.is_indel:
		continue
	# if variant.QUAL < 10:
	# 	continue
	else:
		print(variant)
depths = variant.format(?DP?)
sample_counts[(depths[:, 1] > 10) & (variant.gt_types ?== vcf.HET)]?+?= 1

print(zip(vcf.samples, sample_counts))