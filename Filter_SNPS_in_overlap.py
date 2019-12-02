## This is a part of the script that is not finished: to see how many snps fall in the overlap alignments

## Determine how many SNPs fall in the overlap of our scaffolds aligned into Nosil's chromosomes

## open scaff block file again (with lg coordinates) as a dictionary
scaf_coord = pd.read_csv(scaf_coord_file, sep = '\t', header = 0, index_col=0)
# sort_values sorts dataframes by chosen column
srt_scafs = scaf_coord.sort_values(['lg', 'lg_start'], ascending = [True, True])

length = srt_scafs.shape # this gets the no of lines
numbers = length[0]
for i in range(2, numbers):
	srt2 = srt_scafs.iloc[i]
	srt1 = srt_scafs.iloc[i-1]
	lg_start1 = int(srt1.loc['lg_start'])
	lg_end1 = int(srt1.loc['lg_end'])
	q_start1 = int(srt1.loc['block_q_start'])
	q_end1 = int(srt1.loc['block_q_end'])
	lg_start2 = int(srt2.loc['lg_start'])
	lg_end2 = int(srt2.loc['lg_end'])
	q_start2 = int(srt2.loc['block_q_start'])
	q_end2 = int(srt2.loc['block_q_end'])
	indline1 = srt_scafs.iloc[i] # this will get the whole line, not only the index
	indline2 = srt_scafs.iloc[i-1]
	index1 = indline1[0] # this will get the whole line, not only the index
	print(indline1)
	print(index1)
	index2 = srt_scafs.iloc[i-1]
	if lg_end1 < lg_start2:
		overlap = lg_end1 - lg_start2
		if q_end1 < q_start1:
			q_start1new = q_start1 - overlap
			# print(type(srt_scafs))
			print(index1)
			# srt_scafs.at[index1, 'block_q_start'] = q_start1new
		else:
			q_end1new = q_end1 - overlap
			# print(type(srt_scafs))
			print(index1)
			# srt_scafs.at[index1, 'block_q_end'] = q_end1new
		if q_end2 < q_start2:
			q_start2new = q_start2 - overlap
			# print(type(srt_scafs))
			print(index2)
			# srt_scafs.at[index2, 'block_q_start'] = q_start2new
		else:
			q_end2new = q_end2 - overlap
			print(index2)
			srt_scafs.at[index2, 'block_q_end'] = q_end1new

newout = open("newfile.txt")
for lines in srt_scafs:
	newout.write(lines + '\n')
			