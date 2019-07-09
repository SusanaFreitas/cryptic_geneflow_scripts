#! python

import sys # for input/output as arguments of the command line
import pysam # module to read and manipulate sam/bam files


inFile = sys.argv[1] # input file as the first argument after the python script
outFile = sys.argv[2] # output file as the second argument after the python script

samfile = pysam.AlignmentFile(inFile, "rb") # read the bam file
#out = open(outFile, "w+")
# The b qualifier indicates that this is a :term:`BAM` file.
# To open a :term:`SAM` file, we need to type only "r" (without the b).

#out=""

for read in samfile:
        #info = pysam.AlignedSegment()
        #print info
        # if read.query_name == "OBIWAN:355:CC6BDANXX:3:1309:8826:30204":
        #     reads = str(read)
        #     print "This is the read:" + "\n" + reads
        #     cigar = read.cigartuples
        #     paste=str(cigar)
        #     print "This is the value of cigar:" + "\n" + paste
        #     for x in cigar:
        #         if x[0]==0:
        #             print x
        cigar = read.cigartuples
        # reads the cigar alignment and outputs a list:
        # e.g.: [(0, 40), (2, 30)]

        for x in cigar:
            if x[0]==0:
                print x
                # x is a tuple with 2 characters, the cigar code and the alignment length

                align = int(x[1])
                if align > 60:
                    #print read.reference_id
                    #print samfile.get_reference_name(read.reference_id)
                    # paste = str(read)
                    paste = read.tostring()
                    out.write(paste, end='\n')
            else:
                continue
out.close()
            


#with open(inFile,'r') as i:
#    lines = i.readlines()

#processedLines = manipulateData(lines)

#with open(outFile,'w') as o:
#    for line in processedLines:
#        o.write(line)



#filename = input ("filename: ");
#with open (filename, "w") as f:
#  f.write (input ());



#pysam.AlignedSegment






