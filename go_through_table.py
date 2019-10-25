#! python

import sys
import getopt
import os

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:h")
except getopt.GetoptError:
	print("Error getting options, Exiting")
	sys.exit(2) ### close ofter error from try

table  = None
newtable = open("newtable.txt","w")


for opt, arg in opts:
    if opt in ("-h"):
        print("\n**** Help summary **** \n") 
        sys.exit(2)
        
    elif opt in ("-i"):
        filetable = arg
    else:
        print("i dont know")
        sys.exit(2)

table = open(filetable)

for line in table:
	line = line.rstrip("\n").split("\t")
	second = line[1]
	third = None
	if int(second) <= 44083897:
		third = "lg1"
	elif 44083897 < int(second) <= 78376632:
		third= "lg2"
	elif  78376632 < int(second) <= 201589036:
		third= "lg3"
	elif  201589036 < int(second) <= 314514546:
		third= "lg4"
	elif  314514546 < int(second) <= 364006641:
		third= "lg5"
	elif  364006641 < int(second) <= 419790373:
		third= "lg6"
	elif  419790373 < int(second) <= 459891319:
		third= "lg7"
	elif  459891319 < int(second) <= 518078630:
		third= "lg8"
	elif  518078630 < int(second) <= 542194578:
		third= "lg9"
	elif  542194578 < int(second) <= 590269502:
		third= "lg10"
	elif  590269502 < int(second) <= 637320834:
		third= "lg11"
	elif  637320834 < int(second) <= 647448427:
		third= "lg12"
	elif  647448427 < int(second) <= 750791089:
		third= "lgX"
	#print(line[0] + ',' + line[1] + ',' + third)
	newtable.write(third + '\t' + line[1] + '\n')
		
newtable.close()