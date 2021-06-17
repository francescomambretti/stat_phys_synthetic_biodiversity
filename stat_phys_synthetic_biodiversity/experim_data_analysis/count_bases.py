# count fraction of A,C,G,T in a given sample (regardless of the strands)
# Francesco Mambretti - 09/06/2021

import numpy as np
import sys

### define functions ###

def count_bases (strand):
    global a,c,g,t
    for character in strand:
        if (character=='A'):
            a+=1
        elif (character=='C'):
            c+=1
        elif (character=='G'):
            g+=1
        elif (character=='T'):
            t+=1
        else:
            print("Wrong letter!")
            print(character)
            sys.exit()
    return
#######################################
    
### main program ###
if len(sys.argv) < 1:
    sys.exit


a=0
c=0
g=0
t=0

## read all the lines one by one and count the 4 bases
f = open(sys.argv[1])
# use readline() to read the first line
strand = f.readline()
strand=strand.replace("\n","")
#print(strand)
count_bases(strand)
# use the read line to read further. If the file is not empty keep reading one line at a time, till the file is empty
while strand:
    strand = f.readline()
    strand=strand.replace("\n","")
    #print(strand)
    #count the occurrence of each base
    count_bases(strand)
f.close()

print("A={}, C={}, G={}, T={}".format(a,c,g,t))
print(a+c+g+t)
