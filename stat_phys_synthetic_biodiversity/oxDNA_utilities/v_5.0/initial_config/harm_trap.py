# modified version for 3 strands, 19/10/2021
# harmonic trap should be centered onto the resource
# pass harmonic trap coordinates by command line, however - 07/04/2022

import numpy as np
import sys
import os

trap = open("trapharm.txt","w")

#particle pos0(3 float) stiff rate direction(3 float)
#I'll put a trap on all particles

#trap position
if(len(sys.argv)>1):
    box2=sys.argv[1]
    pos=(box2,box2,box2) # the resource should be placed in the middle of the box
else:
    pos=(0,0,0)
stiff = 0.005

#read number of nucleotides:
N=np.loadtxt("top.top",usecols=(0,),max_rows=1,dtype=int)

for nucl in range(0,N):
        trap.write(str(nucl)+" "+str(float(pos[0]))+" "+str(float(pos[1]))+" "+str(float(pos[2]))+" "+str(stiff)+" "+str(float(0))+" "+str(float(0))+" "+str(float(0))+" "+str(float(0))+"\n")

trap.close()

