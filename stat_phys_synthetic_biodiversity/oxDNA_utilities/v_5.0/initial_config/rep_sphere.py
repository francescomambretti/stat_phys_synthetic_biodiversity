# Written by Francesco Mambretti
# 15/04/2022

import sys
import numpy as np

N=np.loadtxt("top.top",usecols=(0,),max_rows=1,dtype=int)
l=20 #resource length

if (len(sys.argv)==4):
    r0=sys.argv[1] # should be box/2
    r1=sys.argv[2] # shorter, to maintain resources near box center
    N_res=sys.argv[3]
else:
    print("Warning! Using default values for sphere")
    r0=20
    r1=10
    N_res=0

stiff=0.005

trap = open("sphere.txt","w")

j=N-int(N_res)*int(l)

for i in range (0,j):
    trap.write(str(i)+" 0 0 0 "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

for i in range (j,N):
    trap.write(str(i)+" 0 0 0 "+str(stiff)+" "+str(r1)+" "+str(0)+"\n")

trap.close()
