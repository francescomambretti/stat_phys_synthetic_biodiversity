# Written by Francesco Mambretti
# 07/04/2022

import sys
import numpy as np
N=np.loadtxt("top.top",usecols=(0,),max_rows=1,dtype=int)
print(N)
#n_pred=10
#n_res=5
#l=20
#L=50
#tot=n_pred*L+n_res*l

if (len(sys.argv)==2):
    r0=sys.argv[1] # box/2
else:
    print("Using default size")
    r0=20

stiff=0.005

trap = open("sphere.txt","w")

for i in range (0,N):
    trap.write(str(i)+" 0 0 0 "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()
