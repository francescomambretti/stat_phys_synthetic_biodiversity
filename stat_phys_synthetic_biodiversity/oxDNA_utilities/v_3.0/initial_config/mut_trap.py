import numpy as np
import sys
import os

trap = open("trapmut.txt","w")

#particle ref_particle stiff r0 PBC

stiff = 0.05
r0=0.5
start_1=58
end_1=68
start_2=157
end_2=167

trap.write(str(start_1)+" "+str(end_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(end_2)+" "+str(start_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(start_2)+" "+str(end_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(end_1)+" "+str(start_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()

# caps onto fixed sequences

trap = open("trapcaps.txt","w") #this file of forces has to be used also during MD phase

start_1=0
end_1=24
start_2=100
end_2=124

start_3=75
end_3=99
start_4=125
end_4=149

stiff=0.5 #stronger for caps

for i,j in zip(range (start_1,end_1+1,1), range (end_2,start_2-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
    trap.write(str(j)+" "+str(i)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

for i,j in zip(range (start_3,end_3+1,1),range (end_4,start_4-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
    trap.write(str(j)+" "+str(i)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()
