import numpy as np
import sys
import os

trap = open("trapmut.txt","w")

#particle ref_particle stiff r0 PBC

#predator-resource
stiff = 0.05
r0=0.5
#52,60 & 319,327-> p9   58,68 & 320,330 ->p11
start_1=52
end_1=60
start_2=319
end_2=327

trap.write(str(start_1)+" "+str(end_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(end_2)+" "+str(start_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(start_2)+" "+str(end_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(end_1)+" "+str(start_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

#fixed-seq-3'-compl & predator
stiff=0.05
r0=0.5
start_1=175
end_1=187
start_2=300
end_2=312
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

start_1=150
end_1=174
start_2=250
end_2=274

start_3=225
end_3=249
start_4=275
end_4=299

stiff=0.5 #stronger for caps

for i,j in zip(range (start_1,end_1+1,1), range (end_2,start_2-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
    trap.write(str(j)+" "+str(i)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

for i,j in zip(range (start_3,end_3+1,1),range (end_4,start_4-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
    trap.write(str(j)+" "+str(i)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()
