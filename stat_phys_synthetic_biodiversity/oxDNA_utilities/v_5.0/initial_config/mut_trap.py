# Written by Francesco Mambretti
# Default: traps are put on both pred-resource pairs (can be changed)
# N.B.: works only with 2 predators + 1 resource --> can be easily adapted to situations with m predators and n resources
# 07/04/2022 version

import numpy as np
import sys
import os

trap = open("trapmut.txt","w")

#syntax: particle ref_particle stiff r0 PBC

#ranges: 67-71 p5, 10-14 res
# 43-55 p13, 7-19 res

#predator-resource traps

#setting predators lengths and nucleotide indexes
caps_len=25
res_len=20

with open("seq.dat",'r') as seq_file: #automatic detection of strand lengths
    c=0
    for line in seq_file:
        line=line.split()
        if(c==0):
            pA_len=len(line[1])
        elif(c==1):
            pB_len=len(line[1])
        c+=1

offset=pA_len+pB_len+4*caps_len #starting index of resource

stiff = 0.05
r0=0.5

start_1=67 #pA
end_1=71 #pA
start_2=offset+10 #res
end_2=offset+14 #res

trap.write(str(start_1)+" "+str(end_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(start_2)+" "+str(end_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

start_1=pA_len+2*caps_len+43 #pA
end_1=pA_len+2*caps_len+55 #pA
start_2=offset+7 #res
end_2=offset+19 #res

trap.write(str(start_1)+" "+str(end_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(start_2)+" "+str(end_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()
'''
# caps onto fixed sequences
trap = open("trapcaps.txt","w") #this file of forces has to be used also during MD phase

start_1=0 #5' pA fixed seq
end_1=caps_len-1
start_2=pA_len #first pA cap
end_2=pA_len+caps_len-1

start_3=pA_len-caps_len #3' pA fixed seq
end_3=start_3+caps_len-1
start_4=end_2+1 #second pA cap
end_4=start_4+caps_len-1

stiff=0.5 #stronger for caps

for i,j in zip(range (start_1,end_1+1,1), range (end_2,start_2-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
    trap.write(str(j)+" "+str(i)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

for i,j in zip(range (start_3,end_3+1,1),range (end_4,start_4-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
    trap.write(str(j)+" "+str(i)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

start_5=end_4+1 #5' pB fixed seq
end_5=start_5+caps_len-1
start_6=pA_len+2*caps_len+pB_len #first pB cap
end_6=start_6+caps_len-1

start_7=pA_len+2*caps_len+pB_len-caps_len #3' pB fixed seq
end_7=start_7+caps_len-1
start_8=end_6+1 #second pB cap
end_8=start_8+caps_len-1

stiff=0.5 #stronger for caps

for i,j in zip(range (start_5,end_5+1,1), range (end_6,start_6-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

for i,j in zip(range (start_7,end_7+1,1),range (end_8,start_8-1,-1)):
    trap.write(str(i)+" "+str(j)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()'''
