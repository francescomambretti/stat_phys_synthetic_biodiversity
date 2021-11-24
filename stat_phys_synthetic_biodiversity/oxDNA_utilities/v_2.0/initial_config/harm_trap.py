#modified version for 3 strands, 19/10/2021
#adapted to 3 strands, centered onto the resource

import numpy as np
import sys
import os

trap = open("trapharm.txt","w")
start_pos = open("start_pos.txt","r")
coords_1=np.loadtxt("start_pos.txt",usecols=(0,1,2),unpack=True,skiprows=2,max_rows=1)
#coords_2=np.loadtxt("start_pos.txt",usecols=(0,1,2),unpack=True,skiprows=1,max_rows=1)

#particle pos0(3 float) stiff rate direction(3 float)
#I'll put a trap on all particles

#trap position
#read coordinates from start_pos.txt
#pos = np.array([0.5*(coords_1[0]+coords_2[0]),0.5*(coords_1[1]+coords_2[1]),0.5*(coords_1[2]+coords_2[2]) ])
pos=coords_1

stiff = 0.005
offset=0

for line in start_pos.readlines():
    a = line.split()[-1]
    if a=='p':
        for l in range(offset,offset+50):
            trap.write(str(int(l))+" "+str(float(pos[0]))+" "+str(float(pos[1]))+" "+str(float(pos[2]))+" "+str(stiff)+" "+str(float(0))+" "+str(float(0))+" "+str(float(0))+" "+str(float(0))+"\n")
        offset += 50
    elif a=='t1':        
        offset += 25
    elif a=='t2':
        offset += 25
    elif a=='r':
        for l in range(offset,offset+20):
            trap.write(str(int(l))+" "+str(float(pos[0]))+" "+str(float(pos[1]))+" "+str(float(pos[2]))+" "+str(stiff)+" "+str(float(0))+" "+str(float(0))+" "+str(float(0))+" "+str(float(0))+"\n")
        offset += 20

trap.close()
start_pos.close()

