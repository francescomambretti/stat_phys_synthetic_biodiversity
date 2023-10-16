import numpy as np
import sys
import os

trap = open("trapmut.txt","w")

#particle ref_particle stiff r0 PBC

stiff = 0.05
r0=0.5
start_1=44
end_1=47
start_2=108
end_2=111

start_3=69
end_3=78
start_4=109
end_4=118


trap.write(str(start_1)+" "+str(end_2)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(start_2)+" "+str(end_1)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.write(str(start_3)+" "+str(end_4)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")
trap.write(str(start_4)+" "+str(end_3)+" "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()
 
