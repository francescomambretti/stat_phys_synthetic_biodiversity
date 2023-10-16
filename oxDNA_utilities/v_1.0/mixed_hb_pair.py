### Written by Francesco Mambretti, 27/9/2021
###
### Counts the total HB and the MCO between 2 strands: either 1 predator + 1 resource, or 1+1 predators or 1+1 resources as a function of time

### 21/10/2021 version

import numpy as np
import sys
import os

############################################################## initialization ##############################################################

if len(sys.argv) < 5:
    print("please insert the hb_list file, the number of predators, the number of resources and the output file name")
    sys.exit()

n_pred=int(sys.argv[2])
n_res=int(sys.argv[3])
output_file_name=sys.argv[4]

CTRL=100000 # periodically print sanity check
l=20 # resource
L=50 # predator

if (n_res==0 and n_pred==2):
    watershed=L-1
else if (n_res==1 and n_pred==1):
    watershed=n_pred*L-1
else if(n_res==2 and n_pred==0):
        watershed=l-1

#print(watershed) # can be uncommented

#set variables
mco=0
temp_mco=0
prec_p=-1 #track index of previous predator nucleotide
prec_r=-1 #track index of previous resource nucleotide
tot=0
count_line=0
step=0

os.system("rm "+output_file_name)

############################################################## read and process ##############################################################


with open(sys.argv[1],"r") as f:
    for line in f:
        #if (count_line%CTRL==0): #uncomment to know that it is doing something
        #print("line: ", count_line)
        #print(line,line[0],tot)

        if (line[0]=="#"): # a new step starts
            step+=1
            tot=0
            mco=0
            temp_mco=0
            prec_p=-1
            prec_r=-1
    
        elif(line[0]=="\n"):
            if (temp_mco>mco):
                mco=temp_mco
            #the timestep is over, print file
            f = open(output_file_name, "a")
            f.write(str(step)+' '+str(tot)+' '+str(mco)+'\n')
            f.close()
            
        elif(line[0]!="#" and line[0]!="\n"): #equivalently, elif (line[0]!="\n" and line[0]!="#"):
            p = int(line.split()[0]),int(line.split()[1])
            if (p[0]<=watershed and p[1]>watershed): #p[0] always < p[1]
                tot+=1
                if (prec_p==-1 and prec_r==-1):
                    temp_mco=1 # we have a pair of matching bases!
                    mco=1
                    prec_p=p[0] #save indexes for next loop
                    prec_r=p[1]

                elif (prec_p!=-1 and prec_r!=-1): #they are initialized to some nucleotide index
                    if (np.abs(p[0]-prec_p)==1 and np.abs(prec_r-p[1])==1):  #N.B.: this works only if hb_list.dat's second column is ordered in ascending order and the first one is not, otherwise one needs to use np.abs
                       temp_mco+=1

                    else:
                        # equivalently - if (np.abs(p[0]-prec_p)>1 or np.abs(prec_r-p[1])>1): #there is a gap, a non-consecutive HB or similar
                        #beginning of a new set of possibly overlapping bases, preceded by other non-consecutive overlapping bases
                        if (temp_mco>mco):
                            mco=temp_mco
                        temp_mco=1

                prec_p=p[0]
                prec_r=p[1]
                


        count_line = count_line+1    
