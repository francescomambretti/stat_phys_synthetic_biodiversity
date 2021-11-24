### Written by Francesco Mambretti, 21/10/2021
###
### Useful tools for mixed_hb_pair_evo.py code

### 21/10/2021 version

import numpy as np
import sys
import os

### set parameters
CTRL=100000 # periodically print sanity check
l=20 # resource
L=50 # predator

##############################################################

def set_watershed(n_pred,n_res,l,L):

    if (n_res==0 and n_pred==2):
        watershed=L-1
    elif (n_res==1 and n_pred==1):
        watershed=n_pred*L-1
    elif(n_res==2 and n_pred==0):
        watershed=l-1

    return watershed
    
##############################################################

def belongs_to(index, s_low, s_top): #check whether a nucleotide belongs to a strand
    if ( index >= s_low and index < s_top):
        return True
    else:
        return False

##############################################################

def read_and_process (mother_file,output_folder,output_file_name,n_pred,n_res,i,j):

    step=0
    tot=0
    mco=0
    temp_mco=0
    prec_p=-1
    prec_r=-1
    count_line=0

    pred_offset=L*n_pred-1 #index of the last nucleotide belonging to the predator with the highest strand index

    # I am comparing strands i and j
    if (i<n_pred): #it is a predator
        low_i=i*L #index of the first nucleotide of i
        top_i=(i+1)*L
    
    else: #it is a resource
        low_i=pred_offset+1+(i-n_pred)*l #index of the first nucleotide of i
        top_i=pred_offset+1+(i-n_pred+1)*l
        
    if (j<n_pred): #it is a predator
        low_j=j*L #index of the first nucleotide of i
        top_j=(j+1)*L
    
    else: #it is a resource
        low_j=pred_offset+1+(j-n_pred)*l #index of the first nucleotide of i
        top_j=pred_offset+1+(j-n_pred+1)*l

    with open(mother_file,"r") as f: #read file with the full hb list
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
                f = open(output_folder+"/"+str(i)+"_"+str(j)+"_"+output_file_name, "a")
                f.write(str(step)+' '+str(tot)+' '+str(mco)+'\n')
                f.close()

            else: #equivalently, elif (line[0]!="\n" and line[0]!="#"):
                p = int(line.split()[0]),int(line.split()[1])
                #print(p[0],p[1],low_i,top_i,low_j,top_j,i,j,belongs_to(p[1], low_j, top_j))
                
                if ( (belongs_to(p[0], low_i, top_i) and belongs_to(p[1], low_j, top_j)) or ( ( belongs_to(p[1], low_i, top_i) and belongs_to(p[0], low_j, top_j)) ) ): # I am currently not interested in the self-interactions of a strand with itself ; moreover, p[0] always < p[1] so the second condition after "or" should never be met
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

######################################################

def Error (a,a_2,norm):
    return np.sqrt ((a_2/norm-(a/norm)*(a/norm))/norm)

######################################################

def accumulate(ave_tot_hb,ave_mco_hb,appo_tot_hb,appo_mco_hb,ave_tot_hb_2,ave_mco_hb_2):
    ave_tot_hb+=appo_tot_hb
    ave_mco_hb+=appo_mco_hb
    ave_tot_hb_2+=appo_tot_hb*appo_tot_hb
    ave_mco_hb_2+=appo_mco_hb*appo_mco_hb
    
    return

######################################################

def initialize(ave_tot_hb,ave_mco_hb,appo_tot_hb,appo_mco_hb):
    ave_tot_hb=appo_tot_hb
    ave_mco_hb=appo_mco_hb
    ave_tot_hb_2=appo_tot_hb*appo_tot_hb
    ave_mco_hb_2=appo_mco_hb*appo_mco_hb
    
    return ave_tot_hb, ave_mco_hb, ave_tot_hb_2, ave_mco_hb_2

######################################################

def reset(ave_tot_hb,ave_mco_hb,all_tot_hb,all_mco_hb,ave_tot_hb_2,ave_mco_hb_2):
    ave_tot_hb=np.zeros(0)
    ave_mco_hb=np.zeros(0)
    all_tot_hb=np.zeros(0)
    all_mco_hb=np.zeros(0)
    ave_tot_hb_2=np.zeros(0)
    ave_mco_hb_2=np.zeros(0)
    
    return ave_tot_hb,ave_mco_hb,all_tot_hb,all_mco_hb,ave_tot_hb_2,ave_mco_hb_2
