### Written by Francesco Mambretti, 21/10/2021
###
### Based on mixed_hb_pair.py - generalize to count the interactions between any pair of strands in the system
### Counts the total HB and the MCO between all strand pairs as a function of time. The code creates a separate file for each pair
### Also self-interactions of each strand with itself are now included
###
### The code can analyze either a single run or a set of simulations
###
### A file of input options must be provided
###
### This version also includes the HBs between fixed sequences at the ends of predators and blocking caps

### 11/04/2022 version

import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import functions_mixed_hb as func
from input_params import *  #external file where input params are defined

######################################## initialization ###################################################

tot_strands=n_pred+n_res+caps_per_pred*n_pred # +caps_per_pred*n_pred --> caps_per_pred blocking caps for each predator
ave_tmo_hb=np.zeros(0)
ave_mco_hb=np.zeros(0)
all_tmo_hb=np.zeros(0)
all_mco_hb=np.zeros(0)
ave_tmo_hb_2=np.zeros(0)
ave_mco_hb_2=np.zeros(0)

######################################## HB analysis ###################################################

# analyze each simulation independently
for r in range (0,runs):
    mother_folder="run_{}/MD".format(r)
    if not os.path.exists(mother_folder+"/"+output_folder):
        os.makedirs(mother_folder+"/"+output_folder)
    # loop throughout all the strand pairs
    print("Analyze run {}...".format(r))
    for i in range(0,tot_strands): #each strand is compared with all the others, itself included
        for j in range (i,tot_strands):
            print(i,j)
            #create the corresponding files
            os.system("rm "+mother_folder+"/"+output_folder+"/"+str(i)+"_"+str(j)+"_"+output_file_name)
                
            #watershed=functions_mixed_hb.set_watershed(n_pred,n_res,l,L)
            #print(watershed) # can be uncommented
            func.read_and_process(mother_folder+"/"+hb_list_file,mother_folder+"/"+output_folder,
                                    output_file_name,n_pred, n_res,i,j,tot_strands) #iterate through all hb_list_file

######################################## plots ###################################################

print("Plot...")
#loop over folders
for a in range (0,tot_strands):
    for b in range(a,tot_strands):
        #set to 0 all values
        ave_tmo_hb,ave_mco_hb,all_tmo_hb,all_mco_hb,ave_tmo_hb_2,ave_mco_hb_2=func.reset(ave_tmo_hb,ave_mco_hb,all_tmo_hb,all_mco_hb,ave_tmo_hb_2,ave_mco_hb_2)

        for r in range (0,runs):
            print (a,b,r)
            mother_folder="run_{}/MD".format(r)
            timesteps,appo_tmo_hb,appo_mco_hb=np.loadtxt(mother_folder+"/"+output_folder+"/"+str(a)+"_"+str(b)+"_"+output_file_name, usecols=(0,1,2), unpack=True)
            
            if (r!=0):
                func.accumulate(ave_tmo_hb,ave_mco_hb,appo_tmo_hb,appo_mco_hb,ave_tmo_hb_2,ave_mco_hb_2)

            else:
                ave_tmo_hb, ave_mco_hb, ave_tmo_hb_2, ave_mco_hb_2=func.initialize(ave_tmo_hb,ave_mco_hb,appo_tmo_hb,appo_mco_hb)

            all_tmo_hb=np.append(all_tmo_hb,appo_tmo_hb)
            all_mco_hb=np.append(all_mco_hb,appo_mco_hb)

        #plot

        freq=100  #do not plot all the data, just 1 every freq

        plt.errorbar(x=timesteps[::freq],y=ave_tmo_hb[::freq]/runs,yerr=func.Error(ave_tmo_hb[::freq],ave_tmo_hb_2[::freq],runs),marker='o',markersize=4,color='darkorange',label='ave total mixed HB',capsize=4)
        plt.errorbar(x=timesteps[::freq],y=ave_mco_hb[::freq]/runs,yerr=func.Error(ave_mco_hb[::freq],ave_mco_hb_2[::freq],runs),marker='o',markersize=4,color='royalblue',label='ave MCO HB',capsize=4)

        plt.title(members+" ave over {} runs, strand {} vs {}".format(runs,a,b), fontsize=16)

        plt.xlabel("$1000 \\times $ MD steps ",fontsize=14)
        plt.ylabel("Overlap",fontsize=14)

        plt.legend(fontsize=14)
        plt.savefig('ave_HB_analysis_{}_{}.png'.format(a,b),dpi=300)

        plt.clf()
        plt.hist(all_tmo_hb,bins=np.arange(0,51,1),rwidth=0.9,density=True,align='left',color='darkorange',edgecolor='black',label='TMO HB')
        plt.hist(all_mco_hb,bins=np.arange(0,51,1),rwidth=0.7,alpha=1.0,density=True,align='left',color='royalblue',edgecolor='black',label='MCO HB')

        hist, bin_edges = np.histogram(all_mco_hb,bins=np.arange(0,51,1), density=True)
        hist2, bin_edges2 = np.histogram(all_tmo_hb,bins=np.arange(0,51,1), density=True)
        np.savetxt("all_HB_histo_{}_{}.txt".format(a,b), np.column_stack((hist, hist2)), fmt='%1.2f')

        plt.xlabel("Overlap",fontsize=14)
        plt.xlim(-1,30)
        plt.xticks(np.arange(0,31,3))

        plt.title(members+" ave over {} runs, strand {} vs {}".format(runs,a,b),fontsize=16)

        plt.legend(fontsize=14)
        plt.savefig('all_HB_histo_{}_{}.png'.format(a,b),dpi=300)
        plt.clf()
