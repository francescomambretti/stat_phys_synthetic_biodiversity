# Plot the data obtained in mimic_cycles_constraint, divided by the corresponding experimental data
# If the ratio is around 1, the model is predictive for that bin occupancy
#
# Written by Francesco Mambretti, 13/07/2021

import numpy as np
import os
import matplotlib.pyplot as plt

#define 

l=20
L=50
bases=4
index=0
colorlist=('gold','red','blue')
widthlist=(0.9,0.7,0.5)

##load simulated data for comparison
for c in range(2,9,3):
    cycle,errors=np.loadtxt("cycle.{}.txt".format(c+1),unpack=True,usecols=(1,2))
    experim=np.loadtxt("../experim_data/analysis/R{}-06_ready/histogram.txt".format(c+1),unpack=True) #should be changed with the right path to the file
    experim/=np.sum(experim)
    print(cycle,experim)

    ## do the plot ##
    plt.plot(cycle/experim,marker='o',linestyle='none',color=colorlist[index],label="cycle #{}".format(c+1))
    plt.errorbar(np.arange(0,l+1,1), cycle/experim, errors/experim, fmt="o", color='black', capsize=5, markersize=0)
    index+=1

plt.xticks(np.arange(0,l+1,1),fontsize=11)

plt.xlabel("$a$", fontsize=13)
plt.ylabel('mimicked/experimental probability'.format(c), fontsize=13)

plt.title("Mimicked vs experim. protocol 06") # change, in case, the protocol name

plt.yscale('log')

plt.legend(fontsize=14)

plt.savefig('mimic_vs_experim_06.png',dpi=300)
os.system("open mimic_vs_experim_06.png")
