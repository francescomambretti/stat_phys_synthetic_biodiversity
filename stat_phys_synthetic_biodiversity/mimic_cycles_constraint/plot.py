# Plot the data obtained in mimic_cycles_constraint
# Written by Francesco Mambretti, 28/06/2021

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
    print(c+1)

    ## do the plot ##
    plt.bar(np.arange(0,l+1,1),height=cycle,color=colorlist[index],edgecolor='black',alpha=1.0,width=widthlist[index],label="cycle #{}".format(c+1))
    plt.errorbar(np.arange(0,l+1,1), cycle, errors, fmt="o", color='black', capsize=5, markersize=0)
    index+=1

plt.xticks(np.arange(0,l+1,1),fontsize=11)

plt.xlabel("$a$", fontsize=13)
plt.ylabel('bin occupancy probability'.format(c), fontsize=13)

plt.yscale('log')

plt.legend(fontsize=14)

#plt.show()
plt.savefig('mimic_cycles_constraint.png',dpi=300)
os.system("open mimic_cycles_constraint.png")
