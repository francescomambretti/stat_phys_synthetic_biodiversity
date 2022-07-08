## Written by Francesco Mambretti
## 14/08/2021
## To analyze Genetic Algorithm results

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.pyplot import cm

#### define variables ####

gens=30 #generations of Genetic Algorithm
l=20  #maximum allowed overlap
max_seed=1 #number of independent runs

colorlist=('firebrick','orangered','gold','lime','green','steelblue','gray')

height_in_time=np.zeros((gens,l+1))

### average data over the various seeds
for g in range (0,gens,1):
    ave_histo=np.zeros(l+1)
    ave_histo2=np.zeros(l+1)
    for s in range(1,max_seed+1,1):
        data=np.loadtxt("alpha_0.25/seed{}/histogram{}.dat".format(s,g),unpack=True)
        weights = np.ones_like(data)/len(data)
        histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
        ave_histo+=histo
        ave_histo2+=histo*histo

        height_in_time[g]+=histo

	#to plot histograms
    ave_histo/=max_seed
    ave_histo2/=max_seed
    error=np.sqrt((ave_histo2-ave_histo*ave_histo)/max_seed)
    
    if (g%5==0):
        plt.bar(np.arange(0,l+1,1),ave_histo,label='gen. {}'.format(g),alpha=1.0,align='center',color=colorlist[int(g/5)],yerr=error,capsize=4,edgecolor='black',ecolor=colorlist[int(g/5)])

plt.xlabel("Maximum consecutive overlap - species index",fontsize=12)
plt.xticks(np.arange(0,l+2,4),fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=13)
plt.title("$FS$ criterion - $N_{res}=1E4$")
plt.savefig('ave_histo_RSA_FS_1E4_1E6.png',dpi=300)

plt.yscale('log')
plt.savefig('ave_histo_RSA_FS_1E4_1E6_log.png',dpi=300)

#plot heights of bins in time
plt.clf()
plt.figure(figsize=(12,9))

#normalize
height_in_time/=max_seed
#transpose
height_in_time=height_in_time.T

color = iter(cm.brg(np.linspace(0, 1, 11)))

for omega in range (4,14,1):
    c=next(color)
    plt.plot(height_in_time[omega],marker='o',linestyle='dashed', color=c,label='$\omega={}$'.format(omega))
    
plt.xlabel("Cycles",fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
#plt.yscale('log')
plt.legend(fontsize=15)
plt.savefig('omega_in_time.png',dpi=300)
