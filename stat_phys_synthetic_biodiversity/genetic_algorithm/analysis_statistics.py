## Written by Francesco Mambretti
## 14/08/2021
## To analyze Genetic Algorithm results

import matplotlib.pyplot as plt
import numpy as np
import os

#### define variables ####

gens=30 #generations of Genetic Algorithm
l=20  #maximum allowed overlap
max_seed=10 #number of independent runs

colorlist=('firebrick','orangered','gold','lime','green','steelblue')

### average data over the various seeds
for g in range (0,gens,5):
	ave_histo=np.zeros(l+1)
	ave_histo2=np.zeros(l+1)
	for s in range(1,max_seed,1):
		data=np.loadtxt("./seed{}/histogram{}.dat".format(s,g),unpack=True)
		weights = np.ones_like(data)/len(data)
		histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
		ave_histo+=histo
		ave_histo2+=histo*histo

	#to plot histograms
	ave_histo/=max_seed
	ave_histo2/=max_seed
	error=np.sqrt((ave_histo2-ave_histo*ave_histo)/max_seed)
	
	plt.bar(np.arange(0,l+1,1),ave_histo,label='gen. {}'.format(g),alpha=1.0,align='center',color=colorlist[int(g/5)],yerr=error,capsize=4,edgecolor='black',ecolor=colorlist[int(g/5)])

	#plt.yscale('log')
	plt.xlabel("Maximum consecutive overlap - species index",fontsize=12)
	plt.xticks(np.arange(0,l+2,4),fontsize=13)
	plt.yticks(fontsize=13)
	plt.legend(fontsize=13)
plt.savefig('ave_histo_RSA.png',dpi=300)
