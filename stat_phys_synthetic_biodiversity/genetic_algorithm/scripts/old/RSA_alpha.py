## Written by Francesco Mambretti
## 04/09/2021
## To analyze the effect of the alpha parameter on the RSA histogram after the GA simulation
## plot the RSA histogram for a single chosen cycle of the algorithm and for all the alpha values specified

import matplotlib.pyplot as plt
import numpy as np
import os

#### define variables ####

gens=30 #generations of Genetic Algorithm
l=20  #maximum allowed overlap
max_seed=10 #number of independent runs

colorlist=('black','blue','gold','seagreen','purple') #one for each alpha

### average data over the various seeds
for alpha in np.arange (0,1.25,0.25): #loop over alpha
	ave_histo=np.zeros(l+1)
	ave_histo2=np.zeros(l+1)
	for s in range(1,max_seed,1):
		data=np.loadtxt("./alpha_{}/asymp/seed{}/histogram30.dat".format(alpha,s),unpack=True)
		weights = np.ones_like(data)/len(data)
		histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
		ave_histo+=histo
		ave_histo2+=histo*histo

	#to plot histograms
	ave_histo/=max_seed
	ave_histo2/=max_seed
	error=np.sqrt((ave_histo2-ave_histo*ave_histo)/max_seed)
	
	plt.bar(np.arange(0,l+1,1),ave_histo,label='$\alpha=$ {}'.format(alpha),alpha=1.0,align='center',color=colorlist[int(alpha*4)],yerr=error,capsize=4,edgecolor='black',ecolor=colorlist[int(alpha*4)])

plt.xlabel("Maximum consecutive overlap - species index",fontsize=12)
plt.xticks(np.arange(0,l+2,4),fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=13)
plt.title("$FS$ criterion - $N_{res}=5E4$")
plt.yscale('log')
plt.savefig('ave_histo_RSA_change_alpha_log.png',dpi=300)
