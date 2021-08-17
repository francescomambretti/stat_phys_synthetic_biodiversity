## Written by Francesco Mambretti
## 15/08/2021
## To compare Genetic Algorithm results and the null model
## Plot the outcome of cycle 0 for the RSA histogram and the null model (without threshold) prediction

import matplotlib.pyplot as plt
import numpy as np
import os

#### define variables ####

l=20  #maximum allowed overlap
max_seed=10 #number of independent runs

### average data over the various seeds

for g in range (29,30,1):  # modify according to the need
	ave_histoF=np.zeros(l+1)
	ave_histo2F=np.zeros(l+1)
	ave_histoFS=np.zeros(l+1)
	ave_histo2FS=np.zeros(l+1)
	ave_histoR=np.zeros(l+1)
	ave_histo2R=np.zeros(l+1)
	for s in range(1,max_seed,1):
		data=np.loadtxt("max/1e5_1e7/seed{}/histogram{}.dat".format(s,g),unpack=True)
		weights = np.ones_like(data)/len(data)
		histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
		ave_histoF+=histo
		ave_histo2F+=histo*histo

		data=np.loadtxt("asymp/1e5_1e7/seed{}/histogram{}.dat".format(s,g),unpack=True)
		histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
		ave_histoFS+=histo
		ave_histo2FS+=histo*histo

		data=np.loadtxt("constant/1e5_1e7/seed{}/histogram{}.dat".format(s,g),unpack=True)
		histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
		ave_histoR+=histo
		ave_histo2R+=histo*histo
	

	#to plot histograms
	ave_histoF/=max_seed
	ave_histo2F/=max_seed
	errorF=np.sqrt((ave_histo2F-ave_histoF*ave_histoF)/max_seed)
	ave_histoFS/=max_seed
	ave_histo2FS/=max_seed
	errorFS=np.sqrt((ave_histo2FS-ave_histoFS*ave_histoFS)/max_seed)
	ave_histoR/=max_seed
	ave_histo2R/=max_seed
	errorR=np.sqrt((ave_histo2R-ave_histoR*ave_histoR)/max_seed)

	#load null model data
	null=np.loadtxt("null_model.dat",unpack=True)
	null_T=np.loadtxt("null_model_threshold3.dat",unpack=True)

	plt.bar(np.arange(0,l+1,1),ave_histoF,label='gen.30 - selection criterion $F$'.format(g),alpha=0.85,align='center',color='steelblue',yerr=errorF,capsize=4,edgecolor='black',ecolor='steelblue',width=1)

	plt.bar(np.arange(0,l+1,1),ave_histoFS,label='gen.30 - selection criterion $FS$'.format(g),alpha=0.85,align='center',color='green',yerr=errorFS,capsize=4,edgecolor='black',ecolor='green',width=0.8)

	plt.bar(np.arange(0,l+1,1),ave_histoR,label='gen.30 - selection criterion $R$'.format(g),alpha=0.85,align='center',color='crimson',yerr=errorR,capsize=4,edgecolor='black',ecolor='crimson',width=0.6)

	plt.plot(np.arange(0,l+1,1),null,color='black',markersize=6,marker='o',linestyle='none',label='null model')
	plt.plot(np.arange(0,l+1,1),null_T,color='gray',markersize=6,marker='s',linestyle='none',label='null model with $T=3$')	

	plt.yscale('log')
	plt.xlabel("Maximum consecutive overlap - species index",fontsize=12)
	plt.xticks(np.arange(0,l+2,4),fontsize=13)
	plt.yticks(fontsize=13)
	plt.legend(fontsize=11)
plt.savefig('compare_ga_null_1e5_1e7.png',dpi=300)
