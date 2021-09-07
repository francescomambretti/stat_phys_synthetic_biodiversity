## Written by Francesco Mambretti
## 04/09/2021
## Plot the RMS distance between the outcome of each cycle for the RSA histogram and the experimental data (at a given cycle)
## KS test added
## 07/09/2021 version

import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import kstest, ks_2samp

#### define variables ####

l=20  #maximum allowed overlap
max_seed=10 #number of independent runs
gens=30

#load experimental data
print("load experimental data")

experim=np.loadtxt("experim_data_cycle9.dat",unpack=True)
experim/=np.sum(experim) #normalize data

alpha_list=(0.0,0.25,0.5,0.75,1.0)

### average data over the various seeds
for a in alpha_list:
	print("alpha={}".format(a))
	d1=np.zeros(0)
	d2=np.zeros(0)
	ks1=np.zeros((0,2))
	ks2=np.zeros((0,2))
	for g in range (0,gens+1,1):  # modify according to the need
	#	print(g)
		ave_histoF=np.zeros(l+1)
		ave_histo2F=np.zeros(l+1)
		ave_histoFS=np.zeros(l+1)
		ave_histo2FS=np.zeros(l+1)
		#ave_histoR=np.zeros(l+1)
		#ave_histo2R=np.zeros(l+1)
		for s in range(1,max_seed,1):
			data=np.loadtxt("alpha_{}/max/seed{}/histogram{}.dat".format(a,s,g),unpack=True)
			weights = np.ones_like(data)/len(data)
			histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
			ave_histoF+=histo
			ave_histo2F+=histo*histo

			data=np.loadtxt("alpha_{}/asymp/seed{}/histogram{}.dat".format(a,s,g),unpack=True)
			histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
			ave_histoFS+=histo
			ave_histo2FS+=histo*histo

		#	data=np.loadtxt("constant/seed{}/histogram{}.dat".format(s,g),unpack=True)
		#	histo,binEdges=np.histogram(data,bins=np.arange(0,l+2,1),weights=weights)
		#	ave_histoR+=histo
		#	ave_histo2R+=histo*histo
	

		#to plot histograms
		ave_histoF/=max_seed
		ave_histo2F/=max_seed
		errorF=np.sqrt((ave_histo2F-ave_histoF*ave_histoF)/max_seed)
		ave_histoFS/=max_seed
		ave_histo2FS/=max_seed
		errorFS=np.sqrt((ave_histo2FS-ave_histoFS*ave_histoFS)/max_seed)
		#ave_histoR/=max_seed
		#ave_histo2R/=max_seed
		#errorR=np.sqrt((ave_histo2R-ave_histoR*ave_histoR)/max_seed)

		#plt.clf()

		np.savetxt("alpha_{}/max/ave_histoF_cycle{}.dat".format(a,g),ave_histoF)
		np.savetxt("alpha_{}/asymp/ave_histoFS_cycle{}.dat".format(a,g),ave_histoFS)

		d1=np.append(d1,np.sqrt(np.sum((experim-ave_histoF)**2)))
		d2=np.append(d2,np.sqrt(np.sum((experim-ave_histoFS)**2)))
		ks1_s,ks1_p=ks_2samp(experim, ave_histoF)
		ks2_s,ks2_p=ks_2samp(experim, ave_histoFS)

		ks1=np.vstack((ks1,[ks1_s,ks1_p]))
		ks1=np.vstack((ks2,[ks2_s,ks2_p]))

#		print("For alpha={}, cycle {} yields an euclidean distance of the RSA histogram equal to:".format(a,g))
#		print("{:.5} with $F$ criterion and of {:.5} with $FS$ criterion; KS test yields: {:.5},{:.5} and {:.5},{:.5}".format(d1,d2,ks1_s,ks1_p,ks2_s,ks2_p))


		np.savetxt("alpha_{}/max/eucl_dist".format(a), d1)
		np.savetxt("alpha_{}/asymp/eucl_dist".format(a),d2)
		np.savetxt("alpha_{}/max/ks_dist".format(a), ks1)
		np.savetxt("alpha_{}/asymp/ks_dist".format(a),ks2)
		

		#plt.plot(np.arange(0,l+1,1),np.sqrt((experim-ave_histoF)*(experim-ave_histoF)),color='steelblue',markersize=6,marker='o',linestyle='none',label='RMS error; F criterion')
		#plt.plot(np.arange(0,l+1,1),np.sqrt((experim-ave_histoFS)*(experim-ave_histoFS)),color='green',markersize=6,marker='o',linestyle='none',label='RMS error; FS criterion')

		#plt.yscale('log')
		#plt.xlabel("Maximum consecutive overlap - species index",fontsize=12)
		#plt.xticks(np.arange(0,l+2,4),fontsize=13)
		#plt.yticks(fontsize=13)
		#plt.legend(fontsize=11)
		#plt.savefig('alpha_{}/compare_ga_exp_5e4_5e6_{}.png'.format(a,g),dpi=200)
