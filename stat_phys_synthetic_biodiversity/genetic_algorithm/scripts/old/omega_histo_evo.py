## Written by Francesco Mambretti
## 14/08/2021
## To analyze Genetic Algorithm results and plot histo(omega)
## 11/07/2022  version

import matplotlib.pyplot as plt
import numpy as np
import os

#### define variables ####

gens=30 #generations of Genetic Algorithm
l=20  #maximum allowed overlap
max_seed=4 #number of independent runs

colorlist=('firebrick','orangered','gold','lime','green','steelblue','gray')

### average data over the various seeds
for g in (1,12,24):#range (0,gens+1,5):
    ave_histo=np.zeros(l+1)
    ave_histo2=np.zeros(l+1)
    for s in range(0,max_seed,1):
        histo=np.loadtxt("../fitness_FS/seed{}/histogram{}.dat".format(s,g),unpack=True)
        histo/=np.sum(histo)
        ave_histo+=histo
        ave_histo2+=histo*histo

    #to plot histograms
    ave_histo/=max_seed
    ave_histo2/=max_seed
    error=np.sqrt((ave_histo2-ave_histo*ave_histo)/max_seed)

    plt.bar(np.arange(0,l+1,1),ave_histo,label='gen. {}'.format(g),alpha=1.0,align='center',color=colorlist[int(g/5)],yerr=error,capsize=4,edgecolor='black',ecolor=colorlist[int(g/5)])

plt.xlabel("$\omega$",fontsize=12)
plt.xticks(np.arange(0,l+2,4),fontsize=13)
plt.yticks(fontsize=13)
plt.legend(fontsize=13)
plt.title("$f = 10 \\times \Theta(10-\omega) + \\dfrac{\omega}{\Sigma \omega} \\times \Theta(\omega-10) $ - $N_{p}=1E5, N_{r}=1E4$")
plt.savefig("../fitness_FS/ave_histo_omega_F_1E4_1E5.png",dpi=300)

plt.yscale('log')
plt.savefig("../fitness_FS/ave_histo_omega_F_1E4_1E5_log.png",dpi=300)
