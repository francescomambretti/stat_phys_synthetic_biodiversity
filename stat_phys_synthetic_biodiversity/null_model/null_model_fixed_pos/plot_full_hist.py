# Python script by Francesco Mambretti
# load and plot marginal probability distributions for species abundance
# each plot represents the distribution of R simulations, each of which yields l+1 values for the relative species abundance (which coincides with the number of strings having a given overlap a)
# OK version - 26/05/2021 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
	
l=20 #maximum possible consecutive overlap, i.e. maximum possible species index

full_hist=np.loadtxt("full_hist.txt",unpack=True)
full_hist=full_hist.T

plt.rcParams['figure.figsize'] = [20, 20]
cmap = cm.brg
figure, axis = plt.subplots(7,3)

weights = np.ones_like(full_hist[0]) / len(full_hist[0]) # in order to normalize histogram
norm=float(l+1) #l + species 0
nbins=25

for j in range (0,3):
	for i in range (0,7):
		index=i*3+j
		mean=np.mean(full_hist[index])
		if(mean!=0):
			axis[i][j].hist(full_hist[index]/mean,color=cmap(float(index)/float(l)),edgecolor='black', alpha=1,label='a='+str(index),bins=nbins,weights=weights,align='left')
		else:
			axis[i][j].hist(full_hist[index],color=cmap(float(index)/float(l)),edgecolor='black', alpha=1,label='a='+str(index),bins=nbins,weights=weights,align='left')
        	
		axis[i][j].legend(fontsize=14)
		plt.setp(axis[i][j].get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor",fontsize='14')
		plt.setp(axis[i][j].get_yticklabels(), fontsize='14')
		axis[i][j].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
		axis[i][j].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        
plt.tight_layout()
plt.savefig("all_marginals.png",dpi=150)

