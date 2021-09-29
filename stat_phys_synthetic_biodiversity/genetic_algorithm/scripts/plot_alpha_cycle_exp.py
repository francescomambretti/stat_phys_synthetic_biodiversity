# Written by Francesco Mambretti - 07/09/2021

# Edit 28/09/2021: now it plots the KL div and not the euclidean distance, commented option
# Edit 29/09/2021: added the entropy plot as a function of cycle

import numpy as np
import matplotlib.pyplot as plt

listA=("black", "brown", "red", "tomato" , "dimgray")
listB=("aquamarine", "seagreen", "lime", "green", "gold")

alpha_list=(0.0,0.25,0.5,0.75,1.0)

index=-1

for a in alpha_list:
	index+=1
	#cost_euclid1=np.loadtxt("alpha_{}/max/eucl_dist".format(a),unpack=True)
	#cost_euclid2=np.loadtxt("alpha_{}/asymp/eucl_dist".format(a),unpack=True)

	cost_kl1=np.loadtxt("alpha_{}/max/kl_div".format(a),unpack=True)
	cost_kl2=np.loadtxt("alpha_{}/asymp/kl_div".format(a),unpack=True)

	plt.plot(cost_kl1,marker='P',linestyle='none',color=listA[index],label="$\\alpha=${}, $F$".format(a))
	plt.plot(cost_kl2,marker='D',linestyle='none',color=listA[index],label="$\\alpha=${}, $FS$".format(a))

plt.xlabel("GA cycle",fontsize=12)
plt.ylabel("KL divergence",fontsize=12)
plt.legend(ncol=2,fontsize=10)
plt.savefig('KL_divergences.png',dpi=300)
plt.clf()

########## entropy ###########

index=-1

for a in alpha_list:
	index+=1

	entropy1=np.loadtxt("alpha_{}/max/entropy".format(a),unpack=True)
	entropy2=np.loadtxt("alpha_{}/asymp/entropy".format(a),unpack=True)

	plt.plot(entropy1,marker='o',linestyle='none',color=listA[index],label="$\\alpha=${}, $F$".format(a))
	plt.plot(entropy2,marker='s',linestyle='none',color=listA[index],label="$\\alpha=${}, $FS$".format(a))

plt.xlabel("GA cycle",fontsize=12)
plt.ylabel("entropy",fontsize=12)
plt.legend(ncol=2,fontsize=10)
plt.savefig('entropy.png',dpi=300)
plt.clf()

'''index=-1
for a in alpha_list:
	index+=1
	cost_ks1=np.loadtxt("alpha_{}/max/ks_dist".format(a),unpack=True,usecols=(1,))
	cost_ks2=np.loadtxt("alpha_{}/asymp/ks_dist".format(a),unpack=True,usecols=(1,))

	plt.plot(cost_ks1,marker='o',color=listA[index],label="$\\alpha=${}, $F$".format(a))
	plt.plot(cost_ks2,marker='s',color=listA[index],label="$\\alpha=${}, $FS$".format(a))

plt.xlabel("GA cycle",fontsize=12)
plt.legend(ncol=2,fontsize=10)
plt.savefig('KS_distances.png',dpi=300)
'''
