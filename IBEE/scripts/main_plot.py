# Python code to plot Evolutive Algorithm results
# Written by Francesco Mambretti, 11/07/2022
#
# Checkpoints for performance monitoring are currently commented
# Strongly modified wrt previous version --> essentialy no quantities are computed here, they are all computed in analysis.py and then here they are only plotted!
#
# 25/07/2022 version

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from mpl_toolkits.mplot3d import axes3d
from matplotlib.pyplot import cm
import module_plots as mop
from input_params import *
import module_functions as func
from scipy.optimize import curve_fit

################  variables ################

maximum=l+1 #go from 0 to l

list_files=["ShEn_omega","ShEn_RSA","omega_time","unique_sequences","top_{}_percentage_over_cycles".format(n),"zip ratio"]
list_mean_files=["mean_ShEn_omega","mean_ShEn_RSA","mean_omega_time","mean_unique_sequences","mean_top_{}_percentage_over_cycles".format(n),"mean_zip_ratio"]

################  plots ################

for s in range (1,max_seed+1,1):

    c=-1 #index for cycles where I plot p(omega) histograms
    seed_folder=mother_folder+"/seed{}/".format(s)
    for t in cycles_IDs:#loop through cycles and cumulate data from parallel simulations
        #useful variables
        cycle_folder=mother_folder+"/seed{}/cycle_{}/".format(s,t)
        
        if (t in cycles_plot_h_omega): #separate plots for each seed
            c+=1
            mop.plot_histo_omega(t,c)
        
        #MCO cross-matrix among top-n
        if (t%freq_cross_MCO==0):
            mop.plot_cross_MCOs(cycle_folder+"/cross_MCOs.txt",cycle_folder+"/cross_MCOs.png",t,True)
            #do it also for averaged!
            mop.plot_cross_MCOs(results_folder+"/cycle_{}/mean_cross_MCOs.txt".format(t),results_folder+"/cycle_{}/mean_cross_MCOs.png".format(t),t,False)
            
        if (t==max_cycle-1):
            # plots regarding dominant strands at last cycle
            mop.top_n_plot(n,seed_folder) #plot the time evolution of the dominants of the last cycle only
            mop.top_n_fraction(n,seed_folder,False) #plot the fraction of the total population they occupy
            #do it also for averaged!
            mop.top_n_fraction(n,results_folder,True)

    mop.plot_ShEn(seed_folder+"/"+list_files[i_shen_omega]+".txt",seed_folder+"/"+list_files[i_shen_omega]+".png",'navy')
    mop.plot_ShEn(seed_folder+"/"+list_files[i_shen_RSA]+".txt",seed_folder+"/"+list_files[i_shen_RSA]+".png",'darkturquoise')


################ plot quantities avaraged over parallel simulations ################

list_mean_labels=["mean Shannon Entropy ($\omega$)","mean Shannon Entropy (RSA)","$< \omega >$","% unique sequences","mean top-{}".format(n), "zip ratio"]
colors=["navy","darkturquoise","firebrick","gray","forestgreen","black"]

for o in range (0,n_obs):
    # plot mean Shannon Entropies, mean omega of the population and % of unique sequences as a function of time
    mop.plot_obs_ave_parall(results_folder+"/"+list_mean_files[o]+".txt",results_folder+"/"+list_mean_files[o]+".png",list_mean_labels[o],colors[o])

'''
# fit decay
def expo_decay(a,x):
    return 100*np.exp(-a*x)
    
def power_law(a,b,x):
    return b*np.power(x,a)

ydata=np.loadtxt(results_folder+"/unique_sequences.txt",unpack=True,usecols=(0,))
xdata=np.arange(0,len(ydata),1)
popt, pcov = curve_fit(expo_decay, xdata, ydata)
plt.plot(xdata, expo_decay(xdata, *popt),label='Evo Alg model, expo fit: a=%5.3f' % (popt))
plt.plot(xdata,ydata,marker='o')

popt2, pcov2 = curve_fit(power_law, xdata, ydata)
plt.plot(xdata, power_law(xdata, *popt2),label='Evo Alg model, power_law fit: a=%5.3f' % (popt))

ydata2=np.loadtxt("/Users/francescomambretti/Documents/Work/Synthetic_biodiversity_statistical_physics/experim_data/results_all/oligo1_proto06/R1_all_key_filter/min_Q_10/tot_valid.txt",unpack=True,usecols=(0,))
ydata2=ydata2*100
xdata2=np.arange(0,3*len(ydata2),3)
popt3, pcov3 = curve_fit(expo_decay, xdata2, ydata2)
plt.plot(xdata2, expo_decay(xdata2, *popt3),label='Exp data, expo fit: a=%5.3f' % (popt))
plt.plot(xdata2,ydata2,marker='o')

popt4, pcov4 = curve_fit(power_law, xdata2, ydata2)
plt.plot(xdata2, power_law(xdata2, *popt4),label='Exp data, power_law fit: a=%5.3f' % (popt))

plt.legend()
#plt.yscale('log')
plt.show()
'''
