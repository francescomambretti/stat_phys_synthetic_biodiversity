# Written by Francesco Mambretti, 11/07/2022
# custom module for plotting:
# i) omega histogram - Relative Species Abundance &
# ibis) Shannon entropy associated to omega histogram; Simpson concentration and Gini-Simpson index currenly commented
# ii) 2D histograms of MCO and 2nd MCO
# iii) distribution of strand lengths
# iv) top-n trends
# v) abundance histogram &
# vbis) Shannon Entropy of RSA, species is a unique sequence
# vi) # of unique sequences as a function of cycles
#
# 25/07/2022 version

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from input_params import *

############################################  plot_histo_omega ############################################

def plot_histo_omega (t,c): #do it for a given cycle t, choosing color c

    colorlist=('firebrick','orangered','gold','lime','green','steelblue','gray','red','pink','brown','magenta')

    ave_histo=np.zeros(l+1)
    ave_histo2=np.zeros(l+1)
    for s in range(1,max_seed+1,1):
        histo=np.loadtxt(mother_folder+"seed{}/cycle_{}/histogram.dat".format(s,t),unpack=True)
        histo/=np.sum(histo)
        ave_histo+=histo
        ave_histo2+=histo*histo

    #to plot histograms
    ave_histo/=max_seed
    ave_histo2/=max_seed
    error=np.sqrt((ave_histo2-ave_histo*ave_histo)/max_seed)
    
    plt.clf()

    plt.bar(np.arange(0,l+1,1),ave_histo,label='cycle {}'.format(t),alpha=1.0,align='center',color=colorlist[c],yerr=error,capsize=4,edgecolor='black')

    plt.xlabel("$\omega$",fontsize=12)
    plt.xticks(np.arange(0,l+2,4),fontsize=13)
    plt.yticks(fontsize=13)
    plt.legend(fontsize=13)
    #plt.title(title_histo) #input parameter
    plt.savefig(results_folder+"/ave_h_omega_{}_{}_{}.png".format(Np,Nr,t),dpi=300)

    plt.yscale('log')
    plt.savefig(results_folder+"/ave_h_omega_{}_{}_{}_log.png".format(Np,Nr,t),dpi=300)
    
    plt.close()

    return

############################################ plot_histo_ShEn_RSA ############################################
'''
def plot_histo_ShEn_RSA(cycle_folder, create_png_abund, output_ShEn_file): # build histogram of strand abundances - useful to monitor variations of number of species present

    df=pd.read_csv(cycle_folder+"/unique_all_sorted.txt",sep='\t',dtype=int, usecols=('abundance',))
    abundances=df['abundance'].to_numpy()

    #compute entropy of the RSA of the whole population only
    total=np.sum(abundances) #total sequences
    # compute Shannon Entropy of the RSA - species:= sequence identity
    entropy=0
    for counts in abundances: #bin_height is the probability density associated to each bin
        prob=counts/total
        entropy-=prob*np.log2(prob)
    #and update the relative file
    ShEn_file = open(output_ShEn_file,"a") #add the value for each experimental cycle
    string_to_write="{:.4f}".format(entropy)
    ShEn_file.write(string_to_write+'\n')
    ShEn_file.close()

    try:
        # build histogram counting the number of sequences with 1,2,3 copies...
        weights = np.ones_like(abundances)/len(abundances)
        upper_lim = max(abundances)
        histo,bins=np.histogram(abundances, bins=np.arange(0,upper_lim+2,1))
        norm=np.sum(histo)

        np.savetxt(cycle_folder+"/abund_histo.txt", histo, fmt='%d')
        
        fig, ax = plt.subplots(figsize = (9,6))
        
        ax.hist(x=abundances, color="lightcoral",alpha=1.0,rwidth=3,align='left',bins=np.arange(0,upper_lim+1,1),weights=weights)# Main plot
        ax.plot(np.arange(0,upper_lim+1,1),histo/np.sum(histo),marker='o',color='grey',markersize=5,linestyle='none')
       
        ax.set_title(str(fastq_file),fontsize=17)
        ax.set_xlabel("number of copies",fontsize=19)#X axis title
        ax.set_ylabel("counts",fontsize=19)#Y axis title
        
        ax.tick_params(axis = 'both', which = 'major', labelsize = 15)
        ax.set_yscale('log')
        
        if (create_png_abund==True):
            fig.savefig(results_folder+"/"+fastq_file+"/abund_histo.png",dpi=300)
        
        plt.close(fig)
        
    except: #i.e., abundances array is empty
        pass
    return'''

############################################ top_n_fraction ############################################

#determine which fraction of the overall number of strands is constituted by the top-n strands

def top_n_fraction(n,folder,key):

    plt.xlabel("Cycle",fontsize=14)
    plt.ylabel("fraction of top-{} over total".format(n))
    
    if key==False:
        fractions=np.loadtxt(folder+"/top_{}_percentage_over_cycles.txt".format(n),unpack=True)
        plt.plot(cycles_IDs,fractions,marker='s',color='purple')
    else:
        fractions,errors=np.loadtxt(folder+"/mean_top_{}_percentage_over_cycles.txt".format(n),unpack=True,usecols=(0,1))
        #do the plot
        plt.errorbar(cycles_IDs,fractions,yerr=errors,marker='s',color='purple',capsize=5)
    plt.savefig(folder+"/top_{}_percentage_over_cycles.png".format(n),dpi=300)
        
    plt.close()
        
    return

############################################ top_n_plot ############################################

#different from previous versions, do a separate top_n_plot for each cycle, i.e. track the evolution in the past and in the following cycles for the top-n of each cycle.

def top_n_plot(n,cycle_folder):

    a=0
    
    colors=("black","blue","cornflowerblue","gray","teal","lime","seagreen","darkorchid","crimson","maroon")
    
    for c in cycles_IDs: #loop over cycles - select for each of them the top-n strands
        plt.clf()
        for i in range(0,n): # for the n-th most abundant sequence at last cycle
            try:
                counts=np.loadtxt(cycle_folder+"/cycle_{}/top_{}_seq_history.txt".format(c,i),unpack=True,usecols=(1,))
                plt.plot(cycles_IDs,counts,color=colors[i],label=str(i+1),linestyle='-',marker='o')
        
            except: #in case the file does not exist
                pass
        
        plt.title("Top-{} sequences at cycle {}".format(n,c))
      
        plt.xlabel("Cycle",fontsize=13)
        plt.ylabel("Fraction of the total population", fontsize=13)
        plt.xticks(np.arange(0,max_cycle,10))
        plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
               
        plt.legend()
        plt.yscale('log')
        
        plt.tight_layout()
        plt.savefig(cycle_folder+"/top_{}.png".format(n),dpi=300)

        plt.close()
        
        a+=1

    return


############################################ plot_ShEn ############################################
def plot_ShEn(output_ShEn_file,ShEn_png,my_color):

    try:
        entropy_array=np.loadtxt(output_ShEn_file,unpack=True,usecols=(0,))
            
        plt.clf()
        plt.plot(cycles_IDs,entropy_array,marker='o',color=my_color,label="Shannon")
        plt.xlabel("Cycles",fontsize=14)
        plt.ylabel("Shannon entropy",fontsize=14)
        plt.legend()
        plt.savefig(ShEn_png,dpi=300)
        plt.close()
        
    except:
        pass

    return


############################################ plot_obs_ave_parall ############################################

def plot_obs_ave_parall(file,pngfile,my_label,my_color): #OK for all the observables obtained as time series, averaged over parallel runs

    data,errors=np.loadtxt(file,unpack=True,usecols=(0,1))
    
    plt.clf()
    plt.errorbar(cycles_IDs,data,errors,marker='o',capsize=5,linestyle='dotted',color=my_color,label=my_label)
    plt.xlabel("Cycles",fontsize=14)
    plt.ylabel(my_label,fontsize=14)
    plt.legend()
    plt.savefig(pngfile,dpi=300)
    plt.close()
    
    return

############################################ plot_cross_MCOs ############################################

def plot_cross_MCOs(file,pngfile,cycle,key):

    if key==True:
        labels=np.loadtxt(file,unpack=True,max_rows=1,dtype=int)
        cross_MCOs=np.loadtxt(file,unpack=True,skiprows=1)
    else:
        cross_MCOs=np.loadtxt(file,unpack=True,skiprows=1,max_rows=n)
    
    plt.clf()
    cax = plt.imshow(cross_MCOs,cmap='viridis')
    plt.title("Cycle {}".format(cycle),fontsize=14)
    plt.xlabel("$\omega$",fontsize=14)
    plt.ylabel("$\omega$",fontsize=14)
    cbar = plt.colorbar(cax, ticks=[0,5,10,15])
    plt.clim(0,15)
    if key==True:
        plt.xticks(np.arange(0,n,1),label=labels)
        plt.yticks(np.arange(0,n,1),label=labels)
    plt.savefig(pngfile,dpi=300)
    plt.close()
    
    return
