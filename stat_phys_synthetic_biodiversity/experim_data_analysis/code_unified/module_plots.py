# Written by Francesco Mambretti, 01/02/2022
# custom module for plotting:
# i) RSA histogram - Relative Species Abundance &
# ibis) Shannon entropy associated to RSA histogram, and also Simpson concentration and Gini-Simpson index
# ii) 2D/3D scatterplots of MCO and 2nd MCO
# iii) distribution of strand lengths
# iv) top-n trends
# v) abundance histogram
#
# 09/02/2022 version

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from input_params import *

############################################  plot_RSA ############################################

def plot_RSA_and_ShannonEntropy(input_file, upper_lim, output_file, output_ShEn_file, output_png, mylabel, colorname):

    df=pd.read_csv(input_file,sep='\t',dtype=np.uint8, usecols=('MCO',))
    data=df['MCO'].to_numpy()

    weights = np.ones_like(data)/len(data)

    histo,bins=np.histogram(data, bins=np.arange(0,upper_lim+2,1))
    np.savetxt(output_file, histo, fmt='%d')
    norm=np.sum(histo)
    
    #compute Shannon Entropy, Simpson concentration
    entropy=0
    concentration=0
    for bin_height in histo: #bin_height is the probability density associated to each bin
        if(bin_height!=0):
            prob=bin_height/norm
            entropy-=prob*np.log2(prob)
            concentration+=1./(prob*prob)
    #and update the relative file
    ShEn_file = open(output_ShEn_file,"a") #add the value for each experimental cycle
    string_to_write="{:.4f} {:.3e}".format(entropy,concentration)
    ShEn_file.write(string_to_write+'\n')
    ShEn_file.close()
        
    #plot
    plt.clf()
    
    fig, ax = plt.subplots(figsize = (9,6))
    axins = ax.inset_axes([0.68,0.68,0.3,0.3] )#Left, Bottom, Width, Height
    ax.hist(x=data, color=colorname,alpha=0.9, edgecolor='black',rwidth=0.75,align='left',bins=np.arange(0,upper_lim+1,1),weights=weights)# Main plot
    axins.hist(x=data, color=colorname,alpha=0.9, edgecolor='black',rwidth=0.75,align='left',bins=np.arange(0,upper_lim+1,1),weights=weights)# Inset plot
    ax.set_xlabel("MCO - species index",fontsize=19)#X axis title
    ax.set_ylabel("relative abundance",fontsize=19)#Y axis title
    ax.set_xticks(np.arange(0,21,2))
    ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
    
    axins.set_xticks(np.arange(0,21,5))
    axins.set_ylabel("log-scale",fontsize=14)
    axins.set_yticks((1e-06,1e-04,1e-02,1e-01,2e-01))
    ax.set_title(mylabel,fontsize=16)
    
    axins.set_yscale('log')
    ax.set_xlim(0,22)

    plt.savefig(output_png,dpi=300)
    
    plt.close(fig)
    
    return

############################################  scatter_MCO ############################################

def scatter_MCO (file, maximum, file_out_2D, word,specific_results,fastq_file):

    df=pd.read_csv(file,sep='\t',dtype=np.uint8)
    MCO=df['MCO'].to_numpy()
    MCO_2nd=df['MCO_2nd'].to_numpy()
    
    ndata=len(MCO)

    color_2D=np.zeros((maximum,maximum)) #set to 0 2D histogram entries

    for i in range(0,ndata): #count simultaneous occurrencies of MCO, MCO_2nd
        color_2D[int(MCO[i])][int(MCO_2nd[i])]+=1

    file_out = open(file_out_2D,"w")

    for i in range(0,ndata): #can be optimized, loop on much less data
        line="{} {} {:.4f} \n".format(MCO[i],MCO_2nd[i],color_2D[int(MCO[i])][int(MCO_2nd[i])]/ndata)
        file_out.write(line)

    file_out.close()

    fig = plt.figure(figsize = (9,6))

    color_array_2D=np.loadtxt(file_out_2D,usecols=(2,),unpack=True)
    
    fig=plt.figure(figsize = (9,6))
    plt.scatter(MCO, MCO_2nd, alpha=0.8, c=color_array_2D, edgecolors='none',cmap='coolwarm')
    plt.xlabel("MCO",fontsize=19)
    plt.ylabel("2nd MCO",fontsize=19)

    plt.colorbar(extend='both')
    
    plt.xticks(np.arange(0,21,2),fontsize=18)
    plt.yticks(np.arange(0,15,2),fontsize=18)
    plt.xlim(0,20)
    plt.ylim(0,14)
    
    plt.title("2nd MCO vs MCO, "+str(fastq_file))

    plt.savefig(specific_results+"/"+word+"_2D_scatter.png",dpi=300)

    plt.clf()
    plt.close(fig)

    return


############################################  dist_lengths ############################################

def dist_lengths(specific_results,fastq_file): # build distribution of oligomers lengths

    #dist_lengths is already an histogram
    dist_lengths=np.loadtxt(specific_results+"/dist_lengths.txt", unpack=True)
    plt.clf()
    
    fig, ax = plt.subplots(figsize = (9,6))
    ax.bar(x=np.arange(0,len(dist_lengths),1),height=dist_lengths, color='orangered',alpha=0.7, edgecolor='black')# Main plot
    ax.set_xlabel("strand length",fontsize=19)#X axis title
    ax.set_ylabel("counts",fontsize=19)#Y axis title
    ax.set_xticks(np.arange(0,len(dist_lengths)+1,5))
    ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
    
    plt.yscale('log')
    plt.xlim(lower_bound-1,L+1)
    
    plt.title(str(fastq_file),fontsize=19)
    
    plt.savefig(specific_results+"/dist_lengths.png",dpi=300)
    
    plt.close(fig)
    
    return

############################################ abund_histo ############################################

def abund_histo(results_folder, fastq_file, param, create_png_abund): # build histogram of strand abundances - useful to monitor variations of number of species present

    if (param!=-1):
        try:
            df=pd.read_csv(results_folder+"/"+fastq_file+"/MCO_{}/unique_MCO_sorted.txt".format(param),sep='\t',dtype=int, usecols=('abundance',))
            abundances=df['abundance'].to_numpy()
        except:
            pass
            
    else:
        df=pd.read_csv(results_folder+"/"+fastq_file+"/unique_all_sorted.txt",sep='\t',dtype=int, usecols=('abundance',))
        abundances=df['abundance'].to_numpy()

    try:
        weights = np.ones_like(abundances)/len(abundances)
        upper_lim = max(abundances)

        histo,bins=np.histogram(abundances, bins=np.arange(0,upper_lim+2,1))

        if (param!=-1):
            np.savetxt(results_folder+"/"+fastq_file+"/MCO_{}/abund_histo.txt".format(param), histo, fmt='%d')
        else:
            np.savetxt(results_folder+"/"+fastq_file+"/abund_histo.txt", histo, fmt='%d')
        
        #plt.clf()
        
        fig, ax = plt.subplots(figsize = (9,6))
        
        ax.hist(x=abundances, color="firebrick",alpha=1.0,rwidth=3,align='left',bins=np.arange(0,upper_lim+1,1),weights=weights)# Main plot
        ax.plot(np.arange(0,upper_lim+1,1),histo/np.sum(histo),marker='o',color='grey',markersize=5,linestyle='none')
       
        ax.set_title(str(fastq_file),fontsize=17)
        ax.set_xlabel("number of copies",fontsize=19)#X axis title
        ax.set_ylabel("counts",fontsize=19)#Y axis title
        
        ax.tick_params(axis = 'both', which = 'major', labelsize = 15)

        ax.set_yscale('log')
        
        if (create_png_abund==True and param!=-1):
            fig.savefig(results_folder+"/"+fastq_file+"/MCO_{}/abund_histo.png".format(param),dpi=300)
        else:
            fig.savefig(results_folder+"/"+fastq_file+"/abund_histo.png",dpi=300)
        
        plt.close(fig)

    except: #we have a 0-length array
        pass
        
    return

############################################ top_n_fraction ############################################

#determine which fraction of the overall number of strands is constituted by the top-n strands

def top_n_fraction(n,results_folder,fastq_file_list):
    length=len(fastq_file_list)
    fractions=np.zeros(0)
    
    for c in range (0,length): #loop over cycles
        fraction=0
        for i in range (0,n): #accumulate the n most abundant fractions
            fraction+=np.loadtxt(results_folder+"/"+fastq_file_list[c]+"/top_strand_{}_cycle_{}.txt".format(i,c*3),skiprows=c,max_rows=1,usecols=(1,)) # column 0 of the file contains cycle id
        
        fractions=np.append(fractions,fraction)
    np.savetxt(results_folder+"/top_{}_percentage_over_cycles.txt".format(n),fractions,fmt="%.4f")
        
    #do the plot
    plt.plot(np.arange(0,length*3,3),fractions,marker='s',color='purple')
    plt.xlabel("Cycle",fontsize=14)
    plt.ylabel("fraction of top-{} over total".format(n))
    plt.savefig(results_folder+"/top_{}_percentage_over_cycles.png".format(n),dpi=300)
    plt.close()
        
    return

############################################ top_n_plot ############################################

#different from previous versions, do a separate top_n_plot for each cycle, i.e. track the evolution in the past and in the following cycles for the top-n of each cycle.

def top_n_plot(n,results_folder, param, fastq_file_list):

    length=len(fastq_file_list)
    cycles=np.arange(0,length*3,3)
    
    colors=("black","blue","cornflowerblue","gray","teal","lime","seagreen","darkorchid","crimson","maroon")
    
    for c in range (0,length): #loop over cycles - select for each of them the top-n strands
        plt.clf()
        full_folder=results_folder+"/"+fastq_file_list[c]
        
        for i in range(0,n):
            if (param!=-1):
                try:
                    counts=np.loadtxt("/MCO_{}/top_strand_{}_cycle_{}.txt".format(param,i,c*3),unpack=True,usecols=(1,))
                except:
                    pass
            else:
                counts=np.loadtxt(full_folder+"/top_strand_{}_cycle_{}.txt".format(i,c*3),unpack=True,usecols=(1,))

            try:
                plt.plot(cycles,counts,color=colors[i],label=str(i+1),linestyle='-',marker='o')
            
            except:
                pass
        
        plt.title("Top-{} strands at cycle {}".format(n,c*3))
      
        plt.xlabel("Cycle",fontsize=17)
        
        plt.xticks(cycles)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
               
        plt.legend()
        
        if (param!=-1):
            plt.savefig(full_folder+"/MCO_{}/top_{}.png".format(param,n),dpi=300)
        else:
            plt.savefig(full_folder+"/top_{}.png".format(n),dpi=300)

        plt.close()

    return


############################################ plot_ShEn ############################################
def plot_ShEn(output_ShEn_file,word):

    entropy_array,simpson_array=np.loadtxt(output_ShEn_file,unpack=True)
        
    plt.clf()
    plt.plot(cycles_IDs,entropy_array,marker='o',color='navy',label="Shannon")
    plt.plot(cycles_IDs,simpson_array/np.mean(simpson_array),marker='s',color='brown',label="Simpson")
    plt.plot(cycles_IDs,1-simpson_array/np.mean(simpson_array),marker='s',color='gray',label="Gini-Simpson")
    plt.xlabel("Cycles",fontsize=14)
    plt.ylabel("Diversity index",fontsize=14)
    plt.legend()
    plt.savefig(results_folder+"/Shannon_Entropy_"+word+".png",dpi=300)
    plt.close()

    return
