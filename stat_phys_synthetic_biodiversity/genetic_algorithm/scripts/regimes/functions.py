## Useful functions for regime_omega.py

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from input_params import *

####################################### define methods #######################################

def generate_counts_list(bound_1,bound_2,MCO,skip,folder):
    a=0
    counts_list=np.zeros(len(cycles_IDs))
    histo_norm=np.zeros(len(cycles_IDs))

    for c in cycles_IDs:

        curr_folder=folder+"cycle_{}/".format(c)
            
        counts=np.empty(0)
        
        temp_MCO_list,abund_list=np.loadtxt(curr_folder+"unique_all_sorted.txt",skiprows=skip,usecols=(1,2,),unpack=True,ndmin=2)
        total=np.sum(abund_list)
        
        #now select only the data relative to the current MCO
        for i,j in zip (temp_MCO_list,abund_list):
            if i==MCO:
                counts=np.append(counts,j)
                
        counts[::-1].sort() #sort in decreasing abundance order
            
        if bound_2==-1:
            counts_list[a]=np.sum(counts)/total
           
        else:
            counts_list[a]=np.sum(counts[bound_1:bound_2])/total #skip some sequences, outside boundaries
    
        a+=1
        
    return counts_list
    
#####################################################################################################################

def plot_single_omega(MCO,n,folder,key):
    list1=generate_counts_list(0,-1,MCO,1,folder) #all - bound "-1" means: until the last one
    list2=generate_counts_list(n,-1,MCO,1,folder) #remove the first n most abundant ones
   
    f = open(folder+"dominant_{}.txt".format(MCO), "w")
    f.write("all all-dominant dominant \n")
    for el1,el2 in zip(list1,list2):
        f.write(str(el1)+" "+str(el2)+" "+str(el1-el2)+"\n")
    f.close()

    if key==True:
        plt.clf()

        plt.figure(figsize=(9,6))
        plt.title("$\omega=${}".format(MCO),fontsize=20)
        plt.plot(cycles_IDs,list1,marker='o',markersize=12,linestyle='none',color='gold',label='all')
        plt.plot(cycles_IDs,list2,marker='s',markersize=12,linestyle='none',color='teal',label='exclude dominant')
        plt.plot(cycles_IDs,list1-list2,marker='D',markersize=12,linestyle='none',color='navy',label='only dominant - $p-{}$'.format(MCO))

        plt.xlabel("Cycles",fontsize=14)
        plt.ylabel("Fraction of sequences with MCO={}".format(MCO),fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)

        plt.yscale('log')
        plt.savefig(folder+"excl_domin_omega_log={}.png".format(MCO),dpi=300)

        plt.yscale('linear')
        plt.savefig(folder+"excl_domin_omega_lin={}.png".format(MCO),dpi=300)
    
    #else --> do not plot
    
    return list1, list2, list1-list2

#####################################################################################################################

def plot_only_winner(MCO_list,folder,key):

    plt.figure(figsize=(12,9))
    
    n=len(MCO_list)
    color = iter(cm.turbo(np.linspace(0, 0.95, n)))
        
    for MCO in MCO_list:
        c = next(color)
        all_list,nodom_list,onlydom_list=np.loadtxt(folder+"dominant_{}.txt".format(MCO),unpack=True,
                                                    skiprows=1,usecols=(0,1,2))
        if key==True:
            plt.plot(cycles_IDs,onlydom_list,marker='s',markersize=12,linestyle='dashed',
                 color=c,label='$p-{}$'.format(MCO))
    if key==True:
        plt.xlabel("Cycles",fontsize=15)
        plt.ylabel("Sequences with a given MCO",fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.legend(fontsize=13,ncol=2,loc='lower left')
        plt.yscale('log')
        
        plt.savefig(folder+"dominant_log.png",dpi=200)
        
        plt.legend(fontsize=13,ncol=2,loc='upper left')
        plt.yscale('linear')
        plt.savefig(folder+"dominant_lin.png",dpi=200)
    
    #else --> do not plot!
        
    return

########################################################################################################################

def plot_regimes(folder,start,end,mycolor,regime,key):
    
    colorlist1=('black','brown','orangered')
    colorlist2=('purple','blue','green','mediumseagreen','turquoise')
    colorlist3=('orange','peru','red','firebrick','gray','navy')

    omega_list=np.zeros(max_cycle)
    omega_list2=np.zeros(max_cycle)
    norm=0
    
    plt.clf() # set plot
    plt.figure(figsize=(9,6))
    plt.title("${} \leq \omega \leq {}$".format(start,end),fontsize=20)

    for MCO in range (start,end+1):
        try:
            all_list,nodom_list,onlydom_list=np.loadtxt(folder+"/dominant_{}.txt".format(MCO),unpack=True,
                                                    skiprows=1,usecols=(0,1,2))
        except:
            all_list=np.zeros(max_cycle)
            nodom_list=np.zeros(max_cycle)
            onlydom_list=np.zeros(max_cycle)
            
        #plot single curves
        if(regime=="small"):
            plt.plot(cycles_IDs,all_list,linestyle='dashed',color=colorlist1[int(MCO-start)],
                 marker='s',label='$\omega$={}'.format(int(MCO)),alpha=0.9)

        elif(regime=="medium"):
            plt.plot(cycles_IDs,all_list,linestyle='dashed',color=colorlist2[int(MCO-start)],
                 marker='s',label='$\omega$={}'.format(int(MCO)),alpha=0.9)

        elif(regime=="large"):
            plt.plot(cycles_IDs,all_list,linestyle='dashed',color=colorlist3[int(MCO-start)],
                 marker='s',label='$\omega$={}'.format(int(MCO)),alpha=0.9)

        omega_list+=all_list
        omega_list2+=all_list*all_list
                
    norm=end-start+1
    omega_list/=norm
    omega_list2/=norm
    
    yerrors=np.sqrt(omega_list2-omega_list*omega_list)
    
    if key==True:
        
        plt.errorbar(cycles_IDs,omega_list,yerr=yerrors,capsize=5,marker='o',markersize=12,mec='black',ecolor='black',
                     linestyle='none',color=mycolor,label='mean')
        plt.fill_between(cycles_IDs, omega_list-yerrors, omega_list+yerrors,color=mycolor,alpha=0.08)
        
        plt.xlabel("Cycles",fontsize=14)
        plt.ylabel("Fraction of the total sequences",fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        #plt.show()
        plt.savefig(folder+"/ave_omega_behaviour_{}.png".format(regime),dpi=300)

    #else --> do not plot!
    
    return omega_list, omega_list2

#############################################  accu ##############################################

def accu(array,ave_array, ave2_array): #should work as well for scalars and arrays
    ave_array+=array
    ave2_array+=array*array
    
    return ave_array, ave2_array
    
#############################################  norm_and_err ##############################################

def norm_and_err(ave_array,ave2_array,norm_const):

    ave_array=ave_array/norm_const
    ave2_array=ave2_array/norm_const

    return ave_array, np.sqrt((ave2_array-ave_array*ave_array)/(norm_const))
