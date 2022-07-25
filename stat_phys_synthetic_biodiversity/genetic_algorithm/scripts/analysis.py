# Python code to analyze outcomes of  Evolutive (Genetic) Algorithm simulations - based on read_fastq.py and other codes for the analysis of experimental SEDES data
#
# Written by Francesco Mambretti, 11/07/2022
#
# main script. The code processes all the simulations and analyses them, computing:
# 0) the mean histogram p(MCO)
# 1) the associated Shannon Entropy
# 2) the number of unique sequences
# 3) the Shannon Entropy associated to the Relative Species Abundance, where species=sequence
# 4) the average omega of the population
# 5) the zip ratio between the original and the zipped file size for the whole list
# 6) the evolution of most abundant strands at each cycle
# 7) the fraction of population covered by top-n individuals
#
# 25/07/2022 version

############################################## import ##############################################
import numpy as np
import sys
import module_functions as func # import module_functions as func
from input_params import *
import time #Monitor performance

############################################## declare variables ##############################################
start=time.perf_counter()

list_files=["ShEn_omega","ShEn_RSA","omega_time","unique_sequences","top_{}_percentage_over_cycles".format(n),"zip ratio"]
list_mean_files=["mean_ShEn_omega","mean_ShEn_RSA","mean_omega_time","mean_unique_sequences","mean_top_{}_percentage_over_cycles".format(n),"mean_zip_ratio"]

mean=np.zeros((n_obs,max_cycle)) #empty containers for scalar variables
mean2=np.zeros((n_obs,max_cycle))
err=np.zeros((n_obs,max_cycle))

mean_histo_omega=np.zeros((max_cycle,l+1))
mean2_histo_omega=np.zeros((max_cycle,l+1))
err_histo_omega=np.zeros((max_cycle,l+1))

mean_cross_MCOs = np.zeros((max_cycle,n+1,n))
mean2_cross_MCOs = np.zeros((max_cycle,n+1,n))
err_cross_MCOs = np.zeros((max_cycle,n+1,n))

############################################## loop through cycles ##############################################

for s in range (1,max_seed+1,1):
    print("analyzing simulation {}".format(s))
    seed_folder=mother_folder+"/seed{}".format(s)

    for t in cycles_IDs: #loop through cycles
        cycle_folder=mother_folder+"/seed{}/cycle_{}/".format(s,t)

        if t==0: #clean
            os.system("rm "+seed_folder+"/"+list_files[i_shen_omega]+".txt")
            os.system("rm "+seed_folder+"/"+list_files[i_shen_RSA]+".txt")
    
        histo,entropy=func.compute_ShEn_omega(cycle_folder+"/histogram.dat",seed_folder+"/"+list_files[i_shen_omega]+".txt")
        
        #accumulate - DO NOT CHANGE ORDER OF THIS FUNCTIONS, THE CODE WILL NOT WORK
        mean_histo_omega[t],mean2_histo_omega[t]=func.accu(histo,mean_histo_omega[t],mean2_histo_omega[t])
        
        mean[i_shen_omega][t],mean2[i_shen_omega][t]=func.accu(entropy,mean[i_shen_omega][t],mean2[i_shen_omega][t])
        
        mean[i_unique][t],mean2[i_unique][t] = func.accu(func.detect_unique_seqs(cycle_folder),mean[i_unique][t],mean2[i_unique][t])

        func.sort_unique_strands(cycle_folder) #sort by abundance
        
        mean[i_omega][t],mean2[i_omega][t] = func.accu(func.ave_omega_pop(cycle_folder),mean[i_omega][t],mean2[i_omega][t])
        
        mean[i_shen_RSA][t],mean2[i_shen_RSA][t] =func.accu(func.compute_ShEn_RSA(cycle_folder,seed_folder+"/"+list_files[i_shen_RSA]+".txt"),mean[i_shen_RSA][t],mean2[i_shen_RSA][t])
        
        mean[i_zip][t],mean2[i_zip][t] =func.accu(func.zip_ratio(cycle_folder+"/list_seqs_all.dat",cycle_folder+"/list_seqs_all.zip"),mean[i_zip][t],mean2[i_zip][t])

    func.top_n_full_evo(n,seed_folder,s) #this is specific to each seed
    
    for t in range(0,max_cycle,freq_cross_MCO): #necessarily after top_n_full_evo
        mean_cross_MCOs[t],mean2_cross_MCOs[t]=func.accu(np.loadtxt(seed_folder+"/cycle_{}/cross_MCOs.txt".format(t)),mean_cross_MCOs[t],mean2_cross_MCOs[t])
    
    #now, use data produced by top_n_full_evo to accumulate estimation for mean[i_top_n]
    mean[i_top_n],mean2[i_top_n] = func.accu(func.compute_top_n_fraction(seed_folder),mean[i_top_n],mean2[i_top_n])
    
################ compute observables averages over parallel simulations ################
for t in cycles_IDs: #loop through cycles
    mean_histo_omega[t],err_histo_omega[t]=func.norm_and_err(mean_histo_omega[t],mean2_histo_omega[t],max_seed)
    if not os.path.exists(results_folder+"/cycle_{}".format(t)): # For each cycle, create the corresponding folder
        os.makedirs(results_folder+"/cycle_{}".format(t))
    np.savetxt(results_folder+"/cycle_{}".format(t)+"/mean_histo_omega.txt",np.column_stack((mean_histo_omega[t],err_histo_omega[t])),fmt='%.2f')
    
    if (t%freq_cross_MCO==0):
        mean_cross_MCOs[t],err_cross_MCOs[t]=func.norm_and_err(mean_cross_MCOs[t],mean2_cross_MCOs[t],max_seed)
        np.savetxt(results_folder+"/cycle_{}".format(t)+"/mean_cross_MCOs.txt",np.vstack((mean_cross_MCOs[t],err_cross_MCOs[t])),fmt='%.2f')

for o in range (0,n_obs):
    mean[o],err[o]=func.norm_and_err(mean[o],mean2[o],max_seed) #normalize data coming from all seeds and compute errors
    
    np.savetxt(results_folder+"/"+list_mean_files[o]+".txt",np.column_stack((mean[o],err[o])),fmt='%.2f')

############################################## plots ##############################################

cp1=time.perf_counter()
print (str(cp1-start)+" seconds spent since beginning")
#do the plots
os.system("python main_plot.py")
cp2=time.perf_counter()
print (str(cp2-cp1)+" seconds spent plotting")

