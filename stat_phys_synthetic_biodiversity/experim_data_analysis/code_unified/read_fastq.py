# Python code to analyze experimental FASTQ files coming from SEDES experiments
#
# Written by Francesco Mambretti, 01/02/2022
#
# Unified version of previous codes, no need to change source code for different datasets (i.e. different primers, target strand, output options, R1/R2, fw/rev). Now all is included in input_params.py
#
# The code processes all the files included in fastq_file_list and analyses them, conditioned to the quality score of the reads
# 0) extracts the strand sequence from FASTQ files
# 1) computes the distribution of the strands length
# 2) computes MCO and MCO_2nd for each strand
# 3) computes the RSA histogram, where a species is defined according to its max(MCO)=\omega & compute their Shannon entropy
# 4) create a separate folder for each MCO where to save MCO-specific results
# 5) detect abundance of unique strands
# 6) track the evolution of most abundant strands at each cycle
# 7) compute the fraction of population covered by top-10 individuals
#
# Sequences which are not recognized either as forward or as reverse, are not analyzed and counted as non valid ones.
#
# 09/02/2022 version

##### import #####
import numpy as np
import os
import sys
from Bio import SeqIO
import module_functions as func
from input_params import *
import time #Monitor performance

##### read FASTQ #####
start=time.perf_counter()

t=0 # keep track of current FASTQ file analyzed - incremented by 1 after each one

if (key2=="R1" or key2=="R2"):
    print(key1,key2,key3)
    tot_valid=np.zeros(len(fastq_file_list))
    unique_time_series=np.zeros(len(fastq_file_list))
    
    for fastq_file in fastq_file_list:
        tot_valid[t]=func.process_FASTQ(fastq_file,mother_folder=mother_folder,results_folder=results_folder,stop=subset_steps,use_stop=use_stop,t=t)
        #open FASTQ file, read, extract strand sequences, distinguish reverse and forward. Then compute MCO and MCO_2nd for each strand
        specific_results=results_folder+"/"+str(t)+"_"+fastq_file

        # now, we have MCO, MCO_2nd and MCO_3rd for all the strands
        # create a separate folder, one for each MCO, where to store MCO-specific results
        func.create_separate_MCO_folders(l,specific_results)
    
        # then, analyze all the strands and identify dominating sequences (sort by number of copies)
        # the following function also saves separated results for each MCO (so to have MCO-labelled families)
        # also plots the distribution of strand lengths
        unique_time_series[t]=func.detect_unique_strands(specific_results)
        t+=1
    
        #print percentages of unique strands as a function of time
        func.sort_unique_strands(specific_results)

    print(unique_time_series/tot_valid)
    
    func.top_n_full_evo (fastq_file_list,n=n,results_folder=results_folder,param=-1,tot_valid=tot_valid) #overall

elif (key2=="R1R2"):
    tot_valid=np.zeros(len(fastq_file_list_R1))
    unique_time_series=np.zeros(len(fastq_file_list_R1))
    for fastq_file_R1,fastq_file_R2 in zip(fastq_file_list_R1,fastq_file_list_R2):
        tot_valid[t]=func.process_FASTQ(fastq_file_R1,fastq_file_R2,mother_folder=mother_folder,results_folder=results_folder,stop=subset_steps,use_stop=use_stop,t=t)
        
        specific_results=results_folder+"/"+str(t)+"_"+fastq_file_R1+"R2"
        func.create_separate_MCO_folders(l,specific_results)
        unique_time_series[t]=func.detect_unique_strands(specific_results)
        t+=1
        func.sort_unique_strands(specific_results)
    
    #print percentages of unique strands as a function of time
    print(unique_time_series/tot_valid)

    func.top_n_full_evo (fastq_file_list_R1,fastq_file_list_R2,n=n,results_folder=results_folder,param=-1,tot_valid=tot_valid) #overall
    #for temp_mco in np.arange(0,l+1): #for each MCO - uncomment iff really needed
    #   func.top_n_full_evo (n,results_folder,fastq_file_list,-1,tot_valid)
    #   func.top_n_trend (10,results_folder,fastq_file_list,temp_mco,tot_valid)

cp1=time.perf_counter()
print (str(cp1-start)+" seconds spent since beginning")
#do the plots
os.system("python main_plot.py")
cp2=time.perf_counter()
print (str(cp2-cp1)+" seconds spent plotting")
