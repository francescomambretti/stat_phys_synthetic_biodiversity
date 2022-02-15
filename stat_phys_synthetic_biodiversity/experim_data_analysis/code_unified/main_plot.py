# Python code to plot experimental SEDES data
# Written by Francesco Mambretti, 01/02/2022
#
# Checkpoints for performance monitoring are currently commented
#
# 08/02/2022 version

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from mpl_toolkits.mplot3d import axes3d
import module_plots
from input_params import *
import time

maximum=l+1 #go from 0 to l
t=-1

for fastq_file in fastq_file_list_plot:   #loop through cycles
    start=time.perf_counter()
    t+=1
    specific_results=results_folder+"/"+fastq_file
    module_plots.dist_lengths(specific_results,fastq_file)
    
    #cp1=time.perf_counter()
    #print (str(cp1-start)+" seconds spent until dist_lengths")

    if(key3=="all"):
        for word,colorname in zip(keywords,colornames): #separate analysis for all,fw,rev

            #plot frequency of (MCO and MCO_2nd) pairs
            file=specific_results+"/"+word+"_MCO.txt"
            module_plots.scatter_MCO(file, maximum, specific_results+"/"+word+"_MCO_2D_color.txt", word, specific_results,fastq_file)
                    
            #cp2=time.perf_counter()
            #print (str(cp2-cp1)+" seconds spent after scatter")

            # plot RSA histogram - species defined according to its affinity to target strand
            # and associated Shannon entropy
            output_RSA_file=specific_results+"/histo_RSA_"+word+".txt"
            output_RSA_png=specific_results+"/histo_RSA_"+word+".png"
            output_ShEn_file=results_folder+"/Shannon_Entropy_"+word+".txt"
            if (t==0):
                os.system("rm "+output_ShEn_file)
                
            #cp3=time.perf_counter()
            #print (str(cp3-cp2)+" seconds spent after RSA")
            
            module_plots.plot_RSA_and_ShannonEntropy(file, maximum, output_RSA_file, output_ShEn_file, output_RSA_png, word+" strands - cycle "+fastq_file, colorname)
            
            # abundance of strands
            module_plots.abund_histo(results_folder,fastq_file,-1,create_png_abund) #overall
            for temp_mco in np.arange(0,maximum): #for each MCO - default option: only generate txt file, do not plot
                module_plots.abund_histo(results_folder,fastq_file,temp_mco,create_png_abund)
                
            #cp4=time.perf_counter()
            #print (str(cp4-cp3)+" seconds spent after abund histo")

    else:
        #plot frequency of (MCO and MCO_2nd) pairs
        file=specific_results+"/all_MCO.txt"
        module_plots.scatter_MCO(file, maximum, specific_results+"/all_MCO_2D_color.txt", "all", specific_results,fastq_file)
        
        #cp2=time.perf_counter()
        #print (str(cp2-cp1)+" seconds spent after scatter")

        # plot RSA histogram - species defined according to its affinity to target strand
        # and associated Shannon entropy
        output_RSA_file=specific_results+"/histo_RSA_all.txt"
        output_RSA_png=specific_results+"/histo_RSA_all.png"
        output_ShEn_file=results_folder+"/Shannon_Entropy_all.txt"
        if (t==0):
            os.system("rm "+output_ShEn_file)
            
        #cp3=time.perf_counter()
        #print (str(cp3-cp2)+" seconds spent after RSA")
        
        module_plots.plot_RSA_and_ShannonEntropy(file, maximum, output_RSA_file, output_ShEn_file, output_RSA_png, "all strands - cycle "+fastq_file, colornames[0])
        
        # abundance of strands
        module_plots.abund_histo(results_folder,fastq_file,-1,create_png_abund) #overall
        for temp_mco in np.arange(0,maximum): #for each MCO - default option: only generate txt file, do not plot
            module_plots.abund_histo(results_folder,fastq_file,temp_mco,create_png_abund)
                
        #cp4=time.perf_counter()
        #print (str(cp4-cp3)+" seconds spent after abund histo")
#after all the cycles have been analyzed...
# now plot the evolution trend of the dominant strands
module_plots.top_n_plot(n,results_folder,-1,fastq_file_list_plot) #overall
module_plots.top_n_fraction(n,results_folder,fastq_file_list_plot)

if(key3=="all"):
    for word,colorname in zip(keywords,colornames):
        output_ShEn_file=results_folder+"/Shannon_Entropy_"+word+".txt"
        module_plots.plot_ShEn(output_ShEn_file,word) #plot the Shannon entropy as a function of the cycle index

        #for temp_mco in np.arange(0,maximum): #for each MCO - commented by default
        #    module_plots.top_n_plot(n,results_folder,temp_mco)
                
    #cp5=time.perf_counter()
    #print (str(cp5-cp4)+" seconds spent after top")
else:
    output_ShEn_file=results_folder+"/Shannon_Entropy_all.txt"
    module_plots.plot_ShEn(output_ShEn_file,"all") #plot the Shannon entropy as a function of the cycle index
    #cp5=time.perf_counter()
    #print (str(cp5-cp4)+" seconds spent after top")
