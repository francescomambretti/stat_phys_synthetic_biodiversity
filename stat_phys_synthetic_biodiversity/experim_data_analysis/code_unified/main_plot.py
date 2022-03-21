# Python code to plot experimental SEDES data
# Written by Francesco Mambretti, 01/02/2022
#
# Checkpoints for performance monitoring are currently commented
#
# 21/03/2022 version

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from mpl_toolkits.mplot3d import axes3d
import module_plots
from input_params import *
#import time

maximum=l+1 #go from 0 to l
t=-1

for fastq_file in fastq_file_list_plot:   #loop through cycles
    #start=time.perf_counter()
    t+=1
    specific_results=results_folder+"/"+fastq_file
    module_plots.dist_lengths(specific_results,fastq_file)
    
    #MCO cross-matrix among top-n
    module_plots.plot_cross_MCOs(specific_results+"/cross_MCOs.txt",specific_results+"/cross_MCOs.png",t)
    
    #cp1=time.perf_counter()
    #print (str(cp1-start)+" seconds spent until dist_lengths")

    if(key3=="all"):
        for word,colorname in zip(keywords,colornames): #separate analysis for all,fw,rev

            #plot frequency of (MCO and MCO_2nd) pairs
            file=specific_results+"/"+word+"_MCO.txt"
            module_plots.MCO_MCO2nd(file, maximum, specific_results+"/"+word+"_MCO_2D_color.txt", word, specific_results,fastq_file)
                    
            #cp2=time.perf_counter()
            #print (str(cp2-cp1)+" seconds spent after scatter")

            # plot omega histogram - species defined according to its affinity to target strand
            # and associated Shannon entropy
            output_omega_file=specific_results+"/histo_omega_"+word+".txt"
            output_omega_png=specific_results+"/histo_omega_"+word+".png"
            ShEn_omega_file=results_folder+"/Shannon_Entropy_omega_"+word+".txt"
            ShEn_RSA_file=results_folder+"/Shannon_Entropy_RSA_"+word+".txt"
            if (t==0):
                os.system("rm "+ShEn_omega_file)
                os.system("rm "+ShEn_RSA_file)
                
            #cp3=time.perf_counter()
            #print (str(cp3-cp2)+" seconds spent after omega")
            
            module_plots.plot_omega_ShEn(file, maximum, output_omega_file, ShEn_omega_file, output_omega_png, word+" strands - cycle "+fastq_file, colorname)
            
            # abundance of strands
            module_plots.abund_histo_ShEn_RSA(results_folder,fastq_file,-1,create_png_abund,ShEn_RSA_file) #overall
            for temp_mco in np.arange(0,maximum): #for each MCO - default option: only generate txt file, do not plot
                ShEn_RSA_file=results_folder+"/"+fastq_file+"/MCO_{}".format(temp_mco)+"/Shannon_Entropy_RSA_"+word+".txt"
                module_plots.abund_histo_ShEn_RSA(results_folder,fastq_file,temp_mco,create_png_abund,ShEn_RSA_file)
                
            #cp4=time.perf_counter()
            #print (str(cp4-cp3)+" seconds spent after abund histo")

    else:
        #plot frequency of (MCO and MCO_2nd) pairs
        file=specific_results+"/all_MCO.txt"
        module_plots.MCO_MCO2nd(file, maximum, specific_results+"/all_MCO_2D_color.txt", "all", specific_results,fastq_file)
        
        #cp2=time.perf_counter()
        #print (str(cp2-cp1)+" seconds spent after 2D hist")

        # plot omega histogram - species defined according to its affinity to target strand
        # and associated Shannon entropy
        output_omega_file=specific_results+"/histo_omega_all.txt"
        output_omega_png=specific_results+"/histo_omega_all.png"
        ShEn_omega_file=results_folder+"/Shannon_Entropy_omega_all.txt"
        ShEn_RSA_file=results_folder+"/Shannon_Entropy_RSA_all.txt"
        if (t==0):
            os.system("rm "+ShEn_omega_file)
            os.system("rm "+ShEn_RSA_file)
            
        #cp3=time.perf_counter()
        #print (str(cp3-cp2)+" seconds spent after omega")
        
        module_plots.plot_omega_ShEn(file, maximum, output_omega_file, ShEn_omega_file, output_omega_png, "all strands - cycle "+fastq_file, colornames[0])
        
        # abundance of strands
        module_plots.abund_histo_ShEn_RSA(results_folder,fastq_file,-1,create_png_abund,ShEn_RSA_file) #overall
        for temp_mco in np.arange(0,maximum): #for each MCO - default option: only generate txt file, do not plot
            ShEn_RSA_file=results_folder+"/"+fastq_file+"/MCO_{}/".format(temp_mco)+"/Shannon_Entropy_RSA_all.txt"
            module_plots.abund_histo_ShEn_RSA(results_folder,fastq_file,temp_mco,create_png_abund,ShEn_RSA_file)
                
        #cp4=time.perf_counter()
        #print (str(cp4-cp3)+" seconds spent after abund histo")
        
#after all the cycles have been analyzed...
# now plot the evolution trend of the dominant strands
module_plots.top_n_plot(n,results_folder,-1,fastq_file_list_plot) #overall
module_plots.top_n_fraction(n,results_folder,fastq_file_list_plot)

if(key3=="all"):
    for word,colorname in zip(keywords,colornames):
        ShEn_omega_file=results_folder+"/Shannon_Entropy_omega_"+word+".txt"
        ShEn_omega_png=results_folder+"/Shannon_Entropy_omega_"+word+".png"
        module_plots.plot_ShEn(ShEn_omega_file,word,ShEn_omega_png) #plot the Shannon entropy of the omega as a function of the cycle index
        
        ShEn_RSA_file=results_folder+"/Shannon_Entropy_RSA_"+word+".txt"
        ShEn_RSA_png=results_folder+"/Shannon_Entropy_RSA_"+word+".png"
        module_plots.plot_ShEn(ShEn_RSA_file,word,ShEn_RSA_png) #plot the Shannon entropy of the species as a function of the cycle index
        for temp_mco in np.arange(0,maximum):
            ShEn_RSA_file=results_folder+"/MCO_{}".format(temp_mco)+"/Shannon_Entropy_RSA_"+word+".txt"
            ShEn_RSA_png=results_folder+"/MCO_{}".format(temp_mco)+"/Shannon_Entropy_RSA_"+word+".png"
            module_plots.plot_ShEn(ShEn_RSA_file,word,ShEn_RSA_png)

        #for temp_mco in np.arange(0,maximum): #for each MCO - commented by default
        #    module_plots.top_n_plot(n,results_folder,temp_mco)
                
    #cp5=time.perf_counter()
    #print (str(cp5-cp4)+" seconds spent after top")
else:
    ShEn_omega_file=results_folder+"/Shannon_Entropy_omega_all.txt"
    ShEn_omega_png=results_folder+"/Shannon_Entropy_omega_all.png"
    module_plots.plot_ShEn(ShEn_omega_file,word,ShEn_omega_png) #plot the Shannon entropy of the omega as a function of the cycle index
    ShEn_RSA_file=results_folder+"/Shannon_Entropy_RSA_all.txt"
    ShEn_RSA_png=results_folder+"/Shannon_Entropy_RSA_all.png"
    module_plots.plot_ShEn(ShEn_RSA_file,word,ShEn_RSA_png) #plot the Shannon entropy of the species as a function of the cycle index
    for temp_mco in np.arange(0,maximum):
        ShEn_RSA_file=results_folder+"/MCO_{}/Shannon_Entropy_RSA_all.txt".format(temp_mco)
        ShEn_RSA_png=results_folder+"/MCO_{}/Shannon_Entropy_RSA_all.png".format(temp_mco)
        module_plots.plot_ShEn(ShEn_RSA_file,word,ShEn_RSA_png)
    #cp5=time.perf_counter()
    #print (str(cp5-cp4)+" seconds spent after top")

# plot # of unique sequences as a function of time
module_plots.plot_unique_sequences(results_folder+"/tot_valid.txt",results_folder+"/tot_valid.png")
