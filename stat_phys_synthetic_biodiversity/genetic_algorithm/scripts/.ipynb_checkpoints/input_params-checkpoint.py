# set input options
# 11/07/2022 version

import numpy as np
import os

#################################### plots options, to be set ####################################
### parameters for plots
create_png_abund=False # create png files for abundance of strands inside MCO folders

cycles_plot_h_omega=(1,10,20,30,40,50)

#################################### analysis parameters ####################################

l=20 #resource length
n=10 #top-n strands
max_cycle=60

Np=1E5 #predators
Nr=1E3 #resources

max_seed=4 #number of independent simulations

freq_cross_MCO=5 # compute cross-MCO matrix among dominants every freq_cross_MCO cycles

#################################### variables set inside if conditions ####################################

fastq_file_list=()
fastq_file_list_plot=()
cycles_IDs=np.arange(0,max_cycle,1)
mother_folder="../fitness_F_short/"
results_folder=mother_folder+"results"

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
                    
