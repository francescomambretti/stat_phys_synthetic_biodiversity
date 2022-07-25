# set input options
# 25/07/2022 version

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

n_obs=6 #observables to be computed (and plotted): ShEn(p(w)), ShEn(RSA), <w>, # unique, top-n fraction of the total pop., zip ratio

i_shen_omega=0
i_shen_RSA=1
i_omega=2
i_unique=3
i_top_n=4
i_zip=5
# N.B. p(w) and cross-MCO are treated separatly, since they are not scalars

#################################### variables set inside if conditions ####################################

fastq_file_list=()
fastq_file_list_plot=()
cycles_IDs=np.arange(0,max_cycle,1)
mother_folder="../fitness_F_1E5_1E3/"
#mother_folder="../fitness_F_1E5_1E3/"
results_folder=mother_folder+"results"

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
                    
