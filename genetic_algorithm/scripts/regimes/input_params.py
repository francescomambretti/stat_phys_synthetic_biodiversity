# set input options
# 17/07/2022 version

import numpy as np
import os

l=20 #resource length
n=10 #top-n strands
max_cycle=60

cycles_IDs=np.arange(0,max_cycle,1)

Np=1E5 #predators
Nr=1E3 #resources

min_w=4
max_w=13
num_omegas=max_w-min_w

max_seed=4 #number of independent simulations

freq_cross_MCO=5 # compute cross-MCO matrix among dominants every freq_cross_MCO cycles

#################################### variables set inside if conditions ####################################

fastq_file_list=()
fastq_file_list_plot=()
cycles_IDs=np.arange(0,max_cycle,1)
#mother_folder="../../fitness_FS8_1E5_1E3/"
mother_folder="../../fitness_F_1E5_1E3/"
results_folder=mother_folder+"results"

if not os.path.exists(results_folder):
    os.makedirs(results_folder)

######################
small_b=4
small_e=5
medium_b=6
medium_e=7 #8
large_b=8 #9
large_e=10 #12
