# set input options
# 25/07/2022 version

import numpy as np
import os
import matplotlib.pyplot as plt

#################################### plots options, to be set ####################################
### parameters for plots
create_png_abund=False # create png files for abundance of strands inside MCO folders

cycles_plot_h_omega=(1,5,10,15,20,25,30,35,40,45,49)

#################################### analysis parameters ####################################

l=20 #resource length
n=10 #top-n strands
max_cycle=25

Np=1E6 #predators
Nr=1E4 #resources
N_p=str(Np)
N_r=str(Nr)

asymp_omega=10
asymp=True

max_seed=5 #number of independent simulations

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
if asymp==True:
    mother_folder="../size_effects/FS_11_alpha_0.1/Np_1000000_Nr_10000/"
#    mother_folder="../FS11_alpha_0.02/"
    print(mother_folder)
else:
    mother_folder="../fitness_F_5E6_5E4/"
    print(mother_folder)
results_folder=mother_folder+"results"

if not os.path.exists(results_folder):
    os.makedirs(results_folder)

#if asymp==True:
#    title_histo=#"$f = 10 \\times \Theta({}-\omega) + \\dfrac{\omega}{\Sigma \omega} \\times \Theta(\omega-{}) $ - $N_{p}={}, N_{r}={}$".format(asymp_omega,Np,Nr)
#else:
#    title_histo="$\\dfrac{\omega}{\Sigma \omega}$ - $N_{p}={}, N_{r}={}$".format(Np,Nr)
