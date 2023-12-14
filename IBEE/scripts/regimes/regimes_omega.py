## Written by Francesco Mambretti
## 15/07/2022
## To analyze Genetic Algorithm results, identifying three different regimes (in omega) for predators abundances in time
## 03/01/2023  version

import os
import functions as func
from input_params import *

#empty containers for scalar variables
#mean=np.zeros((3,small_e-small_b,max_cycle))
#mean2=np.zeros((3,medium_e-medium_b,max_cycle))
#err=np.zeros((3,large_e-large_b,max_cycle))

for s in range (1,max_seed+1,1):
    seed_folder=mother_folder+"/seed{}/".format(s)
    #clean
    os.system("rm "+seed_folder+"/*domin*".format(s))
    for MCO in range (min_w,max_w,1):
        all,nodom,onlydom = func.plot_single_omega(int(MCO),1,seed_folder,False) #delete only the first most abundant one
    
    # compare trends of winners of the various omega families
    func.plot_only_winner(range (min_w,max_w,1),seed_folder,False)
    
# generate files with averages
for MCO in range (min_w,max_w,1):
    ave_vec=np.zeros((3,max_cycle))
    ave_vec2=np.zeros((3,max_cycle))
    for s in range (1,max_seed+1,1):
        ave_vec, ave_vec2=func.accu(np.loadtxt(seed_folder+"/dominant_{}.txt".format(MCO),unpack=True,skiprows=1),ave_vec, ave_vec2)
        
    ave_vec,err_vec=func.norm_and_err(ave_vec,ave_vec2,max_seed)
    np.savetxt(results_folder+"/dominant_{}.txt".format(MCO),np.column_stack((ave_vec[0],err_vec[0],ave_vec[1],err_vec[1],ave_vec[2],err_vec[2])),fmt='%.4f',header="all, error, no_domin, error, only_dom, error ")
        

func.plot_regimes(results_folder,small_b,small_e,"red","small",True)
func.plot_regimes(results_folder,medium_b,medium_e,"springgreen","medium",True)
func.plot_regimes(results_folder,large_b,large_e,"gold","large",True)

func.plot_regimes_vs_exp(results_folder,small_b,small_e,3,5,"darkviolet","small",True,"small $\omega$")
func.plot_regimes_vs_exp(results_folder,medium_b,medium_e,6,9,"hotpink","medium",True, "intermediate $\omega$")
func.plot_regimes_vs_exp(results_folder,large_b,large_e,10,14,"brown","large",True, "high $\omega$")
