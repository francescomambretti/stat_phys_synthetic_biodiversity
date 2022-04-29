# Written by Francesco Mambretti
# declare here all global variables
# 20/04/2022

import numpy as np

def defining_global_variables():
    global steps_per_block, traj_file, steps, blocks, start_nucl, end_nucl, tot_nucl, N, ave_corr, ave2_corr, glob_ave_corr, glob_ave2_corr, error_corr, glob_ave_L_P, glob_ave2_L_P, error_L_P

    steps_per_block=1
    traj_file=""
    steps=0
    blocks=0
    start_nucl=0
    end_nucl=0
    tot_nucl=0
    N=0
    glob_ave_corr=np.zeros(0)
    glob_ave2_corr=np.zeros(0)
    error_corr=np.zeros(0)
    ave_corr=np.zeros(0)
    ave2_corr=np.zeros(0)
    glob_ave_L_P=0
    glob_ave2_L_P=0
    error_L_P=0

    return
