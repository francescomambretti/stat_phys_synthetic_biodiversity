# Written by Francesco Mambretti
# split trajectory in blocks
# save only the coordinates of the selected nucleotides
# compute the correlation array
# extract persistence length L_P by fitting
# L_P estimate is update according to global average of the tangent vectors correlation array
# 29/04/2022

import os
import sys
import numpy as np

import global_var as g
from data_blocking import *
from functions_pers_length import *
from functions_traj import *
from mod_plot import *

def clean (filelist):
    for file in filelist:
        os.system ("rm "+file)

    return
    
########################################
#main
########################################

g.defining_global_variables()

#read command line
if len(sys.argv) < 7:
    print("Error! Correct usage: "+sys.argv[0]+" trajectoryFile steps blocks start_nucl end_nucl tot_nucl")
    sys.exit()

init(sys.argv) #assign input arguments to global variables

if (g.steps%g.blocks!=0):
    print("Error! The number of blocks has to be consistent with the number of steps!")
    sys.exit()
else:
    g.steps_per_block=g.steps/g.blocks
    print("The trajectory of {} steps has been divided into {} blocks with {} steps in each block".format(g.steps,g.blocks,g.steps_per_block))

list_breaks=find_breaks(g.steps,g.blocks)

split_filename="traj_split"
corr_filename="corr_split"
png_filename="corr_blk"
clean((split_filename+"*",corr_filename+"*"))

read_traj(split_filename,list_breaks,g.start_nucl,g.end_nucl)

#now compute the correlation array on each splitted trajectory file
#and fit its decay to extract the correlation length

for b in range (0,g.blocks):
    L_P=persistence_length(b,split_filename,corr_filename)
    #plot correlation function of each block & the progressive average
    plot_blk_glb_ave(b,png_filename+"_{}.png".format(b));
