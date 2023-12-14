## Useful functions for merging files, sum third cols from different variants of the same dominant MCO
# Francesco Mambretti, 21/11/2022

import numpy as np
import sys

num_of_files=len(sys.argv)-2
print("Merging {} files".format(num_of_files))

tot_y=0

for i in range(0,num_of_files):
    x,y=np.loadtxt(sys.argv[i+1],unpack=True,skiprows=1,usecols=(0,2))
    tot_y+=y

np.savetxt(sys.argv[-1],np.column_stack((x,x-tot_y,tot_y)),fmt="%.4f",header="all \t all-dominant \t dominant ")
