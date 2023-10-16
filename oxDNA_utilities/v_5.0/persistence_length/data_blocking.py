#Written by Francesco Mambretti
#29/04/2022

import numpy as np
import global_var as g

########################################
## data blocking ##

def accumulate(blk_norm,array,ave_array,ave2_array): #at the end of each step in a block

    ave_array+=array #update block averages
    ave2_array+=array*array
    
    blk_norm+=1
    
    return blk_norm,ave_array,ave2_array
    
def averages (blk_norm,ave_array,ave2_array,glob_ave_array,glob_ave2_array,i_blk,filename):

    ave_array/=blk_norm #normalization - averaged over frames in the block
    ave2_array/=blk_norm

    glob_ave_array+=ave_array
    glob_ave2_array+=ave_array*ave_array
    
    length=len(ave_array)
    
    glob_error_array=np.zeros(length)
    glob_error_array=error(glob_ave_array,glob_ave2_array,i_blk)
    
    #print on file progressive estimation of global average and related error
    with open (str(filename)+"{}.txt".format(i_blk),"a") as outfile:
        for i in range (0,length):
            outfile.write(str(ave_array[i])+" "+str(glob_ave_array[i]/(i_blk+1))+" "+str(glob_error_array[i])+"\n")
    
    print(outfile.closed)
    
    return

def error(glob_ave_array,glob_ave2_array,i_blk): #error on progressive (global) average
    length=len(glob_ave_array)
    error_array=np.zeros(length)
    
    if i_blk==0:
        return error_array
    else:
        for i in range(0,length):
            error_array[i]=np.sqrt((glob_ave2_array[i]/(i_blk+1) - (glob_ave_array[i]/(i_blk+1))**2)/(i_blk))
            
        return error_array
