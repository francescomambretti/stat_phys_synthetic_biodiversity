#Written by Francesco Mambretti
#29/04/2022

import numpy as np
from data_blocking import *
import global_var as g
from scipy.optimize import curve_fit

########################################
def init(arguments):
    
    g.traj_file=arguments[1]
    g.steps=float(arguments[2])
    g.blocks=int(arguments[3]) #divide whole traj in blocks
    g.start_nucl=int(arguments[4]) #compute L_P from start_nucl to end_nucl
    g.end_nucl=int(arguments[5]) #indexes follow the overall indexing of the whole system nucleotides
    g.tot_nucl=int(arguments[6]) #tot number of nucleotides

    g.N=g.end_nucl-g.start_nucl+1 #might include or not the blocks at 5' and 3' ends - number of nucleotides to consider

    g.ave_corr=np.resize(g.ave_corr,g.N-1) #set size of np arrays
    g.ave2_corr=np.resize(g.ave2_corr,g.N-1)
    g.glob_ave_corr=np.resize(g.glob_ave_corr,g.N-1)
    g.glob_ave2_corr=np.resize(g.glob_ave2_corr,g.N-1)
    g.error_corr=np.resize(g.error_corr,g.N-1)

    return


########################################
def expo(x,l):  # defining the fitting function
    return np.exp(-x/l)

########################################
def persistence_length(b,split_filename,corr_filename):

    coords = np.zeros((g.N, 3))
    c=0
    corr=np.zeros(g.N-1) #overwritten at each step
    blk_norm=0
    sum = 0
    
    infile=open(split_filename+"_{}.txt".format(b),'r') #Input location of the trajectory file
    g.ave_corr=np.zeros(len(g.ave_corr)) #averages within this block
    g.ave2_corr=np.zeros(len(g.ave2_corr))
    
    for line in infile:
        if (line.split()[0]=='t'):
            #reset
            coords = np.zeros((g.N, 3))
            c=0
        else:
            coords[c] = np.array([ [float(line.split()[0]),float(line.split()[1]),float(line.split()[2])] ])
            c+=1
            
            if (c==g.N-1): #compute scalar products for the current step
                for i in range(0,g.N-1): #for each of the possible nucleotide pair distances
                    sum = 0
                    for j in range(0,g.N-1-i):
                        vec1 = coords[j+1]-coords[j] #distance between consecutive nucleotides
                        vec0 = coords[j+i+1]-coords[j+i] #distance between consecutive nucleotides i bases away
                        sum += np.dot(vec1,vec0)/(np.linalg.norm(vec0)*np.linalg.norm(vec1))
                    corr[i]=sum/(g.N-1-i)
                blk_norm,g.ave_corr,g.ave2_corr=accumulate(blk_norm,corr,g.ave_corr,g.ave2_corr) #for data blocking
    print(blk_norm)
                
    #at the end of the file - i.e., at the end of the block
    averages(blk_norm,g.ave_corr,g.ave2_corr,g.glob_ave_corr,g.glob_ave2_corr,b,corr_filename)
        
    #at the end of each block, fit! The fit is currently done on the progressive global average
    x=np.arange(0,g.N-1)
    p_opt,p_cov = curve_fit(expo,x,g.glob_ave_corr/(b+1),maxfev=1000)
    print("optimized parameter [l] =",p_opt,", parameters uncertainty =",np.sqrt(np.diagonal(p_cov)))
      
    return p_opt
    
