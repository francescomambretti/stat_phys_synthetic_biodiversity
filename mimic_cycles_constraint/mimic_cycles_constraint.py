# Python code to mimic the evolution of a N_s-sized sample approximately following the \pi'(a) distribution for the affinities at cycle 0.
# iterate the sampling re-defining the probability distribution at each time
# Written by Francesco Mambretti, 28/06/2021

# Update 05/07/2021: add MPI feature to simulate many repetitions of the simulation

# OK version

import numpy as np
import scipy as sp
from scipy import stats
from scipy.special import binom
import matplotlib.pyplot as plt
from random import choices
from mpi4py import MPI
import os

#define 

l=20
L=50
bases=4

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print("rank {} of {}".format(rank,size))

if not os.path.exists("run_{}".format(rank)):
	os.mkdir("run_{}".format(rank)) # create a separate folder for each rank

########################## FUNCTIONS DEFINITIONs ####################################################

def my_binom (a,b):
    if (a<0 or a<b):
        return 0
    else:
        return binom(a,b)
        
def n1_part (d,a):
    list1=np.empty(0)
    list2=np.empty(0)
    for r in range (1,d+2):
        list1=np.append( list1, (np.power(-1,r))*(my_binom(d+1,r))*(my_binom(l-r*(a+1),d)) )
        list2=np.append( list2, (np.power(-1,r))*(my_binom(d+1,r))*(my_binom(l-r*a,d)) )
    return np.sum(list1)-np.sum(list2)
    
def n1 (a):
    list3=np.empty(0)
    for d in range (0,l+1):
        list3=np.append(list3, ((bases-1)**d) * n1_part(d,a))
    return (L-l+1)*(bases**(L-l))*(np.sum(list3))
    
def n2_part (d,a,p):
    list1=np.empty(0)
    list2=np.empty(0)
    for r in range (1,d+2):
        list1=np.append( list1, np.power(-1,r)*(my_binom(d+1,r))*(my_binom(l+p-r*(a+1),d)) )
        list2=np.append( list2, np.power(-1,r)*(my_binom(d+1,r))*(my_binom(l+p-r*a,d)) )
    return np.sum(list1)-np.sum(list2)

def n2(a):
    list4=np.empty(0)
    for p in range (-l+1,0):
        list3=np.empty(0)
        for d in range (0,l+p+1):
            list3=np.append(list3, ((bases-1)**d) * n2_part(d,a,p))
        list4=np.append(list4, ((bases**(L-l-p))*np.sum(list3)))
            
    return np.sum(list4)
    
def n3_part (d,a,p):
    list1=np.empty(0)
    list2=np.empty(0)
    for r in range (1,d+2):
        list1=np.append( list1, np.power(-1,r)*(my_binom(d+1,r))*(my_binom(L-p-r*(a+1),d)) )
        list2=np.append( list2, np.power(-1,r)*(my_binom(d+1,r))*(my_binom(L-p-r*a,d)) )
    return np.sum(list1)-np.sum(list2)
    
def n3(a):
    list4=np.empty(0)
    for p in range (L-l+1,L):
        list3=np.empty(0)
        for d in range (0,L-p+1):
            list3=np.append(list3, np.power(bases-1,d) * n3_part(d,a,p))
        list4=np.append(list4, (bases**(p))*(np.sum(list3)))
            
    return np.sum(list4)
    
#compute pi(a) - exact formula

appo=0
norm=0
pi_dist=np.zeros(l+1)

#compute normalization
for a in range (0,l+1):
    norm+=n1(a)+n2(a)+n3(a)

for a in range (0,l+1):
    appo=(n1(a)+n2(a)+n3(a))/((bases**L)*(L+l-1))
    pi_dist[a]=appo

# compute pi'(a) - exact formula
appo=0
norm=0
T=3
pi_pr_dist=np.zeros(l+1)

for a in range (T,l+1):
    norm+=pi_dist[a]

for a in range (T,l+1):
    pi_pr_dist[a]=(pi_dist[a]*(np.sign(a-T+1)))/norm

########################## CYCLES ####################################################

### choose among ranks
for s in range (0,size):
    if(rank==s):
        cycles=9 #can be changed, now it is 9 because we currently have 9 cycles in the experiments
        temp_dist=pi_pr_dist #at cycle 0
        affinities=np.arange(0,l+1,1)
        mybins=np.arange(0,l+2,1)
        N_s=int(1E6) #sample size
        custom = stats.rv_discrete(T,l+1,name='custom', values=(affinities, temp_dist))
    
        for c in range (0,cycles):
            print(c)
            #sample from the distribution numerically
            subsample=np.empty(0)
    
            for i in range (0,N_s):
                var=custom.rvs()
                subsample=np.append(subsample,var)
    
            hist_data,bin_edges=np.histogram(subsample, bins=mybins) # create affinity histogram
            # now this histogram becomes the new probability distribution
            temp_dist= hist_data/np.sum(hist_data)
            custom = stats.rv_discrete(T,l+1,name='custom', values=(affinities,temp_dist))
            if ((c+1)%3==0):
                np.savetxt("run_{}/cycle.{}.txt".format(s,c+1),temp_dist)

MPI.Finalize()
