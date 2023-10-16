#!/usr/bin/env python
# coding: utf-8
# Code version 06/06/2021
# Sample R times the pi(a) distribution for species abundance  - no threshold
# Francesco Mambretti
# MPI parallelization: each rank manages a given number of repetions, R/size

# import packages

import numpy as np
import scipy as sp
from scipy import stats
from scipy.special import binom
import matplotlib.pyplot as plt
from random import choices
from mpi4py import MPI
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import os

# General parameters

l=20
L=50
bases=4

# Define functions
#########################################################################################################################################################

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

##########################################################################################################################################################

# Create the distribution

appo=0
distrib=np.zeros(l+1)
for a in range (0,l+1):
	appo=(n1(a)+n2(a)+n3(a))/((bases**L)*(L+l-1))
	distrib[a]=appo

#uncomment for direct check
#print(distrib)
#print(np.sum(distrib))

# Setup MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank() 
size = comm.Get_size() 

print("{} of {}".format(rank,size))

# Repeat the sampling from the original distribution R times
R=200

chunks=int(R/size) #assign workload to each rank
if (R%size!=0):
	exit()

# other parameters for sampling
hist_chunk=np.zeros(0)
full_hist=np.zeros(0)
N=10**6 # number of samples for each R
subsample=np.empty(0)
affinities=np.arange(0,l+1,1)
mybins=np.arange(0,l+2,1)

custom = stats.rv_discrete(0,l+1,name='custom', values=(affinities, distrib))

# sample null model probability distribution with physical constraint
for c in range (0,chunks):
	for s in range (0,size):
		if(rank==s):
			subsample=np.empty(0)
			for i in range (0,N):
                var=custom.rvs()
				subsample=np.append(subsample,var)
    
			hist_data,bin_edges=np.histogram(subsample, bins=mybins)
			hist_chunk=comm.gather(hist_data,root=0)

	full_hist=np.append(full_hist,hist_chunk)

# do the plot
if (rank==0):
    print(full_hist)
    
    full_hist=np.reshape(full_hist,(R,l+1))
    full_hist=full_hist.T
    full_hist=full_hist.astype(np.float64)
    full_hist/=(N)

    np.savetxt("full_hist.txt",full_hist,fmt='%.7f')

    for i in range (0,l+1):
        print (i,np.mean(full_hist[i]),np.std(full_hist[i]))

    os.system("python plot_full_hist.py")

MPI.Finalize()
