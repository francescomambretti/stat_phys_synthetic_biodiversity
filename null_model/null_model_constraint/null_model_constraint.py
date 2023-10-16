#!/usr/bin/env python
# coding: utf-8
# Code version 24/06/2021 - by Francesco Mambretti
# Sample R times the analytic pi'(a) distribution for species abundance - null model + constraint on minimum possible overlap (threshold T)
#
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
import sys

# read number of sampling repetitions R from command line
if len(sys.argv)==2:
    print ("Argument list:", str(sys.argv))
else:
    print("Wrong argument list! Use: ", sys.argv[0]," sampling repetitions R")
    
# Setup MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print("rank {} of {}".format(rank,size))

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

# Create the pi(a) distribution

appo=0
pi_dist=np.zeros(l+1)
for a in range (0,l+1):
	appo=(n1(a)+n2(a)+n3(a))/((bases**L)*(L+l-1))
	pi_dist[a]=appo
 
# observation: normalization can be also computed as:
# norm = 0
# for a in range (0,l+1):
#    norm+=n1(a)+n2(a)+n3(a)
# try to believe, they are equal

# uncomment for direct check
# print(pi_dist)
# print(np.sum(pi_dist))

T=3
norm=0
appo=0
pi_pr_dist=np.zeros(l+1)

for a in range (T,l+1):
    norm+=pi_dist[a]

for a in range (T,l+1):
    pi_pr_dist[a]=(pi_dist[a]*(np.sign(a-T+1)))/norm # note: the sign function in Numpy is 0 in 0
    
# print(pi_pr_dist)
# print(np.sum(pi_pr_dist))

##########################################################################################################################################################

# Repeat the sampling from the original distribution R times
R=int(sys.argv[1])

chunks=int(R/size) #assign workload to each rank
if (R%size!=0):
	exit()

# other parameters for sampling
hist_chunk=np.zeros(0)
hist_chunk_2=np.zeros(0)
full_hist=np.zeros(0)   #for n_i
full_hist_2=np.zeros(0)   # for n_i^2
N=10**6 # number of samples for each R
#threshold=3 # discard draws with a<threshold and re-sample
subsample=np.empty(0)
affinities=np.arange(0,l+1,1)
mybins=np.arange(0,l+2,1)

custom = stats.rv_discrete(0,l+1,name='custom', values=(affinities, pi_pr_dist), seed=rank)

# sample null model probability distribution with physical constraint
for c in range (0,chunks):
    print("chunk # {}".format(c))
    for s in range (0,size):
        if(rank==s):
            subsample=np.empty(0)
            subsample_2=np.empty(0)
            for i in range (0,N):
                var=custom.rvs()  #extract a species index - i.e. a value of the maximum consecutive overlap a
                subsample=np.append(subsample,var)
                # the following lines are meant not to sample directly from pi'(a), but for empirically sampling from pi(a) and rejecting all the variables extracted smaller than T
                #            var=0
                #            subsample=np.empty(0)
                #            for i in range (0,N):
                #                while (var < threshold):
                #                    var=custom.rvs()
                #                subsample=np.append(subsample,var)
                #                var=0
                
            hist_data,bin_edges=np.histogram(subsample, bins=mybins)  #estimation of p_i*N
            hist_chunk=comm.gather(hist_data,root=0)
            hist_data_2=hist_data*(hist_data-1) #estimation of (p_i^2)*N*(N-1)
            hist_chunk_2=comm.gather(hist_data_2,root=0)
            #print(hist_data,hist_data-1,hist_data_2)

    full_hist=np.append(full_hist,hist_chunk)
    full_hist_2=np.append(full_hist_2,hist_chunk_2)

if (rank==0):
    #print(full_hist)
    
    full_hist=np.reshape(full_hist,(R,l+1))
    full_hist=full_hist.T
    full_hist=full_hist.astype(np.float64)
    full_hist_2=np.reshape(full_hist_2,(R,l+1))
    full_hist_2=full_hist_2.T
    full_hist_2=full_hist_2.astype(np.float64)

    print("$\langle p_i \rangle$ with errorbars is estimated as: ")
    for i in range (0,l+1):
        mean=np.mean(full_hist[i])/N  #<p_i>
        mean2=np.mean(full_hist_2[i])/(N*(N-1)) #<p_i^2>
        print (i,mean,np.sqrt(np.fabs(mean2-mean*mean)))
        
    full_hist/=(N)

    np.savetxt("full_hist.txt",full_hist,fmt='%.7f')
    
# do the plot
    os.system("python plot_full_hist_and_cumul.py 0/1") #0: plot all the marginals, 1: cumulative probability

MPI.Finalize()
