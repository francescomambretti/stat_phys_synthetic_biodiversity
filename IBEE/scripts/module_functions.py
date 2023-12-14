# Written by Francesco Mambretti, 11/07/2022
# custom module with functions used during analysis of Evolutive Algorithm simulations
# 31/07/2022 version

import numpy as np
import os
from Bio import SeqIO
from Bio.Seq import Seq
from more_itertools import locate
from itertools import count, chain
import module_plots as mod_pl
import sys
import pandas as pd
from input_params import *
import subprocess
from difflib import SequenceMatcher

#############################################  accu ##############################################

def accu(array,ave_array, ave2_array): #should work as well for scalars and arrays
    ave_array+=array
    ave2_array+=array*array
    
    return ave_array, ave2_array
    
#############################################  norm_and_err ##############################################

def norm_and_err(ave_array,ave2_array,norm_const):

    ave_array=ave_array/norm_const
    ave2_array=ave2_array/norm_const

    return ave_array, np.sqrt((ave2_array-ave_array*ave_array)/(norm_const))

#############################################  rev_and_compl ##############################################
              
def rev_and_compl(strand): #here, strands with non valid bases are left untouched
 
    my_dna = Seq(strand)
    my_dna = my_dna.reverse_complement()
    rev_strand=str(my_dna)
            
    return rev_strand

############################################  detect_unique_seqs ############################################

def detect_unique_seqs(cycle_folder): #count how many distinct sequences are there in the sample

    l=20 #resource
    
    df=pd.read_csv(cycle_folder+"/list_seqs_all.dat",sep='\t',names=["sequence","MCO"],header=0)
    list_seqs=df['sequence'].values.tolist()
    mco_list=df['MCO'].values.tolist()

    # join them in a 2D list
    complete_list=list(zip(list_seqs,mco_list))
    sorted_list = sorted(complete_list, key=lambda x: x[0])  #sort the list alphabetically according to the first column
        
    n_unique=0
    
    #add headings to file
    file1 = open(cycle_folder+"/unique_all.txt","w") # write results
    file1.write("sequence"+'\t'+"MCO"+'\t'+"abundance"+'\n')

    # count the multiplicity of unique elements of list_strands, and save the strand with the relative overlap
    while (len(sorted_list)>0):
        strand,temp_mco = sorted_list[0] #always take the first element
        copies=0
        
        for elem in sorted_list: #sorted_list[0] is included in this loop
            if (elem==sorted_list[0]):
                copies+=1
            else:
                break
     
        temp_mco=int(temp_mco)
        
        file1.write(strand+'\t'+str(temp_mco)+'\t'+str(copies)+'\n')
            
        # delete them, so to avoid overcounting - and also shorten the total list
        del sorted_list[0:copies]
        
        n_unique+=1
        
    file1.close()
#    print ("In this file, there are {} unique strands".format(n_unique))

    return n_unique

############################################  sort_unique_strands ############################################

def sort_unique_strands (cycle_folder):

    df=pd.read_csv(cycle_folder+"/unique_all.txt",sep='\t')
    df.sort_values(by='abundance', ascending=False, inplace=True, na_position='last')
    df['abundance']=df['abundance'].astype(str)
    df.to_csv(cycle_folder+"/unique_all_sorted.txt",sep='\t',index=False,header="sequence"+'\t'+"MCO"+'\t'+"abundance"+'\n')
    
    #once sorted, unsorted files can be removed to save disk space
    os.system("rm "+cycle_folder+"/unique_all.txt")
    
    return
    
############################################  ave_omega_pop ############################################

def ave_omega_pop (cycle_folder): #compute average omega value of the population at each time
    mco_list,weights=np.loadtxt(cycle_folder+"/unique_all_sorted.txt",unpack=True,usecols=(1,2),skiprows=1)
   
    return np.ma.average(mco_list,weights=weights)
    
############################################  top_n_full_evo ############################################

# tracks the evolution of the n dominant strands at each cycle along all the history
# it selects the list of the top-n for each cycle in folders_list, and tracks their evolution in the past and in the future

def top_n_full_evo (n,seed_folder,seed):

    i=0
    df_list=list() #empty list of pandas dfs - contains the full set of strands
    df_top_n=list() #empty list of pandas dfs - contains the set of top-n strands at each cycle

    for cycle in cycles_IDs: #may be improved, in the future, avoiding to write in sort_unique_strands and read again here?
        
        df=pd.read_csv(seed_folder+"/cycle_{}/unique_all_sorted.txt".format(cycle),sep='\t',dtype = {'sequence': str, 'MCO': np.uint8, 'abundance': np.uint32})
        #do this to reduce memory occupation
        df_list.append(df) #save into list
        df_top_n.append(df[0:n]) # select top-n
        i+=1
                
        if (cycle%freq_cross_MCO==0):
            cross_MCOs_top_n(df[0:n], seed_folder+"/cycle_{}".format(cycle)) #call method for cross-MCOs matrix
    
    a=0 #Incremented by 1 at each cycle
    # now track the evolution of the abundance of the dominant strands of each cycle with the other cycles
    for cycle in cycles_IDs:
        #print("The overall most abundant strands, at cycle {} are:".format(cycle))
        for index, row in df_top_n[a].iterrows(): #loop over df_top_n rows - i.e. select one of the top winners of that cycle, in decreasing abundance order
            howmany=np.zeros(0,dtype=np.float64)
            #print("%d %s %d %.5f" % (index+1, row['sequence'], row['MCO'], row['abundance']/Nr))
            
            for j in cycles_IDs: #select cycle (both backwards and in the future)
                found=False
                for index2, row2 in df_list[j].iterrows(): #loop over df rows
                    if (row2['sequence']==row['sequence']):
                        howmany=np.append(howmany,row2['abundance'])
                        found=True
                        break #I found it, go to the next j
                #try:
                #    row2=df_list[j].loc(df_list[j]['sequence']=="CCGAAAAGCACCGTATGGTCGCCTAGTATATCGGGTGGATCCAATCGTGA") #)row['sequence'])
                #    howmany=np.append(howmany,row2['abundance'])
                #    print(cycle,"FOUND")
                
                if(found==False):
               # except:
                    howmany=np.append(howmany,0)
            
            np.savetxt(seed_folder+"/cycle_{}/top_{}_seq_history.txt".format(cycle,index),np.column_stack((cycles_IDs,howmany/Nr)),fmt="%d %.4f")
                
        a+=1
                
    return

############################################  cross_MCOs_top_n ############################################

# computes the matrix of cross-MCOs among the top-n sequences and prints it to a file

def cross_MCOs_top_n(my_df, specific_folder): #call method for cross-MCOs matrix
        
    cross_MCOs = np.zeros((n+1,n))
    
    #my_df contains the info about the top-n strands of that cycle
    for i1, row1 in my_df.iterrows(): #loop over my_df rows - i.e. select one of the top winners of that cycle, in decreasing abundance order
        cross_MCOs[0][i1]=row1['MCO']
        for i2, row2 in my_df[i1+1:n].iterrows():
            #compute MCO between them
            seq1=row1['sequence']
            seq2=rev_and_compl(row2['sequence'])
            MCO_pair=SequenceMatcher(None, seq1, seq2).find_longest_match(0, len(seq1), 0, len(seq2)).size
            cross_MCOs[i1+1][i2]=MCO_pair
            cross_MCOs[i2+1][i1]=MCO_pair
            
    #save results
    np.savetxt(specific_folder+"/cross_MCOs.txt",cross_MCOs,fmt="%d")

    return cross_MCOs

############################################  compute_ShEn_omega ############################################

def compute_ShEn_omega(input_file,output_ShEn_file):

    histo=np.loadtxt(input_file,usecols=(0,),skiprows=0) #load p(omega) - already sums to 1
    
    #compute Shannon Entropy
    entropy=0
    for bin_height in histo: #bin_height is the probability density associated to each bin
        if(bin_height!=0):
            prob=bin_height
            entropy-=prob*np.log2(prob)
            
    #and update the relative file
    ShEn_file = open(output_ShEn_file,"a") #add the value for each experimental cycle
    string_to_write="{:.4f}".format(entropy)
    ShEn_file.write(string_to_write+'\n')
    ShEn_file.close()

    return histo,entropy

############################################  compute_top_n_fraction ############################################

def compute_top_n_fraction(seed_folder):
    
    fractions=np.zeros(0)
    a=0
    for t in cycles_IDs: #loop over cycles
        frac=0
        for i in range (0,n): #accumulate the n most abundant fractions
            try:
                frac+=np.loadtxt(seed_folder+"/cycle_{}/top_{}_seq_history.txt".format(t,i),skiprows=a,max_rows=1,usecols=(1,)) # column 0 of the file contains cycle id
            except: #in case the file does not exist
                pass
        a+=1
        
        fractions=np.append(fractions,frac)
    np.savetxt(seed_folder+"/top_{}_percentage_over_cycles.txt".format(n),fractions,fmt="%.4f")

    return fractions

############################################  compute_ShEn_RSA ############################################

def compute_ShEn_RSA (cycle_folder,output_ShEn_file):

    df=pd.read_csv(cycle_folder+"/unique_all_sorted.txt",sep='\t',dtype=int, usecols=('abundance',))
    abundances=df['abundance'].to_numpy()

    #compute entropy of the RSA of the whole population
    total=np.sum(abundances) #total sequences
    # compute Shannon Entropy of the RSA - species:= sequence identity
    entropy=0
    for counts in abundances: #bin_height is the probability density associated to each bin
        prob=counts/total
        entropy-=prob*np.log2(prob)

    #and update the relative file
    ShEn_file = open(output_ShEn_file,"a") #add the value for each experimental cycle
    string_to_write="{:.4f}".format(entropy)
    ShEn_file.write(string_to_write+'\n')
    ShEn_file.close()
            
    return entropy

############################################  zip_ratio ############################################

def zip_ratio(file, zip_file):

    #first of all, remove the second column
    data = np.loadtxt(file, usecols=(0,),skiprows=1,dtype='str')
    np.savetxt("all_sequences.dat",data,fmt='%s')

    orig_size = os.path.getsize("all_sequences.dat")
    os.system("zip "+zip_file+" all_sequences.dat")
    zip_size = os.path.getsize(zip_file)
    os.system("rm "+zip_file+" all_sequences.dat")

    return zip_size/orig_size
