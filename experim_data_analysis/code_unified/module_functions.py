# Written by Francesco Mambretti, 01/02/2022
# custom module with functions used during experimental FASTQ files analysis
# 21/03/2022 version

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

#############################################  rev_and_compl ##############################################
              
def rev_and_compl(strand): #here, strands with non valid bases are left untouched
 
    my_dna = Seq(strand)
    my_dna = my_dna.reverse_complement()
    rev_strand=str(my_dna)
            
    return rev_strand
    
#############################################  filter ##############################################
              
def filter(sequence_temp): #decide whether to accept or not a strand according to special criterion
#acts either on forward strands or on reverse strands, reverted and complemented

    if(key_filter==False):
        return True
    else: #choose the criterion at will, according to your needs
        #MCO_alien_1=subprocess.check_output("./find_equal_pair.x "+sequence_temp+" "+fixed_end_5prime, shell=True)
        MCO_alien_1=SequenceMatcher(None, fixed_end_5prime, sequence_temp).find_longest_match(0, len(fixed_end_5prime), 0, len(sequence_temp)).size
        
        if(int(MCO_alien_1)>=alien_bases): #here, ">=" means: if you have at least 'alien_bases' bases identifying you as an alien, you are an alien and you are stopped by the filter
            return False
    
        else: # you may as well be an alien including part of fixed_end_3prime
            #MCO_alien_2=subprocess.check_output("./find_equal_pair.x "+sequence_temp+" "+fixed_end_3prime, shell=True)
            MCO_alien_2=SequenceMatcher(None, fixed_end_3prime, sequence_temp).find_longest_match(0, len(fixed_end_3prime), 0, len(sequence_temp)).size
            if(int(MCO_alien_2)<alien_bases):
                return True # both criteria are satisfied, you are not an alien!
            else:
                return False
    
#############################################  analyze_partial ##############################################
              
def analyze_partial(sequence_temp, length, mycount_NRNF, mycount_p, my_tot_valid, myfile, my_start, my_size): #analyze or forward or reverse strands
        
    if (sequence_temp.startswith(my_start)): # select either forward or reverse sequences, according to key3
        sequence_temp=sequence_temp[my_size:]
       
        if(key3=="rev"):
            sequence_temp=rev_and_compl(sequence_temp)
        
        if (filter(sequence_temp)==True): #apply filter only after having reversed and complemented
            mycount_p+=1
            my_tot_valid+=1
            
            myfile.write(sequence_temp + "\n")

    else:
        mycount_NRNF+=1

    return mycount_NRNF, mycount_p, my_tot_valid
    
#############################################  analyze_all ##############################################
              
def analyze_all(sequence_temp,length,mycount_NRNF,mycountF,mycountR,my_tot_valid,myfile1,myfile2,myfile3): #analyze both forward and reverse strands
 
    if (sequence_temp.startswith(fw_primer)): # select either forward or reverse sequences, according to key3
        sequence_temp=sequence_temp[fw_primer_size:]
        
        if (filter(sequence_temp)==True):
            mycountF+=1
            my_tot_valid+=1

            myfile1.write(sequence_temp + "\n")
            myfile2.write(sequence_temp + "\n")
        
    elif (sequence_temp.startswith(rev_primer)): # select reverse sequences
        sequence_temp=sequence_temp[rev_primer_size:]
        
        #reverse and complement: important! To be done before removing primer
        sequence_temp=rev_and_compl(sequence_temp)
        
        if (filter(sequence_temp)==True):
            mycountR+=1
            my_tot_valid+=1

            myfile1.write(sequence_temp + "\n")
            myfile3.write(sequence_temp + "\n")
                
    else:
        mycount_NRNF=mycount_NRNF+1

    return mycount_NRNF,mycountF,mycountR,my_tot_valid

############################################  process_FASTQ ############################################

def process_FASTQ(*args,mother_folder,results_folder,stop,use_stop,t): #args is used because of variable number of arguments

    if key2=="R1R2":
        fastq_file_R1=args[0]
        fastq_file_R2=args[1]
    else:
        fastq_file=args[0]

    use_stop=int(use_stop)

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    
    if key2=="R1R2":
        specific_results=results_folder+"/"+str(t)+"_"+fastq_file_R1+"R2"
    else:
        specific_results=results_folder+"/"+str(t)+"_"+fastq_file
        
    if not os.path.exists(specific_results):
        os.makedirs(specific_results)
        
    file1 = open(specific_results+"/all.txt","w")
    if (key3=="all"):
        file2 = open(specific_results+"/fw.txt","w")
        file3 = open(specific_results+"/rev.txt","w")
    elif (key3=="fw" or key3=="rev"):
        file1 = open(specific_results+"/all.txt","w") # no need for separate files
        
    count_NRNF=0 #sequences neither forward nor reverse
    countN=0 #sequences with non valid entries
    countF=0 #sequences forward
    countR=0 #sequences reverse
    count_lowQ=0 #count number of low-Q reads
    tot=0
    tot_valid=0
    
    count_length_extra=0 # count strands whose length is smaller than lower_bound and larger than L
    
    counter=-1 #use this only to limit the analysis to a subset of the strands
    
    #build the distribution of strand lengths
    dist_lengths=np.zeros(L+1)

    if (key2=="R1R2"):
        seq_list_R1=SeqIO.parse(mother_folder+"/"+fastq_file_R1+".fastq", "fastq")
        seq_list_R2=SeqIO.parse(mother_folder+"/"+fastq_file_R2+".fastq", "fastq")
        container=chain(seq_list_R1,seq_list_R2)
    else:
        container=SeqIO.parse(mother_folder+"/"+fastq_file+".fastq", "fastq")

    for seq_record in container: #SeqIO.parse(mother_folder+"/"+fastq_file+".fastq", "fastq"):
        counter+=use_stop # True=1, False=0
        if (counter<stop): #this check is always true if use_stop==False, because counter is not updated
        
            score=seq_record.letter_annotations["phred_quality"] #list of the quality score of each base
            if(min(score)>=min_Q): #else, the quality of the read is too low
        
                sequence_temp=str(seq_record.seq[skip_bases:]) # write all sequences on this file
                length=len(sequence_temp)
                
                if (lower_bound<=length):
                    dist_lengths[length]+=1 # update the histogram of the distribution of the strand lengths
                    if ('N' not in sequence_temp): #check the possible presence of 'N' nucleotides and remove those strings
                        if(key3=="all"):
                            count_NRNF,countF,countR,tot_valid=analyze_all(sequence_temp,length,count_NRNF,countF,countR,tot_valid,file1,file2,file3)
                        
                        elif(key3=="fw"):
                            count_NRNF,countF,tot_valid=analyze_partial(sequence_temp,length,count_NRNF,countF,tot_valid,file1,fw_primer,fw_primer_size)
                    
                        elif(key3=="rev"):
                            count_NRNF,countR,tot_valid=analyze_partial(sequence_temp,length,count_NRNF,countR,tot_valid,file1,rev_primer,rev_primer_size)
                            
                    else:
                        countN+=1
                    
                else:
                    count_length_extra+=1
                    
                tot+=1
                
            else: #the Q of the string is too low
                count_lowQ+=1
        else:
            break  #interrupt, counter >= N
        
    #final percentages
    if (key2=="R1R2"):
        print ("FASTQ file -->"+fastq_file_R1+" and "+fastq_file_R2)
    else:
        print ("FASTQ file -->"+fastq_file)
    
    print ("Total number of strands: "+str(tot)+", where we have "+str(tot_valid)+" valid strands")
    if(key3=="all" or key3=="fw"):
        print (str(countF)+", {:.2f}".format(countF/tot*100)+"% sequences are forward (i.e. {:.2f}% of the valid ones)".format(countF/tot_valid*100))
        
    if(key3=="all" or key3=="rev"):
        print (str(countR)+", {:.2f}".format(countR/tot*100)+"% sequences are reverse (i.e. {:.2f}% of the valid ones)".format(countR/tot_valid*100))
        
    print (str(count_NRNF)+", {:.2f}".format(count_NRNF/tot*100)+"% sequences are neither reverse nor forward")
    print (str(countN)+", {:.2f}".format(countN/tot*100)+"% sequences have non valid bases/are aliens")
    print (str(count_length_extra)+", {:.2f}".format(count_length_extra/tot*100)+"% sequences have a length smaller than {} or larger than {}".format(lower_bound,L))
    print (str(count_lowQ)+" strings (corresponding to {:.3f} % of the valid ones) have a too low quality".format(count_lowQ/tot_valid))

    file1.close()
    
    if (key3=="all"):
        file2.close()
        file3.close()
    
    #save the histogram of the strands lengths distribution
    np.savetxt(specific_results+"/dist_lengths.txt", dist_lengths, fmt='%d')

    target_rc=rev_and_compl(str(target)); # necessary for the new find_MCO_serial version

    os.system("./find_MCO_serial.x "+specific_results+"/all.txt "+str(tot_valid)+" "+specific_results+"/all_MCO.txt "+str(target_rc))
    
    if (key3=="all"): #do also separate analysis
        os.system("./find_MCO_serial.x "+specific_results+"/fw.txt "+str(countF)+" "+specific_results+ "/fw_MCO.txt "+str(target_rc))
        os.system("./find_MCO_serial.x "+specific_results+"/rev.txt "+str(countR)+" "+specific_results+"/rev_MCO.txt "+str(target_rc))
        
    return tot_valid
    
######################################## create_separate_MCO_folders #########################################

def create_separate_MCO_folders(l,specific_results):
        
    for m in range (0,l+1): #resource length = l
        if not os.path.exists(specific_results+"/MCO_{}".format(m)):
            os.makedirs(specific_results+"/MCO_{}".format(m))
            
    return

############################################  detect_unique_strands ############################################

def detect_unique_strands(specific_results): #count how many distinct sequences are there in the sample

    l=20 #resource
    
    df=pd.read_csv(specific_results+"/all_MCO.txt",sep='\t',dtype=np.uint8,usecols=('MCO','MCO_2nd','TMO','LTO'))
    mco_list=df['MCO'].values.tolist()
    mco_2nd_list=df['MCO_2nd'].values.tolist()
    tmo_list=df['TMO'].values.tolist()
    lto_list=df['LTO'].values.tolist()
    
    with open(specific_results+"/all.txt") as f:
        list_strands = f.read().splitlines()
        
    # join them in a 2D list
    complete_list=list(zip(list_strands,mco_list,mco_2nd_list,tmo_list,lto_list))
    sorted_list = sorted(complete_list, key=lambda x: x[0])  #sort the list alphabetically according to the first column
    f.close()
        
    n_unique=0
    
    #add headings to file
    file1 = open(specific_results+"/unique_all.txt","w") # write results
    file1.write("strand"+'\t'+"MCO"+'\t'+"MCO_2nd"+'\t'+"TMO"+'\t'+"LTO"+'\t'+"abundance"+'\n')
    
    for temp_mco in np.arange(0,l+1):
        file2 = open(specific_results+"/MCO_{}/unique_MCO.txt".format(temp_mco),"a")
        file2.write("strand"+'\t'+"MCO"+'\t'+"MCO_2nd"+'\t'+"TMO"+'\t'+"LTO"+'\t'+"abundance"+'\n')
        file2.close()

    # count the multiplicity of unique elements of list_strands, and save the strand with the relative overlap
    while (len(sorted_list)>0):
        strand,temp_mco,temp_mco_2nd,temp_tmo,temp_lto = sorted_list[0] #always take the first element
        copies=0
        #positions=[(i, j) for i, j in zip(count(), sorted_list) if j[0] == strand]
        for elem in sorted_list: #sorted_list[0] is included in this loop
            if (elem==sorted_list[0]):
                copies+=1
            else:
                break
        # alternative, slower methods
        #positions=list(locate(sorted_list, lambda x: x[0]==strand)) # returns the list of positions
        #positions = [i for i in range(len(sorted_list)) if sorted_list[i] == (strand,temp_mco)]
        
        #copies=len(positions) # we have this number of copies of that strand
        temp_mco=int(temp_mco)
        temp_mco_2nd=int(temp_mco_2nd)
        temp_tmo=int(temp_tmo)
        temp_lto=int(temp_lto)
        
        file1.write(strand+'\t'+str(temp_mco)+'\t'+str(temp_mco_2nd)+'\t'+str(temp_tmo)+'\t'+str(temp_lto)+'\t'+str(copies)+'\n')
        
        # and now split the data among the MCO folders
        file2 = open(specific_results+"/MCO_{}/unique_MCO.txt".format(temp_mco),"a")
        file2.write(strand+'\t'+str(temp_mco)+'\t'+str(temp_mco_2nd)+'\t'+str(temp_tmo)+'\t'+str(temp_lto)+'\t'+str(copies)+'\n')
            
        # delete them, so to avoid overcounting - and also shorten the total list
        del sorted_list[0:copies]
        #sorted_list.remove((strand,temp_mco)) #slower
        
        n_unique+=1
        
    file1.close()
    file2.close()
    print ("In this file, there are {} unique strands".format(n_unique))

    return n_unique

############################################  sort_unique_strands ############################################

def sort_unique_strands (specific_results):
    l=20
    for temp_mco in np.arange(0,l+1):
        df=pd.read_csv(specific_results+"/MCO_{}/unique_MCO.txt".format(temp_mco),sep='\t')
        df['abundance']=df['abundance'].astype(str)
        df.sort_values(by='abundance', ascending=False, inplace=True, na_position='last')
        df.to_csv(specific_results+"/MCO_{}/unique_MCO_sorted.txt".format(temp_mco),sep='\t',index=False,header="strand"+'\t'+"MCO"+'\t'+"abundance"+'\n')
        #once sorted, unsorted files can be removed to save disk space
        os.system("rm "+specific_results+"/MCO_{}/unique_MCO.txt".format(temp_mco))

    df=pd.read_csv(specific_results+"/unique_all.txt",sep='\t')
    df.sort_values(by='abundance', ascending=False, inplace=True, na_position='last')
    df['abundance']=df['abundance'].astype(str)
    df.to_csv(specific_results+"/unique_all_sorted.txt",sep='\t',index=False,header="strand"+'\t'+"MCO"+'\t'+"abundance"+'\n')
    
    #once sorted, unsorted files can be removed to save disk space
    os.system("rm "+specific_results+"/unique_all.txt")
    
    return

############################################  top_n_full_evo ############################################

# it works for unique_all_sorted.txt and for unique_MCO_sorted.txt files in MCO_* folders
# tracks the evolution of the n dominant strands at each cycle along all the history
# it selects the list of the top-n for each cycle in folders_list, and tracks their evolution in the past and in the future

def top_n_full_evo (*args,n,results_folder,param,tot_valid):
    i=0
    df_list=list() #empty list of pandas dfs - contains the full set of strands
    df_top_n=list() #empty list of pandas dfs - contains the set of top-n strands at each cycle

    if (key2=="R1" or key2=="R2"):
        fastq_file_list=args[0]
        length=len(fastq_file_list)
        print(length)
        for fastq_file in fastq_file_list: #may be improved, in the future, avoiding to write in sort_unique_strands and read again here?
            full_folder=results_folder+"/"+str(i)+"_"+fastq_file
            if (param!=-1):
                df=pd.read_csv(full_folder+"/MCO_{}/unique_MCO_sorted.txt".format(param),sep='\t',
                dtype = {'strand': str, 'MCO': np.uint8, 'MCO_2nd': np.uint8, 'TMO': np.uint8, 'LTO': np.uint8, 'abundance': np.uint32})
            else:
                df=pd.read_csv(full_folder+"/unique_all_sorted.txt",sep='\t',dtype = {'strand': str, 'MCO': np.uint8, 'MCO_2nd': np.uint8, 'TMO': np.uint8, 'LTO': np.uint8, 'abundance': np.uint32})
            #do this to reduce memory occupation
            df_list.append(df) #save into list
            df_top_n.append(df[0:n]) # select top-n
            i+=1
                    
            cross_MCOs_top_n(df[0:n], fastq_file, full_folder) #call method for cross-MCOs matrix
        
    else:
        fastq_file_list_R1=args[0]
        fastq_file_list_R2=args[1]
        length=len(fastq_file_list_R1)
        print(length)
        for fastq_file in fastq_file_list_R1: #may be improved, in the future, avoiding to write in sort_unique_strands and read again here?
            full_folder=results_folder+"/"+str(i)+"_"+fastq_file+"R2"
            if (param!=-1):
                df=pd.read_csv(full_folder+"/MCO_{}/unique_MCO_sorted.txt".format(param),sep='\t',
                dtype = {'strand': str, 'MCO': np.uint8, 'MCO_2nd': np.uint8, 'TMO': np.uint8, 'LTO': np.uint8, 'abundance': np.uint32})
            else:
                df=pd.read_csv(full_folder+"/unique_all_sorted.txt",sep='\t',
                dtype = {'strand': str, 'MCO': np.uint8, 'MCO_2nd': np.uint8, 'TMO': np.uint8, 'LTO': np.uint8, 'abundance': np.uint32})
            #do this to reduce memory occupation
            df_list.append(df) #save into list
            df_top_n.append(df[0:n]) # select top-n
            i+=1
            
            cross_MCOs_top_n(df[0:n], fastq_file, full_folder) #call method for cross-MCOs matrix
        
    a=0 #Incremented by 1 at each cycle
    # now track the evolution of the abundance of the dominant strands of each cycle with the other cycles
    for ID in cycles_IDs:
        if (param!=-1):
            print("The most abundant strands, for MCO={}, at cycle {} are:".format(param,ID))
        else:
            print("The overall most abundant strands, at cycle {} are:".format(ID))
    
        for index, row in df_top_n[a].iterrows(): #loop over df_top_n rows - i.e. select one of the top winners of that cycle, in decreasing abundance order
            howmany=np.zeros(0,dtype=np.float64)
            print("%d %s %d %d %d %d %.5f" % (index+1, row['strand'], row['MCO'], row['MCO_2nd'], row['TMO'], row['LTO'], row['abundance']/tot_valid[a]))
            
            for j in range (0,length): #select cycle (both backwards and in the future)
                found=False
                for index2, row2 in df_list[j].iterrows(): #loop over df rows
                    if (row2['strand']==row['strand']):
                        howmany=np.append(howmany,row2['abundance']/tot_valid[j])
                        found=True
                        break #I found it, go to the next j
                if(found==False):
                    howmany=np.append(howmany,0)
            
            
            if key2=="R1R2":
                full_folder=results_folder+"/"+str(a)+"_"+fastq_file_list_R1[a]+"R2"
            else:
                full_folder=results_folder+"/"+str(a)+"_"+fastq_file_list[a]
                
            if (param==-1):
                np.savetxt(full_folder+"/top_strand_{}_cycle_{}.txt".format(index,ID),np.column_stack((cycles_IDs,howmany)),fmt="%.3e")

            else:
                if not os.path.exists(full_folder+"/MCO_{}/".format(param)):
                    os.makedirs(full_folder+"/MCO_{}/".format(param))
                    
                np.savetxt(full_folder+"/MCO_{}/top_strand_{}_cycle_{}.txt".format(param,index,ID),np.column_stack((cycles_IDs,howmany)),fmt="%.3e")
                
        a+=1
                
    return


############################################  cross_MCOs_top_n ############################################

# computes the matrix of cross-MCOs among the top-n sequences and prints it to a file
                    
def cross_MCOs_top_n(my_df, fastq_file, full_folder): #call method for cross-MCOs matrix
        
    cross_MCOs = np.zeros((n+1,n))
    
    #my_df contains the info about the top-n strands of that cycle
    for i1, row1 in my_df.iterrows(): #loop over my_df rows - i.e. select one of the top winners of that cycle, in decreasing abundance order
        cross_MCOs[0][i1]=row1['MCO']
        for i2, row2 in my_df[i1+1:n].iterrows():
            #compute MCO between them
            seq1=row1['strand']
            seq2=rev_and_compl(row2['strand'])
            MCO_pair=SequenceMatcher(None, seq1, seq2).find_longest_match(0, len(seq1), 0, len(seq2)).size
            cross_MCOs[i1+1][i2]=MCO_pair
            cross_MCOs[i2+1][i1]=MCO_pair
            
    #save results
    np.savetxt(full_folder+"/cross_MCOs.txt",cross_MCOs,fmt="%d")

    return
