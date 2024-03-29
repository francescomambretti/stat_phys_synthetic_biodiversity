# set input options
# 28/11/2023 version

import numpy as np

#################################### general options, keys to be set ####################################

key1="oligo1" #"oligo1", "oligo2", "negative" or "seriesN" - "$KEY1$"
key2="R1" # "R1", "R2", "R1R2" - "$KEY2$"
key3="all" # "fw", "rev" or "all" - "$KEY3$"
key_filter=True #if True, call a dedicated function inside module_functions to apply further filtering on data - $KEY_FILTER$
key_no_cut=False #if True, do not remove primer sequences when printing sequences - $KEY_NO_CUT$

#################################### plots options, to be set ####################################
### parameters for plots
keywords=("all", "fw", "rev") #some of them can be commented
colornames=("limegreen","lightskyblue","sandybrown") #comment here accordingly - colors for omega histograms, all,fw,rev

create_png_abund=True # create png files for abundance of strands inside MCO folders

#################################### analysis parameters ####################################

min_Q=10 #lowest phred quality accepted
l=20 #resource length
subset_steps=4e5 # analyze only the first subset_steps
use_stop=True # decide whether to really do it
n=10 #10 #top-n strands

# info about strands
random_seq=50 #50 random nucleotides, by default
cap_size=25
extra_end=1
L=150#random_seq+cap_size+extra_end #max, with cap and last one - length of predators
lower_bound=70#L-6
upper_bound=145

alien_bases=8 #number of bases for the MCO analysis about alien sequences

#################################### variables set inside if conditions ####################################

fastq_file_list=()
fastq_file_list_plot=()
fastq_file_list_R1=()
fastq_file_list_R2=()
cycles_IDs=np.empty(0)
mother_folder=""
results_folder="./" #"../re-submission/results_bootstrap/9old_20x1E5/$MY_REP/" #"../re-submission/results_bootstrap/rep_0/"
fw_primer=""
fw_full_fixed_seq=""
rev_primer=""
rev_full_fixed_seq=""
fixed_end_5prime=""
fixed_end_3prime=""
full_fw_primer_size=0
full_rev_primer_size=0
fw_primer_size=0
rev_primer_size=0
target=""

#################################### folders, fastq list, cycles IDs ####################################

if key1=="oligo1":
    
    target="CGGTATTGGACCCTCGCATG"

    mother_folder="/home/fmambretti@iit.local/Documents/padova/re-submission/data_oligo1"  #"../data/oligo1_proto06/" # where to take FASTQ files from
    if (key_filter==False and key_no_cut==False):
        results_folder=results_folder+key2+"_"+key3+"/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==False):
        results_folder=results_folder+key2+"_"+key3+"_key_filter/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==True):
        results_folder=results_folder+key2+"_"+key3+"_key_filter_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==False and key_no_cut==True):
        results_folder=results_folder+key2+"_"+key3+"_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    
    if key2=="R1" or key2=="R2":
#        fastq_file_list=("oligo1_R8_S25_"+key2+"_001","oligo1_R9_S26_"+key2+"_001","oligo1_R9bis_S27_"+key2+"_001","oligo1_R10_S28_"+key2+"_001") #("oligo1_"+key2,"R1_06_"+key2,"R2_06_"+key2,"R3-06_"+key2,"R6_06_R1","R9_06_"+key2,"R12_06_"+key2,"R15_06_"+key2,"R18_06_"+key2,"R21_06_"+key2,"R24_06_"+key2) #warning! cycle R6_06_R2 does not exist!
        fastq_file_list= ("oligo1_R9bis_S27_"+key2+"_001","oligo1_R10_S28_"+key2+"_001") #("R9_06_R1",) #("oligo1_R9bis_S27_"+key2+"_001",)
        fastq_file_list_plot= ("0_oligo1_R9bis_S27_"+key2+"_001","1_oligo1_R10_S28_"+key2+"_001") # ("0_R9_06_R1",)  #("0_oligo1_R9bis_S27_"+key2+"_001",)
        #fastq_file_list_plot=("0_oligo1_R8_S25_"+key2+"_001","1_oligo1_R9_S26_"+key2+"_001","2_oligo1_R9bis_S27_"+key2+"_001","3_oligo1_R10_S28_"+key2+"_001") #("0_oligo1_"+key2,"1_R1_06_"+key2,"2_R2_06_"+key2,"3_R3-06_"+key2,"4_R6_06_R1","5_R9_06_"+key2,"6_R12_06_"+key2,"7_R15_06_"+key2,"8_R18_06_"+key2,"9_R21_06_"+key2,"10_R24_06_"+key2) #warning! cycle R6_06_R2 does not exist!

        if key2=="R1":
            cycles_IDs=(9,10) #(0,1,2,3,6,9,12,15,18,21,24)
        elif key2=="R2":
            cycles_IDs= (0,1,2,3,6,9,12,15,18,21,24)
            
    elif key2=="R1R2":
        #fastq_file_list_R1=("R9_06_R1",)
        #fastq_file_list_R2=("R9_06_R2",)
        #fastq_file_list_plot= ("0_R9_06_R1R2",)
        fastq_file_list_R1= ("oligo1_R9bis_S27_R1_001", "oligo1_R10_S28_R1_001")
        fastq_file_list_R2= ("oligo1_R9bis_S27_R2_001","oligo1_R10_S28_R2_001")
        fastq_file_list_plot= ("0_oligo1_R9bis_S27_R1_001R2","1_oligo1_R10_S28_R1_001R2")
        #fastq_file_list_R1=("oligo1_R1","R1_06_R1","R2_06_R1","R3-06_R1","R6_06_R1","R9_06_R1","R12_06_R1","R15_06_R1","R18_06_R1","R21_06_R1","R24_06_R1")
       # fastq_file_list_R2=("oligo1_R2","R1_06_R2","R2_06_R2","R3-06_R2","R6_06_R1","R9_06_R2","R12_06_R2","R15_06_R2","R18_06_R2","R21_06_R2","R24_06_R2") #warning! cycle 6 R2 does not exist - copy R1
        #fastq_file_list_plot=("0_oligo1_R1R2","1_R1_06_R1R2","2_R2_06_R1R2","3_R3-06_R1R2","4_R6_06_R1R2","5_R9_06_R1R2","6_R12_06_R1R2","7_R15_06_R1R2","8_R18_06_R1R2","9_R21_06_R1R2","10_R24_06_R1R2")
        cycles_IDs=(9,10) #np.concatenate([0,1,2],np.arange(3,25,3))
 
elif key1=="seriesN":
    
    target="CGGTATTGGACCCTCGCATG"

    mother_folder="../data/seriesN/" # where to take FASTQ files from
    if (key_filter==False and key_no_cut==False):
        results_folder=results_folder+"seriesN/"+key2+"_"+key3+"/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==False):
        results_folder=results_folder+"seriesN/"+key2+"_"+key3+"_key_filter/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==True):
        results_folder=results_folder+"seriesN/"+key2+"_"+key3+"_key_filter_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==False and key_no_cut==True):
        results_folder=results_folder+"seriesN/"+key2+"_"+key3+"_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    
    if key2=="R1" or key2=="R2":
        fastq_file_list=("R1N_"+key2,"R2N_"+key2,"R4N_"+key2,"R5N_"+key2,"R7N_"+key2,"R10N_"+key2,"R11N_"+key2,"R12N_"+key2,"R16N_"+key2)
        fastq_file_list_plot=("0_R1N_"+key2,"1_R2N_"+key2,"2_R4N_"+key2,"3_R5N_"+key2,"4_R7N_"+key2,"5_R10N_"+key2,"6_R11N_"+key2,"7_R12N_"+key2,"8_R16N_"+key2)
        

    elif key2=="R1R2":
        fastq_file_list_R1=("R1N_R1","R2N_R1","R4N_R1","R5N_R1","R7N_R1","R10N_R1","R11N_R1","R12N_R1","R16N_R1")
        fastq_file_list_R2=("R1N_R2","R2N_R2","R4N_R2","R5N_R2","R7N_R2","R10N_R2","R11N_R2","R12N_R2","R16N_R2")
        fastq_file_list_plot=("0_R1N_R1R2","1_R2N_R1R2","2_R4N_R1R2","3_R5N_R1R2","4_R7N_R1R2","5_R10N_R1R2","6_R11N_R1R2","7_R12N_R1R2","8_R16N_R1R2")

    cycles_IDs=(1,2,4,5,7,10,11,12,16)
 
elif key1=="oligo2":

    target="CGTATCACCAGGCAGTTGAG"

    mother_folder="../data/oligo2/" # where to take FASTQ files from
    if (key_filter==False and key_no_cut==False):
        results_folder=results_folder+"oligo2/"+key2+"_"+key3+"/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==False):
        results_folder=results_folder+"oligo2/"+key2+"_"+key3+"_key_filter/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==True):
        results_folder=results_folder+"oligo2/"+key2+"_"+key3+"_key_filter_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==False and key_no_cut==True):
        results_folder=results_folder+"oligo2/"+key2+"_"+key3+"_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    
    if key2=="R1" or key2=="R2":
        fastq_file_list=("oligo2_"+key2,"R3_oligo2_"+key2,"R6_oligo2_"+key2,"R9_oligo2_"+key2,"R12_oligo2_"+key2,"R15_oligo2_"+key2,"R18_oligo2_"+key2)
        fastq_file_list_plot=("0_oligo2_"+key2,"1_R3_oligo2_"+key2,"2_R6_oligo2_"+key2,"3_R9_oligo2_"+key2,"4_R12_oligo2_"+key2,"5_R15_oligo2_"+key2,"6_R18_oligo2_"+key2)
        cycles_IDs=np.arange(0,19,3)
        
    elif key2=="R1R2":
        fastq_file_list_R1=("oligo2_R1","R3_oligo2_R1","R6_oligo2_R1","R9_oligo2_R1","R12_oligo2_R1","R15_oligo2_R1","R18_oligo2_R1")
        fastq_file_list_R2=("oligo2_R2","R3_oligo2_R2","R6_oligo2_R2","R9_oligo2_R2","R12_oligo2_R2","R15_oligo2_R2","R18_oligo2_R2")
        fastq_file_list_plot=("0_oligo2_R1R2","1_R3_oligo2_R1R2","2_R6_oligo2_R1R2","3_R9_oligo2_R1R2","4_R12_oligo2_R1R2","5_R15_oligo2_R1R2","6_R18_oligo2_R1R2")
        cycles_IDs=np.arange(0,19,3)
        
elif key1=="negative":
        
    target="CGTATCACCAGGCAGTTGAG"
    mother_folder="../data/negative/" # where to take FASTQ files from
    
    if (key_filter==False and key_no_cut==False):
        results_folder=results_folder+"negative/"+key2+"_"+key3+"/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==False):
        results_folder=results_folder+"negative/"+key2+"_"+key3+"_key_filter/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==True and key_no_cut==True):
        results_folder=results_folder+"negative/"+key2+"_"+key3+"_key_filter_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
    elif (key_filter==False and key_no_cut==True):
        results_folder=results_folder+"negative/"+key2+"_"+key3+"_key_no_cut/min_Q_{}/".format(min_Q) # where to store results
        
    if key2=="R1" or key2=="R2":
        fastq_file_list=("R3_negative_"+key2,"R6_negative_"+key2)
        fastq_file_list_plot=("0_R3_negative_"+key2,"1_R6_negative_"+key2)
        cycles_IDs=(3,6)

    elif key2=="R1R2":
        fastq_file_list_R1=("R3_negative_R1","R6_negative_R1")
        fastq_file_list_R2=("R3_negative_R2","R6_negative_R2")
        fastq_file_list_plot=("0_R3_negative_R1R2","1_R6_negative_R1R2")
    
#################################### primers - fixed sequences - dir 5' --> 3' ####################################
# oligo 1
if key1=="oligo1" or key1=="seriesN":
    if (key2=="R1"):
        fw_primer="GGATGGGAGTGCTCTTCTTG"
        fw_full_fixed_seq="GGATGGGAGTGCTCTTCTTGAACTC" # not really used
        rev_primer="ACGATCGTCGTATCACCAG" #primer R, complemented and reverted
        rev_full_fixed_seq="ACGATCGTCGTATCACCAGGCAGTT" #full rev seq, complemented and reverted
        fixed_end_5prime="GGATGGGAGTGCTCTTCTTG"
        fixed_end_3prime="CTGGTGATACGACGATCGT"
        full_fixed_end_5prime="GGATGGGAGTGCTCTTCTTGAACTC"
        full_fixed_end_3prime="AACTGCCTGGTGATACGACGATCGT"

    else: #delete the first 3 bases of each of them
        fw_primer="TGGGAGTGCTCTTCTTG"
        fw_full_fixed_seq="TGGGAGTGCTCTTCTTGAACTC" # not really used
        rev_primer="ATCGTCGTATCACCAG" #primer R, complemented and reverted
        rev_full_fixed_seq="ATCGTCGTATCACCAGGCAGTT" #full rev seq, complemented and reverted
        fixed_end_5prime="TGGGAGTGCTCTTCTTG"
        fixed_end_3prime="CTGGTGATACGACGATCGT"

#oligo2
elif key1=="oligo2" or key1=="negative":
    if (key2=="R1"):
        fw_primer="CCCTATGCGACCCTCCG"
        fw_full_fixed_seq="CCCTATGCGACCCTCCGATGTAGAC" # not really used
        rev_primer="CGAGGATGGATCGGCAAT" #primer R, complemented and reverted - not really used
        rev_full_fixed_seq="CGAGGATGGATCGGCAATCTCAAGG"  #full rev seq, complemented and reverted
        fixed_end_5prime="CCCTATGCGACCCTCCG"
        fixed_end_3prime="ATTGCCGATCCATCCTCG"

    else: #delete the first 3 bases of each of them
        fw_primer="TATGCGACCCTCCG"
        fw_full_fixed_seq="TATGCGACCCTCCGATGTAGAC" # not really used
        rev_primer="GGATGGATCGGCAAT" #primer R, complemented and reverted - not really used
        rev_full_fixed_seq="GCCGATCCATCCTCGCTCAAGG"  #full rev seq, complemented and reverted
        fixed_end_5prime="TATGCGACCCTCCG"
        fixed_end_3prime="GCCGATCCATCCTCG"
    
fw_primer_size=len(fw_primer)
rev_primer_size=len(rev_primer)
fw_primer_size_spec=len(fw_primer) #to be used for special filter
rev_primer_size_spec=len(rev_primer)
    
if key_no_cut==True:
    fw_primer_size=0
    rev_primer_size=0
    
if key2=="R1":
    skip_bases=0
else: #set to 3 in case of R2 and R1R2 analysis - skip the first 3 bases of each line
    skip_bases=3
