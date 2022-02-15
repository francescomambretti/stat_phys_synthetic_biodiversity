# Written by Francesco Mambretti, 14/01/2022
# Track the time evolution of substrand abundances in experimental data
# 15/02/2022 version - future versions will receive input parameters from external file

import numpy as np
import os
import string
import matplotlib.pyplot as plt

main_path="/mnt/FRANCESCO/experim_data/results_LB=0/oligo1_proto06/R1_fw_key_no_cut/min_Q_10/"

fastq_file_list_R1=("0_oligo1_R1","1_R3-06_R1","2_R6_06_R1","3_R9_06_R1","4_R12_06_R1","5_R15_06_R1","6_R18_06_R1","7_R21_06_R1","8_R24_06_R1")

family=("TCGAGAGGGCTTCGGATTAGGTGCTGGTGATACGACGATCGTACAGGTGGTGTAAC","ACTCAGGCAGGCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGCC","AACCGCCACACACTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGC","ACTGCCACATCGACTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGG","ACTCGGCCAGCCACTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGG","TGCCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGCCTCGAACTG","TGCCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGCCGTGCTGGT","TGTGCCGAACTGCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGC")

iterator=-1 #loop over strands
titles=("93","75","8","7","31","2","33","6")

colors=("black","brown","red","orange","gold","forestgreen","turquoise","gray")

for strand in family:
	iterator+=1
	iterator2=-3 #loop over cycles
	count_arr=np.empty(0)
	print(strand)
	
	for fastq_file_R1 in fastq_file_list_R1:
		iterator2+=3
		input_file=main_path+fastq_file_R1+"/all.txt"
		input_file=open(main_path+fastq_file_R1+"/all.txt", 'r')
		Lines = input_file.readlines()
		 
		count = 0
		for line in Lines:
			if (strand[:-1] in line.strip()):
				count += 1
		count_arr=np.append(count_arr,count)
		print(fastq_file_R1,count)

	plt.plot(np.arange(0,25,3),count_arr,marker='s',linestyle='dotted',label=titles[iterator],color=colors[iterator])
	
plt.xticks(np.arange(0,25,3))
plt.xlabel("Cycle",fontsize=14)
plt.ylabel("Counts",fontsize=14)	
plt.yscale('log')
plt.legend(fontsize=10,ncol=2,loc='upper left')
plt.savefig("R1_fw_LB=0.png",dpi=300) #to be changed, if you are plotting anything different
plt.close()
