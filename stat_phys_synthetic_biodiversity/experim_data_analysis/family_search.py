# Written by Francesco Mambretti, 14/1/2022
# Track the time evolution of substrand abundances in experimental data
# 22/01/2022 version - future versions will receive input parameters from external file

import numpy as np
import os
import string
import matplotlib.pyplot as plt

#main_path="/mnt/FRANCESCO/experim_data/results/r9_06_r1-fastq-gz/R2_fw_rev/min_Q_10/"
main_path="/mnt/FRANCESCO/experim_data/results/r9_06_r1-fastq-gz/R1_forward/min_Q_10/"

fastq_file_list_R1=("0_oligo1_R1","1_R3-06_R1","2_R6_06_R1","3_R9_06_R1","4_R12_06_R1","5_R15_06_R1","6_R18_06_R1","7_R21_06_R1")
#fastq_file_list_R2=("0_oligo1_R2","1_R3-06_R2","2_R9_06_R2","3_R12_06_R2","4_R15_06_R2","5_R18_06_R2","6_R21_06_R2")

#family=("TCGAGAGGGCTTCGGATTAGGTGCTGGTGATACGACGATCGTACAGGTGGTGTAAC","ACTCAGGCAGGCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGCC","AACCGCCACACACTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGC","ACTGCCACATCGACTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGG","ACTCGGCCAGCCACTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGG","TGCCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGCCTCGAACTG","TGCCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGCCGTGCTGGT","TGTGCCGAACTGCTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGGGC")

#family=("CTGGTGATACGACGATCGT","GGATGGGAGTGCTCTTCTTG")
#family=("CTGGTGATACGACGATCGTACAGGTGGTGTAACATGCGAGGG",)

###these are among top-100 of R1 forward at cycle 21, and do not contain neither primer F nor primer R complemented and reverted. 
family=("ACCACGCCAAGACTTCTGACCATGCGAGGGTCCGTTACGGATGTCGATCGGGTAAC","CGTGCACTGAAAGGACGCGTCGTGGTAGGGGGACGTCATGCGAGGGTTCCTAAGTG","TGCAGCGAGGGCTCAATGCGAGGGTGTCCGGCCTGGAAACTGCCTGGTGATACGAC","ACTCGGCAGGTAGGCGGTCCCTTTGACATGCGAGGGTCCACGTCGGTGTGTCTTA","ACTCGCCGCTGGTCGTGGTTGCGAGGGTCCATCCCTAGTTCAGAGCGTTGGTGCA","CGCGTCGACCATTGCGAGGGACAATACGGGAGTGAACTGCCTGGTGATACGACGAT","ACCGCGCATGGCCTTGACGTAGCGACGGTGTCTGTGGCGCGATGGAGGGCGGTAAC","CCGTGCCAGGCCATGCGAGGGTTAACTGCCTGGTGATACGACGACTGGTGATACG","AACTCGCAGCCAGGGAGCAAGTCTGTTCATGGGAGGGTCCATCGTGCGCGTCGGC","ACTCGAACGACGCACGTGTTGGTGCAGAAGGCGTGCGCCAGTGTGCGAGGGTCCAT","CGTGCGGGGTCCAGGATGCGAGGGTCGAGGCGGGTAAGCTGGGGGCGGGCTAACTG","TGCGGCGAGGGGTTCGGCGGCGAGTGGTCCGTGTCGTCGGATTGGCCTAACTGCCT","TGCGTGCGGGCGAGGGGGCAATCCGGTCTGGGGTCCAACTGCCTGGTGATACGACG","CCACAGGCGAGGGACCGGTGGCGTGGGACAGGGGTGAACGCCTGGTGATACGACG","CGACAGCGGGGGGACAGTCACTGCGGGCATGCAGGGTGGTGCGGGCGTAACTGCCT","CCGTGCCAGGCCATGCGAGGGTTAACTGCCTGGTGATACGACGATCTGCGAGGGGG","ACGGCCCGTGGTGGCGAATGCTGTTGCGGCGCATGATGCGTGGGTGCTTGCTAACT","ACCTGCGAGGGTCCGTGGTTACGATTAGCGGGAGAACAGTGGCATGTCGGCGCAAC","ACTCGGGCGGCGCAGCGTTGTGCGTAGGGTCCTGTCCTAGTTCTGGTGATACGACG")

iterator=-1 #loop over strands
titles=("5","14","20","26","36","39","41","50","53","56","57","58","61","63","64","69","90","94","98") #these numbers are indexes in top-100
#titles=("primer_R","primer_F")
#titles=("93","75","8","7","31","2","33","6")

#colors=("red","navy")
#colors=("black","brown","red","orange","gold","forestgreen","turquoise","gray")
colors=("blue","cyan","red","brown","lime",'#890000','#890000','#5c0000','#2a6b28','#0b4c07','#003206','#4f5a90','#374798','#30375a','#fff4b1','#ffed86','#ffdb00',"pink","purple","magenta") #may be improved

for strand in family:
	iterator+=1
	iterator2=-3 #loop over cycles
	count_arr=np.empty(0)
	print(strand)
	
	for fastq_file_R1 in fastq_file_list_R1:
		iterator2+=3
		input_file=main_path+fastq_file_R1+"/all.txt"
		input_file=open(main_path+fastq_file_R1+"/all.txt", 'r')
#		output_file=open(main_path+fastq_file_R1+"/{}_contain_special_seq.txt".format(iterator2), 'w')
		Lines = input_file.readlines()
		 
		count = 0
		for line in Lines:
			if (strand in line.strip()):
				count += 1
#				output_file.write(line)
		count_arr=np.append(count_arr,count)
		print(fastq_file_R1,count)
#		if(iterator2==0):
#			plt.plot(iterator2,count,marker='s',linestyle='dotted',label=titles[iterator],color=colors[iterator])
#		else:
#			plt.plot(iterator2,count,marker='s',linestyle='dotted',color=colors[iterator])

	plt.plot(np.arange(0,22,3),count_arr,marker='s',linestyle='dotted',label=titles[iterator],color=colors[iterator])
	
plt.xticks(np.arange(0,22,3))
plt.xlabel("Cycle",fontsize=14)
plt.ylabel("Counts",fontsize=14)	
plt.yscale('log')
plt.legend(fontsize=10,ncol=2,loc='upper left')
plt.savefig("not_aliens_R1_fw.png",dpi=300) #to be changed, if you are plotting anything different
plt.close()
