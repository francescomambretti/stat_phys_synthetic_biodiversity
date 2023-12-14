#launch multiple simulations

import os
import numpy as np
import string

key1_list=("seriesN",)#("oligo1", "oligo2","negative","seriesN")
key2_list=("R1","R2","R1R2")#("R1", "R2", "R1R2")
key3_list=("all",)#("fw", "rev","all")
key_filter_list=(False,)#(True, False)
key_no_cut_list=(True,False)#(True, False)

for a in key1_list:
	for b in key2_list:
		for c in key3_list:
			for d in key_filter_list:
				for e in key_no_cut_list:
					if(d==True and e==True):
						folder="./code_{}_{}_{}_key_filter_key_no_cut".format(a,b,c)
					elif(d==True and e==False):
						folder="./code_{}_{}_{}_key_filter".format(a,b,c)
					elif(d==False and e==True):
						folder="./code_{}_{}_{}_key_no_cut".format(a,b,c)
					else:
						folder="./code_{}_{}_{}".format(a,b,c)
					if not os.path.exists(folder):
						os.mkdir(folder)
					os.system("cp code_unified/* "+folder)
					os.chdir(folder)

					s = open("input_params.py").read()
					s = s.replace("$KEY1$", format(a)).replace("$KEY2$", format(b)).replace("$KEY3$", format(c)).replace("$KEY_FILTER$", format(d)).replace("$KEY_NO_CUT$", format(e) )
					f = open("input_params.py", "w")
					f.write(s)
					f.close()

					os.system("nohup python read_fastq.py > log_LB=60_UB=150 &")
					os.chdir("../")
									
