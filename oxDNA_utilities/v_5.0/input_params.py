runs=1
hb_list_file="hb_list.dat"
n_pred=2
n_res=1
caps_per_pred=2 #set to 0 to analyze run without blocking caps
output_folder="results"
output_file_name="hb_counter"
members="p5-p13-res "
pred_indexes=[0,3]
L_list=[101,101] #list of predators lengths, caps included
caps_indexes=[]
for i in pred_indexes:
    caps_indexes.append(i+1)
    caps_indexes.append(i+2)
res_indexes=[6,]
