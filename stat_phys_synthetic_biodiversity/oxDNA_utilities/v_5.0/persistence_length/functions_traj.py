#Written by Francesco Mambretti
#29/04/2022

import global_var as g

########################################
def find_breaks(steps,blocks): #useful to identify those steps which mark the passage from one block to the other
    list_breaks=[0] * blocks
    for i in range(0,blocks):
        list_breaks[i]=g.steps_per_block*i
        
    return list_breaks

########################################
def read_traj(split_filename,list_breaks,start_nucl,end_nucl):
    i_blk=0
    next=i_blk+1
    count=0
    
    file=open(g.traj_file,"r") #open traj_file
    outfile=open(split_filename+"_{}.txt".format(i_blk),"w")
    
    for line in file.readlines():
        if (line.split()[0]=='t'): #new timestep
            count=0
            if (int(line.split()[2])==list_breaks[next]): #check whether it is time to stop the current block
                i_blk+=1
                next=i_blk+1
                if (next==g.blocks):
                    next=0
                outfile.close()
                outfile=open(split_filename+"_{}.txt".format(i_blk),"w")
            outfile.write(line)
        elif (line.split()[0]!='E' and line.split()[0]!='b'):
            count+=1 #keep track of nucleotide index
            #for the selected nucleotides
            if (count>= start_nucl and count <=end_nucl):
                outfile.write(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+"\n") #select only relevant columns
            
    return
