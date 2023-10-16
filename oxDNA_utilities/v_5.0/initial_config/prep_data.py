# Original bash version by Nicolo' Pedrani
# Edited by Francesco Mambretti, 15/04/2022 Python version
# N.B.: this is the first Python version of the previous bash code
# caps traps are no more available, they cause errors
# harmonic traps to confine strands are replaced by repulsive spheres (apparently, they work better)

import os 

print(" ---------------------------------------------")
print("------- preparing data for simulation -------")
print(" ---------------------------------------------")

####################
#input params
####################

#start_pos.txt ready - if necessary
FORCE="~/oxDNA_UTILS/forces"

box=30 
box2=box/2 
box4=box/4
N_res=1

force_opt="sphere_and_mutual" #sphere, sphere_and_mutual, harmonic_and_mutual

print(box,box2,box4,N_res)

####################
#generate config
####################

print(" °°°°°°°°° start configurations °°°°°°°°°°°°°°")

#python2.7 ~/oxDNA/UTILS/generate-ft-caps.py $box seq.dat 25
os.system("python2.7 ~/oxDNA/UTILS/generate-sa-caps.py "+ str(box) +" seq.dat 25")

print("°°°°°°°° end start configurations °°°°°°°°°°°")

os.system("mv generated.dat start.conf")
os.system("mv generated.top top.top")

####################
#forces
####################

print(" °°°°°°°°°° forces files °°°°°°°°°°°°°°°°°°°°°")
#traps to make strands closer one to each other - or to make strands closer to a given point

if (force_opt=="sphere"): #suitable when using generate-sa-caps.py
    os.system("python rep_sphere.py "+str(box2)+" "+str(box4)+" "+str(N_res)) #generate sphere.txt
    os.system("python "+FORCE+"/generate_nico.py sphere.txt rep_sph")
    os.system("cat sphere.txt >forces.txt")

elif (force_opt=="sphere_and_mutual"): #suitable when using generate-ft-caps.py
    os.system("python rep_sphere.py "+str(box2)+" "+str(box4)+" "+str(N_res)) #generate sphere.txt
    os.system("python mut_trap.py") #generate trapmut.txt
    os.system("python "+FORCE+"/generate_nico.py sphere.txt rep_sph")
    os.system("python "+FORCE+"/generate_nico.py trapmut.txt mut_trap") # generate mut_traps.txt
    os.system("cat sphere.txt trapmut.txt > forces.txt")

elif (force_opt=="harmonic_and_mutual"): #N.B.: harmonic forces, apparently, raise issues when used with (partly) double stranded DNA
    os.system("python harm_trap.py "+str(box2)) #generate trapharm.txt
    os.system("python mut_trap.py") #generate trapmut.txt
    os.system("python "+FORCE+"/generate_nico.py trapharm.txt harm_trap") # generate harmonic_traps.txt
    os.system("python "+FORCE+"/generate_nico.py trapmut.txt mut_trap") # generate mut_traps.txt
    os.system("cat trapharm.txt trapmut.txt > forces.txt")


print("°°°°°°°°° end forces files °°°°°°°°°°°°°°°°°°")
print("---------------------------------------------")
print("----------------- DONE ----------------------")
print("---------------------------------------------")
