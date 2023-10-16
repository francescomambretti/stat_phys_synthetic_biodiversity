# Written by Francesco Mambretti
# Modified on 28/10/2021: added possibility to restart after MD

import os
import string
import numpy
import sys

if (len (sys.argv) != 2):
    print("Error! Please specify whether to perform min (0), relax (1), MD (2) or restart from MD(3)")
    exit(-1)

option=int(sys.argv[1])

executable="~/oxDNA/build/bin/oxDNA"
start=0 #first simulation index
end=20 #last simulation index

if (option==0):
    # copy relevant files into main folders
    os.chdir("./initial_config")
    os.system("cp forces.txt start.conf top.top ../min/")
    os.system("cp forces.txt top.top ../relax/")
    os.system("cp top.top ../MD/")
    os.system("cp top.top ../MD_restart/")
    os.chdir("../")

    #copy model folders into single runs folders
    for i in range (start,end):
        folder='run_{}'.format(i)
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        if not os.path.exists(folder+"/min"):
            os.makedirs(folder+"/min")

        if not os.path.exists(folder+"/relax"):
            os.makedirs(folder+"/relax")

        if not os.path.exists(folder+"/MD"):
            os.makedirs(folder+"/MD")

        if not os.path.exists(folder+"/MD_restart"):
            os.makedirs(folder+"/MD_restart")        
            
        os.system("cp -r min/ relax/ MD/ MD_restart/"+folder)

        os.chdir(folder+"/min")
        os.system("nohup "+executable+" Input_min &")
        os.chdir("../../")

elif (option==1):
    for i in range (start,end):
        folder='run_{}'.format(i)

        os.chdir(folder+"/relax")
        os.system("cp ../min/stage1.dat ./") # copy final configuration of the minimization step
        os.system("nohup "+executable+" Input_relax &")
        os.chdir("../../")

elif (option==2):
    for i in range (start,end):
        folder='run_{}'.format(i)

        os.chdir(folder+"/MD")
        os.system("cp ../relax/relax.dat ./") # copy final configuration of the relax step
        os.system("nohup "+executable+" Input_MD &")
        os.chdir("../../")

elif (option==3):
    for i in range (start,end):
        folder='run_{}'.format(i)

        os.chdir(folder+"/MD_restart")        
        os.system("cp ../MD/MD.dat ./") # copy final configuration of the MD step
        os.system("nohup "+executable+" Input_MD_restart &")
        os.chdir("../../")
