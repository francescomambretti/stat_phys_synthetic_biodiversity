# Written by Francesco Mambretti
# 24/11/2021 version

import os
import string
import numpy
import sys

if (len (sys.argv) != 2):
    print("Error! Please specify whether to perform min (0), relax (1) or MD (2)")
    exit(-1)

option=int(sys.argv[1])

executable="~/oxDNA/build/bin/oxDNA"
start=0
end=10

if (option==0):
    # copy relevant files into main folders
    os.system("cp initial_config/forces.txt initial_config/start.conf initial_config/top.top min/")
    os.system("cp initial_config/forces.txt initial_config/start.conf initial_config/top.top relax/")

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

        os.system("cp -r min/ relax/ MD/ "+folder)
        os.chdir(folder+"/min")
        os.system("nohup "+executable+" Input_min &")
        os.chdir("../../")

elif (option==1):
    for i in range (start,end):
        folder='run_{}'.format(i)

        os.chdir(folder+"/relax")
        os.system("cp ../min/stage1.dat ./")
        os.system("nohup "+executable+" Input_relax &")
        os.chdir("../../")

elif (option==2):
    for i in range (start,end):
        folder='run_{}'.format(i)

        os.chdir(folder+"/MD")
        os.system("cp ../relax/relax.dat ./")
        os.system("cp ../relax/top.top ./")
        os.system("nohup "+executable+" Input_MD &")
        os.chdir("../../")
