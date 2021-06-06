import numpy as np
import matplotlib.pyplot as plt
import os

x=np.loadtxt("all_max_cons_overlaps.dat",unpack=True)
y=np.loadtxt("all_tot_overlaps.dat",unpack=True)

#fig = plt.figure()
  
# syntax for 3-D projection
#ax = plt.axes(projection ='3d')

#ax.plot3D(np.arange(0,10000,1), x, y, 'orangered', linestyle='none',marker='o',markersize=4)

#ax.set_xlabel("Throw")
#ax.set_ylabel("largest $a_{max}$")
#ax.set_zlabel("largest $a_{tot}$")

#plt.plot(x/y,marker='o',color='navy',linestyle='none')
num_bin=20
plt.hist(x/y,bins=num_bin,label='$a_{max}/a_{tot}$',edgecolor='black')
#plt.hist(y,bins=num_bin,label='$a_{tot}$')
#plt.show()

plt.legend()

plt.savefig("compare.png",dpi=150)
os.system("open compare.png")
