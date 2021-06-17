# created by Francesco Mambretti on 15/06/2021
# plot the histogram of the relative species abundance for three cycles - can be easily adapted to arbitrary number of cycles
#


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_col
import os
import sys

file1="./R3-03_ready/all_max_cons_overlaps.dat"
file2="./R6-03_ready/all_max_cons_overlaps.dat"
file3="./R9-03_ready/all_max_cons_overlaps.dat"

data1=np.loadtxt(file1, unpack=True)
data2=np.loadtxt(file2, unpack=True)
data3=np.loadtxt(file3, unpack=True)
weights1 = np.ones_like(data1)/len(data1)
weights2 = np.ones_like(data2)/len(data2)
weights3 = np.ones_like(data3)/len(data3)

histo1,bins1=np.histogram(data1, bins=np.arange(0,22,1))
histo2,bins2=np.histogram(data2, bins=np.arange(0,21,1))
histo3,bins3=np.histogram(data3, bins=np.arange(0,21,1))

np.savetxt('./R3-03_ready/histogram.txt',histo1,fmt='%d')
np.savetxt('./R6-03_ready/histogram.txt',histo2,fmt='%d')
np.savetxt('./R9-03_ready/histogram.txt',histo3,fmt='%d')

plt.hist(x=data1, color='blue',alpha=1.0, edgecolor='black',rwidth=0.95,align='left',bins=np.arange(0,21,1),weights=weights1,label='cycle 3')
plt.hist(x=data2, color='red',alpha=1.0, edgecolor='black',rwidth=0.65,align='left',bins=np.arange(0,21,1),weights=weights2,label='cycle 6')
plt.hist(x=data3, color='gold',alpha=0.6, edgecolor='black',rwidth=0.4,align='left',bins=np.arange(0,21,1),weights=weights3,label='cycle 9')
plt.xticks(np.arange(0,21,1))
plt.xlabel("largest $a_{max}$ - species index")
plt.ylabel("relative abundance")

plt.yscale('log')

plt.legend()

plt.tight_layout()
plt.savefig("histogram-03.png",dpi=150)
os.system("open histogram-03.png")
