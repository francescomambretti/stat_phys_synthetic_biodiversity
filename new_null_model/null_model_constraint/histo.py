import numpy as np
import matplotlib.pyplot as plt
import os

data=np.loadtxt('all_max_cons_overlaps.dat', unpack=True)
plt.hist(x=data, color='green',alpha=0.8, edgecolor='black',align='left',bins=np.arange(0,21,1))
plt.xticks(np.arange(0,21,1))
plt.xlabel("overlap")
plt.ylabel("Count (overlap)")

plt.tight_layout()
plt.savefig("histogram.png",dpi=150)
os.system("open histogram.png")
