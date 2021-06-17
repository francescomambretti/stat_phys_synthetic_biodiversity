# plot the 2D-histogram of the largest maximum consecutive overlap vs the corresponding total overlap measured at the position which corresponds to the begin of the chosen largest maximum consecutive overlap subsequence

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_col
import os

x=np.loadtxt("all_max_cons_overlaps.dat",unpack=True)
y=np.loadtxt("all_tot_overlaps.dat",unpack=True)

plt.xlabel("$\\bar{a}$")
plt.ylabel("$a_{tot}$")

plt.hist2d(x, y, bins=(np.arange(0,21,1),np.arange(0,21,1)), cmap=plt.cm.ocean_r, density=True, norm=mpl_col.Normalize())
plt.xlim(0,21)
plt.ylim(0,21)

plt.xticks(np.arange(0,21,1),labels=np.arange(0,21,1),ha='left')
plt.yticks(np.arange(0,21,1),labels=np.arange(0,21,1))

plt.colorbar()
plt.clim(0,0.1)
plt.grid()

plt.savefig("compare.png",dpi=150)
os.system("open compare.png")
