# plot the 2D-histogram of the largest maximum consecutive overlap vs the corresponding total overlap measured at the position which corresponds to the begin of the chosen largest maximum consecutive overlap subsequence

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_col
import os

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_xticks(np.arange(0,21,2))
ax.set_xticklabels(np.arange(0,21,2), fontsize=14) 
ax.set_xticks(np.arange(1,21,2), minor=True)
ax.set_yticks(np.arange(0,21,2))
ax.set_yticklabels(np.arange(0,21,2), fontsize=14)
ax.set_yticks(np.arange(1,21,2), minor=True)

x=np.loadtxt("all_max_cons_overlaps.dat",unpack=True)
y=np.loadtxt("all_tot_overlaps.dat",unpack=True)

plt.xlabel("$\\bar{a}$",fontsize=16)
plt.ylabel("$a_{tot}$",fontsize=16)

plt.hist2d(x, y, bins=(np.arange(0,21,1),np.arange(0,21,1)), cmap=plt.cm.ocean_r, density=True, norm=mpl_col.Normalize())

plt.xlim(0,21)
plt.ylim(0,21)

#plt.xticks(np.arange(0,21,2),labels=np.arange(0,21,2),ha='left',fontsize=14)
#plt.yticks(np.arange(0,21,2),labels=np.arange(0,21,2),fontsize=14)

cb=plt.colorbar()
for t in cb.ax.get_yticklabels():
	t.set_fontsize(14)

plt.clim(0,0.1)
plt.grid(which='both')

plt.savefig("compare.png",dpi=150)
os.system("open compare.png")
