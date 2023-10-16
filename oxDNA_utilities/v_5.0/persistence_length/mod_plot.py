#Written by Francesco Mambretti
#29/04/2022

import numpy as np
import global_var as g
import matplotlib.pyplot as plt

########################################
def plot_blk_glb_ave(i_blk,filename):
    y1=g.ave_corr
    y2=g.glob_ave_corr/(i_blk+1)
    err_y2=g.error_corr
    x=np.arange(0,g.N-1)
    
    plt.clf()
    plt.plot(x,y1,marker='o',linestyle='None',color='orangered',label='blk {}'.format(i_blk))
    plt.errorbar(x,y2,yerr=err_y2,marker='o',capsize=5,linestyle='None',color='gold',label='glob ave')
    
    plt.xlabel("$n$",fontsize=14)
    plt.ylabel("$\\frac{<t_n \cdot t_0>}{<t_0>^2}$",fontsize=14)
    
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename,dpi=250)
