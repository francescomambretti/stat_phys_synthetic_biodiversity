## null_model_fixed_pos
Generate a population of L-long DNA strands (predators) and subsequently check their maximum consecutive overlap with l-long target strand, for fixed relative position between the predator and the target.

Requirements: C++ compiler, C++11 std, Python3

To compile: make

----------------------
nullmodel.py --> Python code, sample from null model distribution; do it R times for doing statistics over the marginal distributions 

To run: python3 nullmodel.py

----------------------

plot_full_hist.py --> Python code, load and plot marginal probability distributions for species abundance

To run: python3 plot_full_hist.py


----------------------

Observation: extracting mean and std. dev. of the entries of full_hist, from full_hist.txt file, simply requires a few Python lines with numpy. A simple example follows, but it's nothing sophisticated and it can be edited upon need.

```
import numpy as np

data=np.loadtxt("full_hist.txt",unpack=True,dtype=np.float64)
data=data.T.  #transpose

for i in range (0,21):
        print(i,np.mean(data[i],),np.std(data[i]),np.mean(data[i])+np.std(data[i]))
```
