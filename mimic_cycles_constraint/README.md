## mimic_cycles constraint
Generate a population of L-long DNA strands (predators) and subsequently check their maximum consecutive overlap with l-long target strand, for fixed relative position between the predator and the target.


Requirements: Python3

1) `python mimic_cycles_constraint.py`
Mimick a null-model driven evolution, by repeatedly sampling from the null model distribution with the physical constraint.
Periodically prints output data.

----------------------

2) `python plot.py` 

Plot the histograms of the chosen cycles for a single run

----------------------

3) `python ave.py` 

Compute the average histogram corresponding to the average among many different runs.

----------------------

4) `python plot_mimic_vs_experim.py` 

Plot the ratio between such an histogram and the corresponding experimental data.
