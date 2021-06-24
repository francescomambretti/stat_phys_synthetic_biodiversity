## null_model_crosscheck
Generate a population of L-long DNA strands (predators) and subsequently check their maximum consecutive overlap with l-long target strand, looking over all the possible relative positions of the target and the predator.

Requirements: C++ compiler, C++11 std, Python3 with MPI

To compile: `make`

`null_model_constraint.cpp`  --> C++ code to generate a population of L-long DNA strands and subsequently check their maximum consecutive overlap with l-long target strand

----------------------
`null_model_constraint.py` --> Python code, sample from (null model + physical constraint) distribution; do it R times for collecting statistics over the marginal distributions. Returns for each species the mean and the error.

To run: `python3 null_model_constraint.py`

----------------------

`histo.py` --> Python code plots the RSA (relative species abundance) histogram

To run: `python3 histo.py`

----------------------

`plot_full_hist_and_cumul.py` --> Python code, load and plot marginal probability distributions for species abundance. Plots also cumulative probabilities for having >= a given species index.

To run: `python3 plot_full_hist_and_cumul.py (0/1)`

----------------------

