## null_model_fixed_pos

Generate a population of L-long DNA strands and subsequently check their (largest) maximum consecutive overlap with target strand. Compute also the largest total overlap with the target strand, aiming to plot them together. "Largest" means that all possible relative positions between the target and the predator are checked (null model + overlap). Built on the basis of null_model_constraint (see code in other folder of this repository).

Requirements: `C++ compiler, C++11 std, Python3`

To compile: `make`

----------------------
`compare.py` --> Python code, do various comparisons/correlations between max consecutive and tot overlaps. Here, the majority of lines are commented, can be uncommented according to the desired functionality

To run: `python3 compare.py`
