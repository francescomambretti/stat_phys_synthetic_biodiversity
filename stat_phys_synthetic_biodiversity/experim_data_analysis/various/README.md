`python count_bases.py filename`

Counts the number of A,C,G,T in a list of DNA sequences

-------------------------------------

C++ code for analyzing the sequences:

`./experim_overlaps.x target_file sequences_file_path output_folder_path`

to compile: simply `make` (C++-17 is used, but C++-11 is enough)

The code is described at the beginning of the .cpp file

-------------------------------------

`histo.py`: does the histogram of the relative species abundance. It has to be internally modified to change the input files

-------------------------------------

`compare.py`: 2D histogram of maximum and total consecutive overlaps

-------------------------------------

`family_search.py`: tracks the abundance of strands containing (or not) a (set of) given subsequence(s) along the experimental cycles (ususally backwards)


-------------------------------------

`nullmodel_vs_exp_data.ipynb`: compare theoretical (null) model for the distribution of the Maximum Consecutive Overlap $\omega$ with experimental data

-------------------------------------

`decay_rate_and_MCO_LTO.ipynb`: MCO distribution trend, time dependence of the height of each single bin with and without dominant individuals. 
Second part: 2D histogram and correlation between MCO and LTO of experimental data
