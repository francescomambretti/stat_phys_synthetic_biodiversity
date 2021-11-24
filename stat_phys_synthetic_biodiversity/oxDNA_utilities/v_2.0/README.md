******************************

All these files were created by Francesco Mambretti, to set up and analyze oxDNA (https://dna.physics.ox.ac.uk/index.php/Main_Page) simulations.

******************************

`find_MCO` computes the Maximum Consecutive Overlap (MCO) between two sequences of any length

C++-11 required (to use vectors)

To compile: `g++ -O3 -o find_MCO.x find_MCO.cpp --std=c++11`

******************************

The following scripts require Python 3.

Note: `mixed_hb_pair.py` is not used anymore.

******************************

`mixed_hb_pair_evo.py` counts the total number of mixed hydrogen bonds (HB) and the MCO between an arbitrary group of strands (how many predaors and how many resources you like) as a function of time, in all the simulations of a given set of simulations. Also self interactions are considered.

To be used like: `python mixed_hb_pair_evo.py`, where all the input parameters are stored in a file like the one provided (`input_params.py`).
The program outputs a file called `output_file_name` (set into `input_params.py`) with the computed quantities.
The code also generates the plots of:
- the histograms of Maximum Consecutive Overlap and Total Mixed Overlap between any pair of strands in the system
- the average over _n_ independent simulations of the MCO and TMO as a function of time, with errorbars

******************************

`loop_run.py` is used to execute several independent simulations. 

Requires, as input parameter, whether to perform min (0), relax (1), MD (2) or restart from MD(3).
According to the input option, the script executes `end`-`start` (specified inside the script) concurrent simulations, each one with a different random seed. 
The code expects that the `prep_data` script has already been run.


******************************

`initial_config` is a folder that contains file which allow to set up the simulations. 
`prep_data` is a Bash script which calls some Python scripts, with the target to generate external forces files, and the initial configuration and the topology files.

******************************

`input_files` is a folder that contains input files for typical oxDNA simulations:
- minimization (`min`)
- relax (both with CPU and GPU-parallelized flavours)
- Molecular Dynamics (`MD`) (both with CPU and GPU-parallelized flavours)
