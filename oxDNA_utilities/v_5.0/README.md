******************************

All these files were created by Francesco Mambretti, to set up and analyze oxDNA (https://dna.physics.ox.ac.uk/index.php/Main_Page) simulations.

******************************

`find_MCO` computes the Maximum Consecutive Overlap (MCO) between two sequences of any length

C++-11 required (to use vectors)

To compile: `g++ -O3 -o find_MCO.x find_MCO.cpp --std=c++11`

N.B.: fixed bugs of the previous versions, affecting MCO and TMO computation.

******************************

Where not differently indicated, the following scripts require Python 3:

******************************

`mixed_hb_pair_evo.py` counts the total number of mixed hydrogen bonds (HB) and the MCO between an arbitrary group of strands (how many predators and how many resources you like) as a function of time, in all the simulations of a given set of simulations. Also self interactions are considered.
Useful functions are implemented in `functions_mixed_hb.py`, and input parameters must be set in `input_params.py`.

To be used like: `python mixed_hb_pair_evo.py`, where all the input parameters are stored in a file like the one provided (`input_params.py`).
The program outputs a file called `output_file_name` (set into `input_params.py`) with the computed quantities.
The code also generates the plots of:
- the histograms of Maximum Consecutive Overlap and Total Mixed Overlap between any pair of strands in the system
- the average over _n_ independent simulations of the MCO and TMO as a function of time, with errorbars

From version 3.0, the code is also capable of working with caps and fixed sequences at the ends of the predator.
From version 4.0, the code manages an arbitrary number of predators with/without caps (the two conditions can not be merged, at the moment) and an arbitrary number of resources.
From version 5.0, predator strands can have different lengths. Moreover, strands indexes and lengths are provided in the input file.

******************************

`loop_run.py` is used to execute several independent simulations. 

Requires, as input parameter, whether to perform min (0), relax (1), MD (2).
According to the input option, the script executes n=`end`-`start` (specified inside the script) concurrent simulations, each one with a different random seed. 
The code requires the `prep_data.py` script to have already been executed.

******************************

`initial_config` is a folder that contains file which allow to set up the simulations. 
`prep_data.py` is a Python script which calls other Python scripts, with the target to generate external forces files, and the initial configuration and the topology files.
`prep_data.py` substitutes the old `prep_data` Bash script.

Python2 is still required for some of the Python scripts generating forces.
See the dedicate README.md inside the folder.

******************************

`input_files` is a folder that contains input files for typical oxDNA simulations:
- minimization (`min`)
- relax (both with CPU and GPU-parallelized flavours)
- Molecular Dynamics (`MD`) (both with CPU and GPU-parallelized flavours)

******************************

`generate-sa-caps.py` is a customized version of the `generate-sa.py` script provided with oxDNA, in the `UTILS` directory. Allows to generate ssDNAs with caps at their ends, randomly displaced in a simulation box. Caps length is specified as input parameter (after box side and sequences filename)

`generate-ft-caps.py` does the same thing but requires a `start_pos.txt` argument, since it works with user-defined positions for the strands (no randomly placing into the box)

N.B.: these scripts work with Python 2!

******************************

Added all the `persistence_length` folder for computation and plot of tangent vectors correlation functions. Customized version, adapted also for strands with caps (it is possible to compute correlation function only on a fraction of the total strand). See dedicated README.md in the folder.
