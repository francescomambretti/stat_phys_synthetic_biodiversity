`find_MCO` computes the Maximum Consecutive Overlap (MCO) between two sequences of any length

C++-11 required (to use vectors)

To compile: `g++ -O3 -o find_MCO.x find_MCO.cpp --std=c++11`

******************************

The following scripts require Python 3.

******************************

`mixed_hb_pair.py` counts the total number of mixed (TMO) hydrogen bonds (HB) and the MCO between 2 strands: either 1 predator + 1 resource, or 1+1 predators or 1+1 resources as a function of time

To be used like: `python mixed_hb_pair.py hb_list.dat N M output_file_name` where `N` and `M` are the numbers of predators and resources, respectively
The program outputs a file called `output_file_name` with the computed quantities

******************************

`mixed_hb_pair_evo.py` counts the total number of mixed hydrogen bonds (HB) and the MCO between an arbitrary group of strands as a function of time. Also self interactions are considered

To be used like: `python mixed_hb_pair_evo.py`, where all the input parameters are stored in a file like the one provided (`input_params.py`).
The program outputs a file called `output_file_name` (set into `input_params.py`) with the computed quantities

******************************

`loop_run.py` is used to execute several independent simulations. 

Requires, as input parameter, whether to perform min (0), relax (1), MD (2) or restart from MD(3).
According to the input option, the script executes `end`-`start` (specified inside the script) concurrent simulations, each one with a different random seed. 
The code expects that the `prep_data` script has already been run.
