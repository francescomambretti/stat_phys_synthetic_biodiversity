`find_MCO` computes the Maximum Consecutive Overlap (MCO) between two sequences of any length

C++-11 required (to use vectors)

To compile: `g++ -O3 -o find_MCO.x find_MCO.cpp --std=c++11`

******************************

`mixed_hb_pair.py` counts the total number of hydrogen bonds (HB) and the MCO between 2 strands: either 1 predator + 1 resource, or 1+1 predators or 1+1 resources as a function of time

To be used like: `python mixed_hb_pair.py hb_list.dat N M output_file_name` where `N` and `M` are the numbers of predators and resources, respectively
The program outputs a file called `output_file_name` with the computed quantities

