## README - Genetic Algorithm v 5.0

************************************************************************************************************************************************

### Requirements

- MPI-C++ compilers (tested with MPICH and OpenMPI flavors)
- Armadillo libraries, http://arma.sourceforge.net/ (tested with 9.900.5 and 10.6.2 versions)
- OpenMP libraries
- optional (but useful): Lapack & BLAS
- Python3

************************************************************************************************************************************************

### Compile and execute
To compile: `make`
To remove object files and executable: `make clean`
To execute: `mpiexec -np PROCS ./evolution        input file name     random generator seed       mutations on/off (1/0)    folder_for_output` where `PROCS` is the chosen number of MPI parallel processes

************************************************************************************************************************************************

### Structure of the code
- `evolution.cpp` is the main file
- functions implementations are mostly gathered inside `functions.cpp`, while `functions.h` is the related header
- `global.h` contains the declaration of global variables with the `extern` keyword
- the random number class is declared and defined in the `random.h` and `random.cpp` files
- `seed.in` contains the initial common seed, `Primes` contains couples of numbers used to distinguish the random numbers sequence generated for each rank, `input.dat` the input information and the `makefile` is self-explaining
-  Python scripts for analysis are provided in the `scripts` folder (see below)``

************************************************************************************************************************************************

### What does this code do?

The code simulates a population of `N_pred` predator DNA strands (`L`-long arrays, where each element of the array could take one of 2/4 values), which interact with a population of `N_res` DNA strands called resources, each one of them having the same length `l`. Note that `l < L` (but this is not mandatory for the code to correctly run).
Note also that it must be `N_pred > N_res`. 
For each generation (the number can be set in the input file), `N_res` predators are selected according to a given criterion. They are subsequently amplified to become again `N_pred`-sized.
Each parallel rank executes its own evolution process, and the results are saved in a dedicated folder, called `seedS`, where `S` is a unique index referring to each rank.

Currently, three selection criteria have been inserted: 
- selection based on the fitness of each predator (defined as the maximum consecutive overlaps between the predator and the resource), - option 0
- selection based on the fitness, with an asymptotic fitness attributed after a threshold, - option 2
- purely random selection, regardless of the fitness of the predator - option 1
The user has to choose one of them, and also has to set the value of the `alpha` parameter, which regulates the amount of purely random selection moves at each cycle (`alpha=0` means no random selection, `alpha=1` means only random selection).

The code outputs a file called `histogram*.dat` for each generation `*`, which contains the histogram of the fitnesses of each of the survived predators, before the amplification process. (NEW FROM V_5.0).
The code also outputs, for each cycle, the full list of the sequences with their fitness (Maximum Consecutive Overlap with the resource, ref. https://www.mdpi.com/1099-4300/24/4/458). The same list is also split in separate files for each MCO value.

************************************************************************************************************************************************

### Analysis

Several scripts are provided here for analysis, in the `scripts` folder. These are much different from the previous versions

- `analysis.py`: main script. The code processes all the simulations and analyses them, computing:
 0) the mean histogram p(MCO)
 1) the associated Shannon Entropy
 2) the number of unique sequences
 3) the Shannon Entropy associated to the Relative Species Abundance, where species=sequence
 4) the average omega of the population
 5) the zip ratio between the original and the zipped file size for the whole list
 6) the evolution of most abundant strands at each cycle
 7) the fraction of population covered by top-n individuals

The functions are implemented in `module_functions.py`.

To execute: `python analysis.py`, after setting all the relevant parameters in `input_params.py`, which stores parameters & global options.

- `main_plot.py` performs all the plots of the quantities listed above; some are plotted for each single seed, others are only plotted after averaging over parallel simulations. The script is also called at the end of  `analysis.py`. The associated methods are implemented in `module_plots.py`

- The `regime` folder contains three scripts (main, input options and functions implementation) to compute and plot the mean trend of the abundance within a population of sequences having a given fitness. The groups are divided in small-, medium- and high- fitness.

The `old` folder contains the script used in the previous versions: `compare_GA_experim.py`, `RSA_alpha.py`, `RSA_histo_evo.py`, `compare_ga_null_exp.py`.`
