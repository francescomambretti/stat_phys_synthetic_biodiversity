## README - Genetic Algorithm v 2.0

************************************************************************************************************************************************

### Requirements

- MPI-C++ compilers (tested with MPICH and OpenMPI flavors)
- Armadillo libraries, http://arma.sourceforge.net/ (tested with 9.900.5 and 10.6.2 versions)
- OpenMP libraries
- optional (but useful): Lapack & BLAS

************************************************************************************************************************************************

### Compile and execute
To compile: `make`
To remove object files and executable: `make clean`
To execute: `mpiexec -np PROCS ./evolution        input file name     random generator seed       mutations on/off (1/0)` where `PROCS` is the chosen number of MPI parallel processes

************************************************************************************************************************************************

### Structure of the code
- `evolution.cpp` is the main file
- functions implementations are mostly gathered inside `functions.cpp`, while `functions.h` is the related header
- `global.h` contains the declaration of global variables with the `extern` keyword
- the random number class is declared and defined in the `random.h` and `random.cpp` files
- `seed.in` contains the initial common seed, `Primes` contains couples of numbers used to distinguish the random numbers sequence generated for each rank, `input.dat` the input information and the `makefile` is self-explaining
- Jupyter Notebooks and Python scripts for analysis are provided in the current version (see below)

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

The code outputs a file called `histogram*.dat` for each generation `*`, which contains the list of the fitnesses of each of the survived predators, before the amplification process.
The program also outputs some files named `pred_histF_cycle_*.dat` where `F` represents the fitness and `*` the cycle. This file contains the list of the `L`-long sequences having a fitness `F` after cycle `*`. This is a way to spot the emergence of dominating sequences during an evolution process.

************************************************************************************************************************************************

### Analysis

Several scripts are provided here for analysis, in the `scripts` folder.

- `compare_GA_experim.py`: this script compares the RSA histogram of an experimental cycle with all the corresponding histograms of all the simulated GA cycles for all the `alpha` values. It computes the RMS distance (Euclidean) between the simulated RSA histogram and the experimental one. It also performs the Kolmogorov-Smirnov test to confront the underlying distributions. Both quantities aim to estimate the best combination of (alpha,cycle index), so to minimize a distance in a given metrics between the GA histogram and the experimental one.
Options for plot (one plot for each (alpha,cycle index)) are commented by default.

- `RSA_alpha.py`: script which analyzes the effect of the alpha parameter on the RSA histogram after the GA simulation. In particular, plot the RSA histogram for a single chosen cycle of the algorithm and for all the alpha values specified.

- `RSA_histo_evo.py`: plots some intermediate steps of the evolution of the RSA histogram during a GA simulation. Every `g` generations, the RSA histogram of all the independent simulations performed is averaged and errorbars are computed. As outcome, chosen a GA simulation with a given selection criterion, and for fixed `alpha`, the code outputs two plots with all the histograms superimposed, in linear and log scale.

-`compare_ga_null_exp.py`: compare Genetic Algorithm results (for all the three selection criteria), the null model and the experimental data. Compare the RSA histogram of the desired set of GA cycles, the experimental histogram (at a given cycle) and the null models with/without threshold. Data are averaged over all the independent runs performed.
