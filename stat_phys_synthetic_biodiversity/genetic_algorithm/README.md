## README - Genetic Algorithm v 1.0

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
- Jupyter Notebooks for analysis are not provided in the current version (still under development)

************************************************************************************************************************************************

### What does this code do?

The code simulates a population of `N_pred` predator DNA strands (`L`-long arrays, where each element of the array could take one of 2/4 values), which interact with a population of `N_res` DNA strands called resources, each one of them having the same length `l`. Note that `l < L` (but this is not mandatory for the code to correctly run).
Note also that it must be `N_pred > N_res`. 
For each generation (the number can be set in the input file), `N_res` predators are selected according to a given criterion. They are subsequently amplified to become again `N_pred`-sized.
Each parallel rank executes its own evolution process, and the results are saved in a dedicated folder, called `seedS`, where `S` is a unique index referring to each rank.

Currently, three selection criteria have been inserted: 
- selection based on the fitness of each predator (defined as the maximum consecutive overlaps between the predator and the resource), 
- selection based on the fitness, with an asymptotic fitness attributed after a threshold, 
- purely random selection, regardless of the fitness of the predator 

The code outputs a file called `histogram*.dat` for each generation `*`, which contains the list of the fitnesses of each of the survived predators, before the amplification process.
The program also outputs some files named `pred_histF_cycle_*.dat` where `F` represents the fitness and `*` the cycle. This file contains the list of the `L`-long sequences having a fitness `F` after cycle `*`. This is a way to spot the emergence of dominating sequences during an evolution process.

