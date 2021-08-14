//global variables - declared here, defined elsewhere

#ifndef __global_
#define __global_

#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

//random numbers
extern int seed_rannyu;

//-- Input parameters --//
extern string input_name;
extern uword it; //index
extern uword cycles; //number of generations
extern int L, l, N_pred, N_res, nbases; //length of predator, target; number of predators and resources; number of possible nucleotides
extern bool mut; //mutations
extern double p_m;
// useful for simulation
extern uword Tot_Fit, predator, option;

//vec for simulation
extern umat filaments;
extern uvec resource;
extern field<string> name;
extern uvec fitnesses;
extern Mat<int> results; // matrix with 2 columns (predator, max overlap) and N_res rows


#endif
