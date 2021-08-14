// Original code by Nicolò Pedrani
// Modified by Francesco Mambretti - 06/07/2021

// Genetic algorithm. Simulates the evolution of a population of ssDNAs, which have to compete for a limited set of resources.

// To run: // mpiexec -np $PROCS ./evolution input.dat $SEED 0

// Currently, mutation is off by default

// 14/08/2021 version

#include "functions.h"
#include "global.h"
#include "mpi.h"
#include <unistd.h>
#include <stdio.h>

//#define PRINT  //for output

using namespace std;
using namespace arma;

int main (int argc, char *argv[]) {

    //MPI setup
    int procs,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status stat;

    if (argc!=4) {
        cerr << "insert input file name, random generator seed and if you want (1) or not (0) mutations during population update" << endl;
        return -1;
    }

    input_name = argv[1];
    seed_rannyu = atoi(argv[2]);
    mut = atoi(argv[3]);

    int seed = seed_rannyu*(rank+1); // run several parallel simulations with different random seeds
    Read_input(); //each rank reads its own input file
    cout << "my name is rank " << rank << " of " << procs << " and my seed is " << seed << endl;

    //each rank creates its own folder, copies there the important files and works inside it
    string a = "mkdir seed"+to_string(seed);
    cout << a << endl;
    system(a.c_str());
    string b = "seed"+to_string(seed);
    chdir(b.c_str());
    cout << "cd "+b << endl;
    //system("cp ../analysis.ipynb ../Primes ../seed.in ."); //copy essential files for the simulation
    system("cp ../Primes ../seed.in ."); //copy essential files for the simulation
    system("pwd");

    //###########################################################
    Random rnd = random_initialization(seed);
    arma_rng::set_seed(2); // fix the seed for this, in order to give the same target to all the parallel ranks
    //###########################################################

    if (rank == 0) {
        cout << "------------------------" << endl;
        cout << "number of different bases: " << nbases <<endl;
        cout << "genes per predator: " << L  << endl;
        cout << "genes per resource: " << l << endl;
        cout << "population size: " <<  N_pred<< endl;
        cout << "resource population size: " <<  N_res << endl;
        cout << "probability mutation: "<< p_m << endl;
        cout << "random seed generator: " << seed << endl;
        cout << "------------------------" << endl;
        cout << endl;

        if (mut)
            cout << "you have chosen to do mutations during update" << endl;
    }
    
    uword iprint =  N_res/10; //print every 10 steps
 
    resource = randi<uvec>(l,distr_param(0,nbases-1)); // randi: generate object with random integer values in specified interval. uvec is a vector of unsigned integers. generate l elements, whose value is chosen between 0 and nbases-1. distr_param lives inside armadillo
    initial_population(rnd);  //randomly generate initial population

    #ifdef PRINT
    cout << "Population" << endl;
    strands.print();
    cout << "resource genes: " << endl;
    resource.print();
    #endif
    results.set_size(N_res,2); //save the outcomes of a selection inside a generation: predator id & its fitness. N_res predators survive to each generation.
    name.set_size(N_res); // save the sequences of the winning 50mers
    fitnesses.set_size(N_pred); // save the fitness of each individual at each generation (overwritten)

    if (rank==0)
        cout << "number of cycles: " << cycles << endl;
 
    wall_clock timer;
    timer.tic();
    fstream out;
    cout << endl;

    uword i = 0;
    for (it = 0; it<cycles; it++) {
        results.fill(-1);
        out.open("histogram"+to_string(it)+".dat",ios::out);
        total_fitness(rnd); //compute the fitness of each individual
        for(i=0;i< N_res;i++) {
         #ifdef PRINT
            cout << "selection and lunch" << endl;
         #endif
			Selection(rnd,option); //decide which predator to be selected at the i-th algorithm step
			Lunch_time(rnd,i,predator); //save data of the winners inside results matrix
			if (i%iprint == 0 and rank==0)
                cout << "genetic algoritm step: " << i << endl;
        }
        results.col(1).print(out);
        prepare_hist(it);
        #ifdef PRINT
        cout << "Update" << endl;
        #endif
        Update(rnd); //amplify winning population for a new algorithm iteration
        if (rank==0)
            cout << "---------# END CYCLE # " << it << "  ----------- " << endl;
        out.close();
    }

    chdir(".."); //exit from directory
    MPI_Finalize();

    double time = timer.toc();
    if (rank==0)
        cout << "time spent: " << time <<" seconds" <<endl;
    rnd.SaveSeed();
    
    return 0;
}
