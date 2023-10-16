// Original code by Nicolò Pedrani
// Modified by Francesco Mambretti - 06/07/2021

// Genetic/Evolutive algorithm. Simulates the evolution of a population of sequences of nucleotides, which compete for a limited set of resources (fixed target sequence, unique). The evolution can be a mix of fitness-driven and random contributions, with the alpha parameter, 0<= alpha <= 1, regulating the amount of random selection process, being 1-alpha the amount of driven selection (i.e. option=0 or option=2)
// Prints also the identities and the fitnesses of the individuals in the initial population.
// The initial population built by any seed is the same, unless key_inpopgen is set, to build the initial population with a uniform sampling in MCO, reading external (experimental) files.
// Sequences with MCO<cut are automatically discarded.

// To run: // mpiexec -np $PROCS ./evolution    input.dat   $SEED   0   results_folder
// Currently, mutation is off by default

// 30/11/2022 version

#include "functions.h"
#include "global.h"
#include "mpi.h"
#include <unistd.h>
#include <stdio.h>

//#define PRINT  //for some output

using namespace std;
using namespace arma;

int main (int argc, char *argv[]) {

    // initialization
    
    //MPI setup
    int procs,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status stat;

    if (argc!=5) {
        cerr << "insert input file name, random generator seed, if you want (1) or not (0) mutations during population update and the master folder for results" << endl;
        return -1;
    }

    input_name = argv[1];
    seed_rannyu = atoi(argv[2]);
    mut = atoi(argv[3]);
    master_folder=argv[4];

    int seed = 1; //for initial population generation - give the same to all ranks
    Read_input(); //each rank reads its own input file
    cout << "my name is rank " << rank << " of " << procs << " and my seed is " << seed << endl;
    
    //variables declaration & sizing
    results.set_size(N_res,2); //save the outcomes of a selection inside a generation: predator id & its fitness. N_res predators survive to each generation.
    name.set_size(N_res); // save the sequences of the winning 50mers
    fitnesses.set_size(N_pred); // save the fitness of each individual at each generation (overwritten)

    //each rank creates its own folder, copies there the important files and works inside it
    create_folders(master_folder,rank,seed);
    
    Random rnd = random_initialization(seed);
    arma_rng::set_seed(123); // fix the seed for this, in order to give the same target resource to all the parallel ranks

    if (rank == 0)
        begin_cout(seed);
    
    uword iprint =  N_res/10; //print every * steps
 
    //generate target strand
    resource = {1,2,2,3,0,3,3,2,2,0,1,1,1,3,1,2,1,0,3,2}; //"CGGTATTGGACCCTCGCATG"; //randi<uvec>(l,distr_param(0,nbases-1)); // randi: generate object with random integer values in specified interval. uvec is a vector of unsigned integers. generate l elements, whose value is chosen between 0 and nbases-1. distr_param lives inside armadillo
    rev_and_compl(resource); //reverse-complement the resource
    
    cout << "resource sequence: " << endl;
    cout << resource << endl;
    
    //randomly generate initial population and print it to file
    initial_names=initial_population(rnd,key_inpopgen);
    
    #ifdef PRINT
    cout << "Population" << endl;
    individuals.print();
    #endif
    
    if (rank==0)
        cout << "number of cycles: " << cycles << endl;
 
    wall_clock timer;
    timer.tic();
    ofstream outfile;

    // added: print histogram before selection cycles
    system("pwd");
    system("mkdir cycle_0");
    total_fitness(rnd,false); //compute the fitness of each individual
    
    //print initial pop. and their fitnesses
    outfile.open("./initial_pop.dat", ios::out);
    
    //sort by fitness
    uvec sorting_indices = sort_index(fitnesses);
    fitnesses=sort(fitnesses);
    
    for (int i=0; i < N_pred; i++)
        outfile << initial_names[sorting_indices[i]] << " " << fitnesses[i] << endl;
    
    outfile.close();
    
    outfile.open("./cycle_0/histogram.dat",ios::out);
    uvec h1 = hist(fitnesses, linspace<vec>(0,l+1,l+1));
    vec h2 = arma::conv_to<vec>::from(h1);
    h2 = h2/sum(h2); //normalize between 0 and 1
    
    if (outfile.is_open()){
        h2.print(outfile);
    }
    else{
        cerr << "Error in opening file \n";
        exit(-1);
    }
    outfile.close();
    
    uword i = 0;
    double r; // temporary random variable
    vec v;
    string command;
    
    //now, reset a different seed for any trajectory
    seed=seed_rannyu*(rank+1); // run several parallel simulations with different random seeds
    rnd = random_initialization(seed);
    
    for (it = 0; it<cycles; it++) {
        results.fill(-1);
        command = "mkdir cycle_"+to_string(it+1);
        system(command.c_str());
        
        outfile.open("./cycle_"+to_string(it+1)+"/histogram.dat",ios::out);
        total_fitness(rnd,false); //compute the fitness of each individual
        for(i=0;i< N_res;i++) {
         #ifdef PRINT
            cout << "selection and lunch" << endl;
         #endif
            r=rnd.Rannyu();
            if (r>alpha) // apply driven selection criterion, usually option!=1
                Selection(rnd,option); //decide which predator to be selected at the i-th algorithm step
            else // perform purely random selection
                Selection(rnd,1); //decide which predator to be selected at the i-th algorithm step
            
			Lunch_time(rnd,i,predator); //save data of the winners inside results matrix
			if (i%iprint == 0 and rank==0)
                cout << "genetic algorithm step: " << i << endl;
        }
        v=conv_to<vec>::from(results.col(1)) ;
        h1 = hist(v, linspace<vec>(0,l+1,l+1));
        h2 = arma::conv_to<vec>::from(h1);
        h2 = h2/sum(h2); //normalize between 0 and 1
        h2.print(outfile);
        prepare_hist(it);
        
        #ifdef PRINT
        if (rank==0)
            cout << "Update" << endl;
        #endif
        
        Update(rnd); //amplify winning population for a new algorithm iteration
        
        if (rank==0)
            cout << "---------# END CYCLE # " << it << "  ----------- " << endl;
        outfile.close();
    }

    chdir(".."); //exit from directory
    MPI_Finalize();

    double time = timer.toc();
    if (rank==0)
        cout << "time spent: " << time <<" seconds" <<endl;
    rnd.SaveSeed();
    
    return 0;
}


