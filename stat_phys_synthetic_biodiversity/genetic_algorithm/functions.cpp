#include "functions.h"

//#define PRINT

using namespace std;
using namespace arma;

// define global variables declared in global.h

//===================================================================

// random numbers
int seed_rannyu;

// -- Input parameters --//
string input_name;
uword it;
uword cycles;
int L, l, N_pred, N_res, nbases;
bool mut;
double p_m;
// useful for simulation
uword Tot_Fit,predator, option;

//vec for simulation
umat filaments;
uvec resource;
field<string> name;
uvec fitnesses;
Mat<int> results; // matrix with 2 columns (predator, max overlap) and  N_res rows

int asymp=12; // beyond this value, saturate the cost function

double alpha;

//=================================================================================================================

// FUNCTIONS

//==========================
// random_initialization
//==========================

Random random_initialization(int rank_index) {

    Random rnd; //this object will be returned at the end of the function
    
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for (int k=0; k < rank_index; k++) //assign different p1 and p2 values to each process
            Primes >> p1 >> p2 ;
        }
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                //cout << seed[0] << '\t' << seed[1] << '\t' << seed[2] << '\t' << seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    }
    else
        cerr << "PROBLEM: Unable to open seed.in" << endl;

   return rnd;
}

//=============
// Read_input
//=============

void Read_input() {

    fstream ReadInput(input_name,ios::in); //open file with input parameters

    ReadInput >> cycles; // generations
    ReadInput >> L; // number of bases for predators
    ReadInput >> l; // number of bases for resource
    ReadInput >> N_pred; // predators population size
    ReadInput >> N_res; // resources population size
	ReadInput >> nbases; // use 0,1 or 0,1,2,3
    ReadInput >> p_m; // mutation probability (ignored if mutations are turned off)
    ReadInput >> option; // which type of cost function to use
    ReadInput >> alpha; // fraction of purely random selection

    ReadInput.close();
}

//==========================
// compute_affinity
//==========================

uword compute_affinity(uword i, Random& rnd) { //compute the value of the max consecutive overlap between the target strand and the i-th predator

    uword Max = 0; //it will be Max overlaps
    uword overlap; //counter for overlap
    
    //first cycle on all possible configurations, just in one direction (not considering possible reverse option)
    //resource is long k genese and pred is l genese, so we have k(from right) +k (from left) + (n-k-1)(central possibilities)  possible overlaps with predator
    //I cycle on these starting from the left and going to right direction
    
    uword limit = L+l-1; //rightmost possible coordinate of a nucleotide
    int p,r;
    
    for(uword j=0;j<limit;j++) { //loop over all the possible positions
        overlap=0; //setting zero the overlap counter
        //cycling on predator and resource genes that are overlapping
        //with two variables.
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the predator reference frame (fixed); r: track the internal resource coordinate
            if (p<0) continue; //safe check
            else {
                if ( match(filaments.col(i)[p],resource(r)) )
                    overlap++;
                else {
                    if (overlap > Max)
                        Max = overlap;
                    overlap=0;
                }
            }
        }//end for p,r
        
        if (overlap>Max)
            Max = overlap; //update max
        
    } //end for j
    
    return Max;
} //end compute_affinity

//==========================
// total_fitness
//==========================

void total_fitness(Random& rnd) {

    for (uword i=0;i<N_pred;i++)
        fitnesses(i) = compute_affinity(i, rnd);

    #ifdef PRINT
    cout << "fitnesses: " << endl;
    fitnesses.print();
    #endif

    uvec q = find(fitnesses > 0);

    if (q.size() == 0) {
        cout << "the population has not enough members" << endl;
        exit(EXIT_FAILURE);
    }

    Tot_Fit = sum(fitnesses); //probably useless
    #ifdef PRINT
    cout << "Total fitness " << Tot_Fit << endl;
    #endif
    
    return;
}


//==========================
// Selection
//==========================

void Selection(Random& rnd, uword option) {

    double x, cost;
    
    if (option==0){
        for(;;) {
            predator = int(rnd.Rannyu(0, N_pred));
            x = rnd.Rannyu();
            if ( ( x <= double(fitnesses(predator))/double(l) ) && ( all(results.col(0) != predator) ) )
                break;
        }
    } // use maximum consecutive overlap as fitness
    
    else if (option==1){ //random uniform choice, regardless of fitness
        for(;;) {
            predator = int(rnd.Rannyu(0,N_pred));
            if ( ( all(results.col(0) != predator) ) )
                break;
        }
    }
    
    else if (option==2){
        for(;;) {
            predator = int(rnd.Rannyu(0,N_pred));
            if (fitnesses(predator)>asymp)
                cost=asymp;
            else
                cost=fitnesses(predator);
            
            x = rnd.Rannyu();
            if ( ( x <= cost/double(l) ) && ( all(results.col(0) != predator) ) )
                break;
        }
    } // use maximum consecutive overlap with saturation as fitness

    #ifdef PRINT
    cout << " the winning sequence is " << predator << endl;
    filaments.col(predator).print();
    #endif

    return;
}

//==========================
// Lunch_time
//==========================


void Lunch_time(Random& rnd, uword o, uword pred) {

// Now we have to save predator and the relative value of the cost function
    results(o,0) = pred;  //o indexes the iteration index between 0 and N_res (i.e. the number of individuals which survive to each selection cycle, before amplification
	results(o,1) = fitnesses(pred);
	//cout << "max overlap: " << fitnesses(pred) << endl;
// and the species name	
	string a="";
  	for (auto el : filaments.col(pred) )
        a += to_string(el);
  	name(o) = a; //save the sequence of 50 nucleotides of the predator
    #ifdef PRINT
	cout << "results" << endl;
	results.print();
    #endif
}

//==========================
// prepare_hist
//==========================

void prepare_hist(uword it) {

	//cout << "------------ HISTOGRAMS --------------" << endl;
	uvec predators;
	vector<string> hist;
	//results.print();
	for (uword i=0;i<=l;i++) {
		predators = find( results.col(1) == i );
		if (predators.size() != 0) {
			//cout << "predators for maximum overlap of: " << i << endl;
			//predators.print();
			for (auto& el : predators) 
				hist.push_back(name(el));
			print_vector(hist,"pred_hist"+to_string(i)+"_cycle_"+to_string(it)+".dat");
			hist.clear();
		}
	}

}

//==========================
// print_vector
//==========================

void print_vector(vector<string> v, string file) {
	fstream fd;
	fd.open(file,ios::out);
	for (auto el : v) fd << el << "\n";
	fd.close();
}


//==========================
// initial_population
//==========================

void initial_population(Random& rnd) {
//generate initial population of predators
    uword i;
    double x;
    
    filaments.set_size(L,N_pred);
    if (nbases==2) {
        for (i=0;i<N_pred;i++){
            for (auto& el : filaments.col(i)) {
                if (rnd.Rannyu() < 0.5)
                    el = 0;
                else
                    el = 1;
            }
        }
    }//end if
    
    else if (nbases==4) {
        for (i=0;i<N_pred;i++){
            for (auto& el : filaments.col(i)) {
                x = rnd.Rannyu();
                if (x < 0.25)
                    el = 0;
                else if (x < 0.5 && x>=0.25)
                    el = 1;
                else if (x < 0.75 && x>=0.5)
                    el = 2;
                else
                    el = 3;
            }
        }
    }//end else if
    
}//end initial_population

//==========================
// match
//==========================

//in DNA object A-T or T-A, and G-C or C-G
//in this case 1-0 or 0-1, and 2-3 or 3-2
bool match(uword x,uword y) {
	if ( (x==0 and y==1) or (x==1 and y==0) or (x==2 and y==3) or (x==3 and y==2))
		return true;
	else
		return false;
}


//==========================
// Update
//==========================

void Update(Random& rnd) {

    //take the winning set of a given algorithm cycle and amplify them for the new generation
    umat new_population(L,N_pred);
	uvec indices = conv_to<uvec>::from(results.col(0));
	new_population.cols(0, N_res-1) = filaments.cols(indices);
	uvec appo;
    uword num_to_create=N_pred-N_res;

	for (uword i=0;i<num_to_create;i++) {
		appo = new_population.col(int( rnd.Rannyu(0, N_res) ));
		new_population.col(N_res+i) = appo;
	}

	filaments = new_population;

    if (mut) {
	// this v is needed since filaments.col(i) does not support .begin() and .end() vector methods
        uvec v;
        for (uword i=0;i<N_pred;i++) {
            v = filaments.col(i);
            Mutation(rnd,v,L);
            filaments.col(i) = v;
        }
    }
	
} //end Update


//==========================
// Mutation
//==========================

void Mutation(Random& rnd,uvec& v, int genes) {
//diversi tipi di mutazione
/*
      // new predator	
      double r = rnd.Rannyu();
      if (r<p_m && nbases==2){
	      for (auto& el : v) {
		      if (rnd.Rannyu() < 0.5)
			      el = 0;
		      else el=1;
	      }
      }
      if (r<p_m && nbases==4){
	      for (auto& el : v) {
		      double x = rnd.Rannyu();
	      	         if (x < 0.25)
				 el = 0;
		 	 else if (x < 0.5 && x>=0.25)
				 el = 1;
			 else if (x < 0.75 && x>=0.5)
				 el = 2;
			 else 
				 el = 3;
	      }
      }
      //random
      r = rnd.Rannyu();
      if (r < p_m)
		v = shuffle(v); //questo non lo modifica
      //shift di l posizioni per m città contigue
      r = rnd.Rannyu();
      if (r < p_m ) {
	        int m = int ( rnd.Rannyu(1,genes-1) );//quanti shiftare, parto sempre dal secondo, non l'ultimo
		int ll = int ( rnd.Rannyu(1, genes - m ) );//lo shift, con un certo limite
		shift_vector(v,m,L);
      }
      r = rnd.Rannyu();
      //singolo scambio di 2 città
      if (r < p_m) {
	        int m = int (rnd.Rannyu(1,genes-1) );
		int ll = int (rnd.Rannyu(1,genes-1) );
		while (L == m) {
            ll = int (rnd.Rannyu(1,genes-1) );
		}	
	        swap(* (v.begin()+m) ,*( v.begin()+L) );
      }
      r = rnd.Rannyu();
      //inversione dell'ordine delle città in un certo range
      if (r < p_m) {
	        int m = int (rnd.Rannyu(1,genes-1));//inizio, non può essere il primo o l'ultimo
	        int ll = int (rnd.Rannyu(m+1,genes));//fine, può essere anche l'ultimo
	        reverse(v.begin()+m,v.begin()+L);
      }*/
}

//==========================
// shift_vector
//==========================

void shift_vector(uvec &v, int m, int n) {
    /*int i,j;
    
	for (i = m; i > 0 ; i-- ) {
		for (j = 0; j < n; j++ ) {
			swap (*(v.begin()+i+j) , *(v.begin()+i+j+1) );
		}
	}*/
}

