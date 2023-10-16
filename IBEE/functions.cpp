#include "functions.h"

//#define PRINT

using namespace std;
using namespace arma;

// define global variables declared in global.h

//===================================================================

// random numbers
int seed_rannyu;

// -- Input parameters --//
string input_name, master_folder;
uword it;
uword cycles;
int L, l, N_pred, N_res, nbases;
bool mut;
double p_m; //mutation
uword Tot_Fit,predator, option;

//vec for simulation
umat individuals;
uvec resource;
field<string> name, initial_names;
vec fitnesses;
Mat<int> results; // matrix with 2 columns (predator, max overlap) and  N_res rows

int asymp=10; // beyond this value, saturate the cost function

bool key_inpopgen=1; // if ==1, read initial population with a w-uniform sampling from experimental data

double alpha; //fraction of random update moves


//=================================================================================================================

// FUNCTIONS

// num to char
char num_to_char(int a){
    if (a==0)
        return 'A';
    else if (a==1)
        return 'C';
    else if (a==2)
        return 'G';
    else if (a==3)
        return 'T';
};

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

//==========================
// create_folders
//==========================

void create_folders(string master_folder, int rank, int seed){
    
    string a = "mkdir "+master_folder+"/seed"+to_string(rank+1);
    cout << a << endl;
    system(a.c_str());
    string b = master_folder+"/seed"+to_string(rank+1);
    system(("cp Primes seed.in "+b).c_str()); //copy essential files for the simulation
    if (key_inpopgen==1){
        system(("cp all_experim.txt "+b).c_str());
    }
    chdir(b.c_str());
    cout << "cd "+b << endl;
    system("pwd");

    return;
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
    
    return;
}

//=============
// begin_cout
//=============
//only rank 0 should call this method
void begin_cout(int seed) {
    
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

    return;
}


//==========================
// compute_affinity
//==========================

uword compute_affinity(uword i, Random& rnd) { //compute the value of the max consecutive overlap between the target strand and the i-th predator

    uword Max = 0; //it will be Max overlap
    uword overlap = 0; //counter for overlap

    //first cycle on all possible configurations, just in one direction (not considering possible reverse option)
    //resource is long k genes and pred is l genes, so we have k(from right) +k (from left) + (n-k-1)(central possibilities)  possible overlaps with predator
    //I cycle on these starting from the left and going to right direction
    
    int p,r;
    uword j;
    
    for(j=0;j<L;j++) { //loop over all the possible positions
        overlap=0; //setting zero the overlap counter
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the predator reference frame (fixed); r: track the internal resource coordinate
            if (individuals.col(i)[p]==resource(r)){  //( match(individuals.col(i)[p],resource(r)) )
                overlap++;
            }
            else {  //if there is a mismatch
                if (overlap > Max)
                    Max = overlap;
                overlap=0;
            }
            //cout <<j << " " << p << " " << r << " " << Max << " " << num_to_char(individuals.col(i)[p]) << " " << num_to_char(resource(r)) << endl; // DEBUG
        }//end for p,r
        
        if (overlap>Max) // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            Max = overlap; //update max
        
    } //end for j
    
    for(j=0;j<l;j++){
        overlap=0; //setting to zero the overlap counter
        for (p=0, r=j+1; p<L && r<l; p++,r++) {
            if (individuals.col(i)[p]==resource(r)){//( match(individuals.col(i)[p],resource(r)) )
                overlap++;
            }
            else {  //if there is a mismatch
                if (overlap > Max)
                    Max = overlap;
                overlap=0;
            }
        }//end for p,r
        
        if (overlap>Max){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            //MCO
            Max = overlap; //update max
        }
        
    } //end for j
    
    return Max;
} //end compute_affinity

//==========================
// total_fitness
//==========================

void total_fitness(Random& rnd, bool print_0) {

    for (uword i=0;i<N_pred;i++){
        fitnesses(i) = compute_affinity(i, rnd);
        if (print_0==true)
            cout << initial_names(i) << " " << fitnesses(i) << endl;
    }
    
    #ifdef PRINT
    cout << "fitnesses: " << endl;
    fitnesses.print();
    #endif

    uvec q = find(fitnesses > 0);

    if (q.size() == 0) {
        cout << "the population has not enough members" << endl;
        exit(EXIT_FAILURE);
    }

    Tot_Fit = sum(fitnesses); //normalization
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
    double max_fit=double(fitnesses.max());
    
    if (option==0){
        for(;;) { //loop until needed
            predator = int(rnd.Rannyu(0, N_pred));
            x = rnd.Rannyu();
            if (fitnesses(predator)>3){ //discard sequences with too low-MCO
                if ( ( x <= double(fitnesses(predator))/max_fit ) && ( all(results.col(0) != predator) ) )
                    //if ( ( x <= double(fitnesses(predator))/double(Tot_Fit) ) && ( all(results.col(0) != predator) ) )
                    break;
            }
        }
    } // use maximum consecutive overlap as fitness
    
    else if (option==1){ //random uniform choice, regardless of fitness
        for(;;) {
            predator = int(rnd.Rannyu(0,N_pred));
            if (fitnesses(predator)>3){ //discard sequences with too low-MCO
                if ( ( all(results.col(0) != predator) ) )
                    break;
            }
        }
    }
    
    else if (option==2){
        for(;;) {
            predator = int(rnd.Rannyu(0,N_pred));
            if (fitnesses(predator)>3){ //discard sequences with too low-MCO
                if (fitnesses(predator)>asymp)
                    cost=asymp;
                else
                    cost=fitnesses(predator);
                x = rnd.Rannyu();
                //if ( ( x <= cost/double(Tot_Fit) ) && ( all(results.col(0) != predator) ) )
                if ( ( x <= cost/max_fit ) && ( all(results.col(0) != predator) ) )
                    break;
            }
        }
    } // use maximum consecutive overlap with saturation as fitness

    #ifdef PRINT
    cout << " the winning sequence is " << predator << endl;
    individuals.col(predator).print();
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
    for (auto el : individuals.col(pred) ){
        if (el==0)
            a += 'A';
        else if (el==1)
            a += 'C';
        else if (el==2)
            a += 'G';
        else if (el==3)
            a += 'T';
    }
    
  	name(o) = a; //save the sequence of 50 nucleotides of the predator
    #ifdef PRINT
	cout << "results" << endl;
	results.print();
    #endif
}

//==========================
// prepare_hist
//==========================

void prepare_hist(uword it) { //cycle index passed as argument

	uvec predators;
	vector<string> list_seqs;
    string folder;
    fstream outfile;
    
	for (uword i=0;i<=l;i++) { //loop over different MCO values
		predators = find( results.col(1) == i ); //select all the predators with that MCO
        
        //for each cycle, create the folder
        if (i==0/* and it==0*/){
            folder="cycle_"+to_string(it);
            system(("mkdir "+folder).c_str());
            outfile.open(folder+"/list_seqs_all.dat",ios::out); //open file for cycle it
            outfile << "sequence \t MCO \n";
        }
        
		if (predators.size() != 0) {
			for (auto& el : predators)  //loop over predators
                list_seqs.push_back(name(el));
            
			print_vector(list_seqs,folder+"/list_seqs_MCO_"+to_string(i)+".dat",0);
            
            //save on another file the full list & the omega of each sequence
            for (auto el : list_seqs)
                outfile << el << "\t" << i <<"\n";
            
            list_seqs.clear();
		}
	}
    
    outfile.close();//close stream

    return;
}

//==========================
// print_vector
//==========================

void print_vector(vector<string> v, string file,int key) {
	fstream outfile;
    if (key==0)
        outfile.open(file,ios::out);
    else
        outfile.open(file,ios::app);
	for (auto el : v) outfile << el << "\n";
	outfile.close();
}


//==========================
// initial_population
//==========================

field<string> initial_population(Random& rnd, bool key_inpopgen) {
//generate initial population of predators
    uword i;
    double x;
    field<string> name (N_pred);
    field<string> full_seq_list (1E6);
    string sequence;
    
    individuals.set_size(L,N_pred);
    
    if (key_inpopgen==0){
        if (nbases==2) {
            for (i=0;i<N_pred;i++){
                for (auto& el : individuals.col(i)) {
                    if (rnd.Rannyu() < 0.5)
                        el = 0;
                    else
                        el = 1;
                }
            }
        }//end if
        
        else if (nbases==4) {
            for (i=0;i<N_pred;i++){
                for (auto& el : individuals.col(i)) {
                    x = rnd.Rannyu();
                    if (x < 0.25){
                        el = 0;
                        name(i)+='A';
                    }
                    else if (x < 0.5 && x>=0.25){
                        el = 1;
                        name(i)+='C';
                    }
                    else if (x < 0.75 && x>=0.5){
                        el = 2;
                        name(i)+='G';
                    }
                    else{
                        el = 3;
                        name(i)+='T';
                    }
                }
            }
        }//end else if
    }
    
    else // key_inpopgen==1
    {
        uword j;
        i=0;
        uword len_exp_file=0;
        ifstream in;
        in.open("all_experim.txt"); //read sequences from external files
        if (in.good()){
            while(!in.eof()){
                in >> full_seq_list(i);
                len_exp_file++;
                i++;
            }
        }
        
        for (i=0;i<N_pred;i++){
            //x=rnd.Rannyu(0,len_exp_file-1);//randomly choose one
            //x=i;
            name(i)=full_seq_list[i];
            
            for (j=0; j<L; j++) { //initialize also the numeric sequence
                if (name[i][j]=='A')
                    individuals.col(i)(j)=0;
                else if (name[i][j]=='C')
                    individuals.col(i)(j)=1;
                else if (name[i][j]=='G')
                    individuals.col(i)(j)=2;
                else if (name[i][j]=='T')
                    individuals.col(i)(j)=3;
            }
            
        }
        in.close();
        
    }
    
    
    return name;

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
	new_population.cols(0, N_res-1) = individuals.cols(indices);
	uvec appo;
    uword num_to_create=N_pred-N_res;

	for (uword i=0;i<num_to_create;i++) {
		appo = new_population.col(int( rnd.Rannyu(0, N_res) ));
		new_population.col(N_res+i) = appo;
	}

	individuals = new_population;

    if (mut) {
	// this v is needed since individuals.col(i) does not support .begin() and .end() vector methods
        uvec v;
        for (uword i=0;i<N_pred;i++) {
            v = individuals.col(i);
            Mutation(rnd,v,L);
            individuals.col(i) = v;
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

//==========================
// rev_and_compl
//==========================

void rev_and_compl(uvec &strand){ //take care! to match with A/T, C/G binding rules --> 0 <==> 3 and 1 <==>2

    uvec strand1 = reverse(strand);
    
    for (uword i=0; i<strand1.n_elem; i++){
        if (strand1(i)==0)
            strand1(i)=3;
        else if (strand1(i)==3)
            strand1(i)=0;
        else if (strand1(i)==2)
            strand1(i)=1;
        else if (strand1(i)==1)
            strand1(i)=2;
    }
    
    strand=strand1;

    return;
}
