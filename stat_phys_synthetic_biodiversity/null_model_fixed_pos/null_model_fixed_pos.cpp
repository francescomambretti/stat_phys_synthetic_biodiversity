//
//  null_model_fixed_pos.cpp
//  
//
//  Created by Francesco Mambretti on 27/05/21.
//  Generate a population of L-long DNA strands and subsequently check their maximum consecutive overlap with target strand at fixed relative position (do not slide one on the top of the other)
//  Matching bases: 0-1 and 2-3
//

#include "null_model_fixed_pos.hpp"

using namespace std;

int main(int argc, char** argv){
    
    if (argc!=2){
        cerr << "Error! Correct usage: " << argv[0] << " number of strands to generate" << endl;
        exit(EXIT_FAILURE);
    }
    
    long int pop_size=atoi(argv[1]);
    
    int l=20; //target length
    int L=50; //predator length
    int a,b; //for loops indexes
    
    std::default_random_engine generator; // start random number generator
    generator.seed(9);
    std::uniform_int_distribution<> distribution(0,3);
    std::uniform_int_distribution<> distribution_pos(-l+1,L-1); //possible relative positions
    
    vector<int> target (l);  //allocate the target
    for (a=0; a < l; ++a){
        target[a]=distribution(generator);
        //cout << a << '\t' << target[a] << endl;
    }
        
    vector<int> predator (L);
    
    //generate population of random strands and check their max consecutive overlap with the target
    
    int appo, pred_orig;
    ofstream outFile("all_max_cons_overlaps_fixed_pos.dat");
    
    for (a=0; a<pop_size; ++a){
        for (b=0; b<L; ++b){ //allocate a random predator
            predator[b]=distribution(generator);
            //cout << b << '\t' << predator[b] << endl;
        }
        //uniformly choose the relative position of the target and the predator
        pred_orig=distribution_pos(generator);
        appo=check_max_cons_overlap_fixed_pos(target,predator,l,L,pred_orig);
        //cout << appo << endl;
        outFile << appo << endl;
    }
    
    //outFile.close();
    
    return 0;
}


// FUNCTIONS ==================================================================================

//======================================
// check_max_cons_overlap_fixed_pos
//======================================

int check_max_cons_overlap_fixed_pos(vector<int> target, vector<int> predator, int l, int L, int pred_orig){
    
    int p; //counter
    int i,j, start;
    int limit_loop;
    vector<int> overlaps (0,0);
    
    if (pred_orig < 0){
        start=0;
        limit_loop = pred_orig+l; //the target ends l-1 positions after its origin, which is in pred_orig
        overlaps.resize(limit_loop,0); // number of bases which may effectively overlap
    }
    
    else if (pred_orig >=0 and pred_orig <= L-l) { //fixed predator, target completely within predator boundaries
        start = pred_orig;
        limit_loop = pred_orig+l;
        overlaps.resize(limit_loop,0); // number of bases which may effectively overlap
    }
    
    else{ //fixed predator, target exceeds at the right
        start = pred_orig;
        limit_loop=L-1;
        overlaps.resize(L-pred_orig,0); // number of bases which may effectively overlap
    }
    
    fill(overlaps.begin(), overlaps.end(), 0);
    p=0; //it counts consecutive overlaps
    for (j=start; j < limit_loop; j++){
        if (match (target[j-i], predator[j])) //they match
            overlaps[p]++;
        else p++;
    }
    
    return *max_element(overlaps.begin(), overlaps.end());

}

//=============
// match
//=============

bool match (int x, int y){
    
    if (x<2){ // x=0 or 1
        if (y<2){ // y=0 or 1
            if (x!=y)
                return true;
            else
                return false;
        }
        else
            return false;
    }
    
    else { //x = 2 or 3
        if (y>1){
            if (x!=y)
                return true;
            else
                return false;
        }
        else
            return false;
    }
}
