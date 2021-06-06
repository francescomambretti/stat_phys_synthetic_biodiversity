//
//  max_and_tot_overlap.cpp
//  
//
//  Created by Francesco Mambretti on 04/06/21.
//  Generate a population of L-long DNA strands and subsequently check their (largest) maximum consecutive overlap with target strand. Compute also the largest total overlap with the target strand, aiming to plot them together. "Largest" means that all possible relative positions between the target and the predator are checked.
//  Built on the basis of null_model_constraint
//  Matching bases: 0-1 and 2-3
//  Accepts as input the number of strands to be generated and a random seed
//

#include "max_and_tot_overlap.hpp"

using namespace std;

int main(int argc, char** argv){
    
    if (argc!=3){
        cerr << "Error! Correct usage: " << argv[0] << " number of strands to generate      myseed" << endl;
        exit(EXIT_FAILURE);
    }
    
    long int pop_size=atoi(argv[1]);
    int myseed = atoi(argv[2]);
    
    int l=20; //target length
    int L=50; //predator length
    int a,b; //for loops indexes
    
    std::default_random_engine generator; // start random number generator
    generator.seed(myseed);
    std::uniform_int_distribution<> distribution(0,3);
    
    vector<int> target (l);  //allocate the target
    for (a=0; a < l; ++a){
        target[a]=distribution(generator);
        //cout << a << '\t' << target[a] << endl;
    }
        
    vector<int> predator (L);
    
    //generate population of random strands and check their max consecutive overlap with the target
    
    int appo, appo2;
    ofstream outFile("all_max_cons_overlaps.dat");
    ofstream outFile2("all_tot_overlaps.dat");
    
    for (a=0; a<pop_size; ++a){
        for (b=0; b<L; ++b){ //allocate a random predator
            predator[b]=distribution(generator);
            //cout << b << '\t' << predator[b] << endl;
        }
        
        appo=check_max_cons_overlap(target,predator,l,L);
        appo2=check_tot_overlap(target,predator,l,L);
        //cout << appo << endl;
        outFile << appo << endl;
        outFile2 << appo2 << endl;
    }
    
    outFile.close();
    outFile2.close();
    
    system("python compare.py");
    
    return 0;
}


// FUNCTIONS ==================================================================================

//======================================
// check_max_cons_overlap
//======================================

int check_max_cons_overlap(vector<int> target, vector<int> predator, int l, int L){
    
    int my_max = 0, temp_max;
    int p,q; //counter
    vector<int> overlaps (l,0);
    vector<int> overlaps_q (l,0);
    int i,j;
    int limit_loop = L-l+1;
    
    // overlap with resource completely over predator
    for (i=0;i< limit_loop;++i) {
        fill(overlaps.begin(), overlaps.end(), 0);
        p=0; //it counts consecutive overlaps
        for (j=i;j<l+i;++j) { // internal target index, always take l values
            if (match (target[j-i], predator[j])) //they match
                overlaps[p]++;
            else p++;
        }
        temp_max=*max_element(overlaps.begin(), overlaps.end());
        if (temp_max > my_max)
            my_max = temp_max;
    }
    
    // overlap with resource partly beyond predator bounds
    // left and right "excess parts" are treated exploiting symmetry
    for (i=0;i<l-1;++i) {
        fill(overlaps.begin(), overlaps.end(), 0);
        fill(overlaps_q.begin(), overlaps_q.end(), 0);
        p=0;q=0;
        for (j=0;j<=i;++j) {
            if ( match(target[l-1-j], predator[i-j]))
                overlaps[p]++;
            else p++;
            
            if ( match(target[j], predator[L-1-i+j]))
                overlaps_q[q]++;
            else q++;
        }
        temp_max=*max_element(overlaps.begin(), overlaps.end());
        if (temp_max > my_max)
            my_max = temp_max;
        
        temp_max=*max_element(overlaps_q.begin(), overlaps_q.end());
        if (temp_max > my_max)
            my_max = temp_max;
     }
    
    return my_max;
}

//======================================
// check_tot_overlap
//======================================

int check_tot_overlap(vector<int> target, vector<int> predator, int l, int L){
    
    int my_max = 0;
    int tot=0, tot_2=0;
    int i,j;
    int limit_loop = L-l+1;
    
    // overlap with resource completely over predator
    for (i=0;i< limit_loop;++i) {
        tot=0;
        for (j=i;j<l+i;++j) { // internal target index, always take l values
            if (match (target[j-i], predator[j])) //they match
                tot++;
        }
        
        if (tot > my_max)
            my_max = tot;
    }
    
    // overlap with resource partly beyond predator bounds
    // left and right "excess parts" are treated exploiting symmetry
    for (i=0;i<l-1;++i) {
        tot=0;
        tot_2=0;
        for (j=0;j<=i;++j) {
            if ( match(target[l-1-j], predator[i-j]))
                tot++;
            
            if ( match(target[j], predator[L-1-i+j]))
                tot_2++;
        }
        
        if (tot > my_max)
            my_max = tot;
        
        if (tot_2 > my_max)
            my_max = tot_2;
     }
    
    return my_max;
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
