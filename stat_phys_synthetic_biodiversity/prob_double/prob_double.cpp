//
//  prob_double.cpp
//  
//
//  Created by Francesco Mambretti on 27/02/21.
//
//  Computes the probability of picking up a repeated individual in a pool of N strings, divided in n species with r=N/n members each, if one samples k elements from this pool.
//  The process is repeated ntrials times

#include "prob_double.hpp"

using namespace std;

int main(){
    long long int N=200E9;//200E9
    long long int n=25E9;//25E9
    long long int r=8;
    int k=1E6;//1E6
    int ntrials=20;
    int s,t, index;
    long long int number;
    long long int label;
    
    //random number
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long long int> distrib(0, N-1);
    
    vector<long long int> sample(k);
    
    for (t=0; t < ntrials; ++t){
        for (s=0; s < k; ++s){
            number=distrib(gen); //random selection of a number between 0 and N
            label=(number/r);
           // cout << number << '\t' << label << '\n';
            sample[s]=distrib(gen)/r;
        }
        cout  <<"trial " << t<< " "<< count_rep(sample,k) <<endl;
    }
    
    return 0;
}

//==================================
//functions

//==================================
// fill_population
//==================================

void fill_population(vector<long int>& population, int N, int r, int n){ //generate initial population
    int i,j;
    for (i=0; i < n; ++i){
        for (j=0; j < r; ++j){
            population[r*i+j]=i; //set r members for each of the n species
        }
    }
    
    return;
}

//==================================
// count_rep
//==================================

int count_rep(vector<long long int> sample, int k){  //find the number of repeated individuals
    int i,j;
    int found=0;
    
    //sort vector
    sort(sample.begin(),sample.end());
    
    for (i=0; i < k-1; ++i){
        for (j=i+1; j < k; ++j){
            if (sample[j]==sample[i]) //found
            {
                //cout << "found! " << sample[i] << "   " << sample[j] << endl;
                found++;
            }
        }
    }
    return found;
}
