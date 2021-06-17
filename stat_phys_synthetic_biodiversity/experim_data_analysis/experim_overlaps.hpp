//
//  experim_overlaps.hpp
//  
//
//  Created by Francesco Mambretti on 07/06/21.
//

#ifndef experim_overlaps_hpp
#define experim_overlaps_hpp

#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <random>
#include <fstream>
#include <map>

using namespace std;

int check_max_cons_overlap(vector<int>,vector<int>,int,int,int&);
int check_tot_overlap(vector<int>,vector<int>,int,int,int);
bool match (int,int);
vector<int> rev_and_comp(vector<int>&);

#endif /* experim_overlaps_hpp */
