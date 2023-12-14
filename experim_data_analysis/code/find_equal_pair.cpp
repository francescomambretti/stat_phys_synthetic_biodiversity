// Written by Francesco Mambretti, 21/02/2022 - simplified version of the find_MCO.cpp code within my oxDNA utilities
// Finds the number of identical bases between two DNA sequences passed by command line
// Note: if there are two overlapping regions with the same overlap, the code only finds one of them
// It corresponds to the MCO between a sequence and the complemented&reverted of a second sequence
//
// 21/03/2022 version

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

//functions
int compute_equal_bases(vector<char>,vector<char>);

int main(int argc, char** argv){

	if (argc!=3)
	{
		cerr << "Please execute like this: " << argv[0] << " sequence_1 sequence_2" << endl;
		exit(EXIT_FAILURE);
	}	

	string first=argv[1];
    string second=argv[2];
	
	//convert strings into vectors
	std::vector<char> first_seq;
    std::copy(first.begin(), first.end(), std::back_inserter(first_seq));

	std::vector<char> second_seq;
    std::copy(second.begin(), second.end(), std::back_inserter(second_seq));

    int equal_bases;
    equal_bases = compute_equal_bases(first_seq,second_seq);
    cout<<equal_bases;
    
	return 0;
}

//***********************************************************//
// FUNCTIONS 

int compute_equal_bases (vector<char> first_seq, vector<char> second_seq) { //compute the value of the max consecutive overlap between two strands

    int Max = 0; //save here current maximum
    int overlap; //counter for overlap
	int L=first_seq.size();
	int l=second_seq.size();
    int p,r,j;
    
    for(j=0;j<L;j++) { //loop over all the possible positions
        overlap=0; //setting to zero the overlap counter
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: tracks the left end of the first sequence in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
            if ( first_seq[p]==second_seq[r] ){
                overlap++;
            }
            else {
                if (overlap > Max){
                    Max = overlap;
                }
                overlap=0;
            }
		}//end for p,r
        
        if (overlap>Max){ //necessary when there is a group of matching bases at the last loop, which is not ended by a non-matching pair
            Max = overlap; //update max
		}
	} //end for j
    
    //now, the second sequence (length l) has the left end bases which are empty, since the first sequence left end is already at their right
    for(j=0;j<l-1;j++){
        overlap=0; //setting to zero the overlap counter
        for (p=0, r=j+1; p<L && r<l; p++,r++) {
            if (first_seq[p]==second_seq[r]){
                overlap++;
            }
            else {  //if there is a mismatch
                if (overlap > Max){
                    Max = overlap;
                }
                overlap=0;
            }
        }
        
        if (overlap>Max){ //necessary when there is a group of matching bases at the last loop, which is not ended by a non-matching pair
            Max = overlap; //update max
        }
    }//end for j
    
	return Max;
}


                            
