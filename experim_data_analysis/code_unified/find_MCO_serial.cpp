// Written by Francesco Mambretti, 03/02/2022
// Based on find_MCO.cpp
// Finds:
// - Maximum Consecutive Overlap (MCO) and 2nd-MCO between all the strands and the given resource (target)
// - Total Mixed Overlap (TMO), computed in the same relative position of the two strands where the MCO is found
// - Largest Total Overlap (LTO), i.e. the largest set of HBs, regardless whether they are consecutive or not
//
// ssDNA target strand (COMPLEMENTED AND REVERTED!) is passed from command line as argument
// In this version, all the computations concerning the subsequences involved in overlaps and their positions
// are commented. 
// Operations have been sped up by not reverting all the input strands and checking matches; instead, only target strand is already passed as reverted and complemented! 
// This means that also the old version of the match function is not useful anymore
// Also fixed issues in TMO computation
//
// 21/03/2022 version

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

//functions
void compute_MCO_TMO_LTO(vector<char>,vector<char>,int,int,int&,int&,int&/*,int&, int&,vector<int>&,vector<int>&*/);
void compute_MCO_2nd_solo(vector<char>,vector<char>,int,int,int&,/*int, int, vector<int>,vector<int>,*/int);
//bool match(char,char);

// global variables
int tot_lines;

int main(int argc, char** argv){

	if (argc!=5)
	{
		cerr << "Please provide: a file with all the sequences to be analyzed as an argument of the executable, the number of lines in that file, the output file name where to store MCO list and the target strand (complemented & reverted)!" << endl;
		exit(EXIT_FAILURE);
	}
        
    //******************************** Settings *******************************//
    
    string filename=argv[1];
    tot_lines=atoi(argv[2]);
    string outfile_name=argv[3];
    
    //declare variables
    string curr_strand; //temporary variable
    std::vector<char> sequence;
    int MCO,MCO_2nd,TMO,LTO;
    string target=argv[4];
    std::vector<char> target_seq;
    
    // convert target strand into vector
    std::copy(target.begin(), target.end(), std::back_inserter(target_seq));

    ifstream in;  //open input file with sequences
    if (!in) {
        cout << "No input file!" << endl;
        exit(-1);
    }
    
    in.open(filename);
    ofstream out;
    out.open(outfile_name,ios::out);
    out << "MCO" << '\t' << "MCO_2nd" << '\t' << "TMO" << '\t' << "LTO" << '\n';

    int l=target_seq.size();
    int L;
 
    //******************************** Loop over input file *******************************//
    //int first_pos_def=0, second_pos_def=0;
    for (int i=0; i<tot_lines; ++i){
        in >> curr_strand;
        
        //convert string into vector
        std::copy(curr_strand.begin(), curr_strand.end(), std::back_inserter(sequence));

        L=sequence.size();
        
        compute_MCO_TMO_LTO(sequence,target_seq,L,l,MCO,TMO,LTO);//,first_pos_def,second_pos_def,first_subseq_def, second_subseq_def);
        //compute_MCO_2nd_solo must be executed after the previous function, since we need to know the MCO, to compute MCO_2nd
        compute_MCO_2nd_solo(sequence,target_seq,L,l,MCO_2nd,/*first_pos_def,second_pos_def,first_subseq_def, second_subseq_def,*/MCO);
        out << MCO << '\t' << MCO_2nd << '\t' << TMO << '\t' << LTO << endl;
        sequence.clear();
    }

    in.close(); // close input file with sequences
    out.close();

	return 0;
}

//***********************************************************//
// FUNCTIONS 


//decide pairing rules
/*bool match (char a, char b){
	if (a==b)//((a=='A' and b=='T') or (a=='T' and b=='A') or (a=='C' and b=='G')  or (a=='G' and b=='C') )
		return true;
	else
		return false;
}*/

void compute_MCO_TMO_LTO (vector<char> first_seq, vector<char> second_seq, int L, int l, int& MCO, int& TMO, int& LTO /*,int& first_pos_def, int & second_pos_def, vector<int>& first_subseq_def, vector<int>& second_subseq_def*/) { //compute the value of the max consecutive overlap between two strands

    int actual_MCO = 0, actual_TMO=0, actual_LTO=0; //it will be Max overlaps
    int overlap, tot_mix_overlap, largest_tot_overlap; //counters for overlap (consecutive) and total overlap (non consecutive)
    int p,r,k,m,j;
    bool changed=false;

    for(j=0;j<L;j++) { //loop over all the possible positions
        overlap=0; //setting to zero the overlap counter
        tot_mix_overlap=0;
        largest_tot_overlap=0;
        changed=false;
        //cycling on predator and resource nucleotides
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
            if (first_seq[p]==second_seq[r]){
                overlap++;
                tot_mix_overlap++;
                largest_tot_overlap++;
            }
            else {  //if there is a mismatch
                if (overlap > actual_MCO){
                    changed=true;
                    actual_MCO = overlap;
                }
                overlap=0;
            }
            if(changed){
                actual_TMO=tot_mix_overlap;
            }
        }//end for p,r
        
        if (overlap>actual_MCO){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            actual_MCO = overlap; //update max
            actual_TMO = tot_mix_overlap;
        }
        //check whether to update LTO
        if (largest_tot_overlap>actual_LTO)
            actual_LTO=largest_tot_overlap;
    } //end for j
    
    //now, the second sequence (length l) has the left end bases which are empty, since the first sequence left end is already at their right
    
    for(j=0;j<l;j++){
        overlap=0; //setting to zero the overlap counter
        tot_mix_overlap=0;
        largest_tot_overlap=0;
        changed=false;
        for (p=0, r=j+1; p<L && r<l; p++,r++) {
            if (first_seq[p]==second_seq[r]){
                overlap++;
                tot_mix_overlap++;
                largest_tot_overlap++;
            }
            else {  //if there is a mismatch
                if (overlap > actual_MCO){
                    changed=true;
                    actual_MCO = overlap;
                }
                overlap=0;
            }
            if(changed)
                actual_TMO=tot_mix_overlap;
        }//end for p,r
        
        if (overlap>actual_MCO){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            actual_MCO = overlap; //update max
            actual_TMO = tot_mix_overlap;
        }
        //check whether to update LTO
        if (largest_tot_overlap>actual_LTO)
            actual_LTO=largest_tot_overlap;
        
    } //end for j

    //values to be used externally
    MCO=actual_MCO;
    TMO=actual_TMO;
    LTO=actual_LTO;
    
    return ;
}

void compute_MCO_2nd_solo (vector<char> first_seq, vector<char> second_seq, int L, int l, int& MCO_2nd, /*int first_pos_def, int second_pos_def, vector<int> first_subseq_def, vector<int> second_subseq_def,*/ int MCO) { //compute the value of the max consecutive overlap between two strands

    int Max = 0; //it will be Max overlaps
    int overlap; //counter for overlap (2nd MCO)
    int j,m,p,r,k;

    for(int j=0;j<L;j++) { //loop over all the possible positions
        overlap=0; //setting to zero the overlap counter
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
            if (first_seq[p]==second_seq[r]){
                overlap++;
            }
            else {  //if there is a mismatch
                if (overlap< MCO and overlap > Max /*and (first_pos!=first_pos_def or second_pos!=second_pos_def) */){
                    Max = overlap;
                }
                overlap=0;
            }
        }//end for
        
        if (overlap< MCO and overlap>Max /*and (first_pos!=first_pos_def or second_pos!=second_pos_def)*/ ){ // necessary when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            Max = overlap; //update max
        }
    }//end for j
    
    for(j=0;j<l-1;j++){
        overlap=0; //setting to zero the overlap counter
        for (p=0, r=j+1; p<L && r<l; p++,r++) {
            if (first_seq[p]==second_seq[r]){
                overlap++;
            }
            else {  //if there is a mismatch
                if (overlap< MCO and overlap > Max){
                    Max = overlap;
                }
                overlap=0;
            }
        }
        
        if (overlap< MCO and overlap>Max){ //necessary when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            Max = overlap; //update max
        }
    }//end for j

    MCO_2nd=Max;
    
    return;
}
