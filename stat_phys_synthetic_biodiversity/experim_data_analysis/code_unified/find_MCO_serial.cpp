// Written by Francesco Mambretti, 03/02/2022
// Based on find_MCO.cpp
// Finds the MCO and the 2nd-MCO between all the strands and the given resource (target)
// Performs also the calculation of Total Mixed Overlap (TMO), computed in the same relative position of the two strands where the MCO is found
// Added the calculation of Largest Total Overlap (LTO), i.e. the maximum set of total HBs
// ssDNA target strand is passed from command line as argument
// 07/02/2022 version

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

//functions
//void compute_MCO_TMO(vector<char>,vector<char>,int&,int&,int&,int&);
void compute_MCO_TMO_LTO(vector<char>,vector<char>,int&,int&,int&,int&, int&,vector<int>&,vector<int>&);
void compute_MCO_2nd_solo(vector<char>,vector<char>,int&,int, int, vector<int>,vector<int>,int);
bool match(char,char);

// global variables
int tot_lines;

int main(int argc, char** argv){

	if (argc!=5)
	{
		cerr << "Please provide: a file with all the sequences to be analyzed as an argument of the executable, the number of lines in that file, the output file name where to store MCO list and the target strand!" << endl;
		exit(EXIT_FAILURE);
	}
    
    
    //******************************** Settings *******************************//
    
    string filename=argv[1];
    tot_lines=atoi(argv[2]);
    string outfile_name=argv[3];
    
    //declare variables
    string curr_strand; //temporary variable
    std::vector<char> sequence;
    int MCO, MCO_2nd,TMO,LTO;
    string target=argv[4];
    std::vector<char> target_seq;
    
    // convert target strand into vector
    std::copy(target.begin(), target.end(), std::back_inserter(target_seq));
    
    //the aim is to compute the max cons overlap between sequences and target (reverse)
    std::reverse(target_seq.begin(), target_seq.end());

    ifstream in;  //open input file with sequences
    if (!in) {
        cout << "No input file!" << endl;
        exit(-1);
    }
    
    in.open(filename);
    ofstream out;
    out.open(outfile_name,ios::out);
    out << "MCO" << '\t' << "MCO_2nd" << '\t' << "TMO" << '\t' << "LTO" << '\n';
    
    //******************************** Loop over input file *******************************//
    int first_pos_def=0, second_pos_def=0;
    for (int i=0; i<tot_lines; ++i){
        in >> curr_strand;
        
        //convert string into vector
        std::copy(curr_strand.begin(), curr_strand.end(), std::back_inserter(sequence));

        first_pos_def=0;
        second_pos_def=0;
        vector <int> first_subseq_def(0);
        vector <int> second_subseq_def(0);
        compute_MCO_TMO_LTO(sequence,target_seq,MCO,TMO,LTO,first_pos_def,second_pos_def,first_subseq_def, second_subseq_def);
        //compute_MCO_2nd_solo must be executed after the previous function, since we need to know the MCO, to compute MCO_2nd
        compute_MCO_2nd_solo(sequence,target_seq,MCO_2nd,first_pos_def,second_pos_def,first_subseq_def, second_subseq_def,MCO);
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
bool match (char a, char b){
	if ((a=='A' and b=='T') or (a=='T' and b=='A') or (a=='C' and b=='G')  or (a=='G' and b=='C') )
		return true;
	else
		return false;
}

void compute_MCO_TMO_LTO (vector<char> first_seq, vector<char> second_seq, int& MCO, int& TMO, int& LTO, int& first_pos_def, int & second_pos_def, vector<int>& first_subseq_def, vector<int>& second_subseq_def) { //compute the value of the max consecutive overlap between two strands

    int actual_MCO = 0, actual_TMO=0, actual_LTO=0; //it will be Max overlaps
    int overlap, tot_mix_overlap, largest_tot_overlap; //counters for overlap (consecutive) and total overlap (non consecutive)
    int L=first_seq.size();
    int l=second_seq.size();
    int limit = L+l-1; //rightmost possible coordinate of a nucleotide
    int p,r;
    vector <int> first_subseq(0);
    vector <int> second_subseq(0);
    int first_pos=0, second_pos=0;

    for(int j=0;j<limit;j++) { //loop over all the possible positions
        overlap=0; //setting to zero the overlap counter
        tot_mix_overlap=0;
        largest_tot_overlap=0;
        //cycling on predator and resource genes that are overlapping
        //with two variables.
        first_subseq.clear();
        second_subseq.clear();
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
            if (p<0)
                continue; //safe check
            else {
                if ( match(first_seq[p],second_seq[r]) ){ //if the base pair matches
                    overlap++;
                    tot_mix_overlap++;
                    largest_tot_overlap++;
                    first_subseq.push_back(first_seq[p]);
                    second_subseq.push_back(second_seq[r]);
                    if (overlap==1){
                        first_pos=p;
                        second_pos=r;
                    }
                }
                else {  //if there is a mismatch
                    if (overlap > actual_MCO){
                        actual_MCO = overlap;
                        first_subseq_def=first_subseq;
                        second_subseq_def=second_subseq;
                        first_pos_def=first_pos;
                        second_pos_def=second_pos;
                        //TMO
                        actual_TMO=tot_mix_overlap;
                    }
                    
                    overlap=0;
                    first_subseq.clear();
                    second_subseq.clear();
                    }
                }
            
        }//end for p,r
        
        if (overlap>actual_MCO){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            //MCO
            actual_MCO = overlap; //update max
            first_subseq_def=first_subseq;
            second_subseq_def=second_subseq;
            first_pos_def=first_pos;
            second_pos_def=second_pos;
            //TMO
            actual_TMO=tot_mix_overlap;
        }
        
        //check whether to update LTO
        if (largest_tot_overlap>actual_LTO){
            actual_LTO=largest_tot_overlap;
        }
        
    } //end for j

    MCO=actual_MCO;
    TMO=actual_TMO;
    LTO=actual_LTO;
    
    return;
}

void compute_MCO_2nd_solo (vector<char> first_seq, vector<char> second_seq, int& MCO_2nd, int first_pos_def, int second_pos_def, vector<int> first_subseq_def, vector<int> second_subseq_def, int MCO) { //compute the value of the max consecutive overlap between two strands

    int Max = 0; //it will be Max overlaps
    int overlap; //counter for overlap (2nd MCO)
    int L=first_seq.size();
    int l=second_seq.size();
    int limit = L+l-1; //rightmost possible coordinate of a nucleotide
    int p,r;
    vector <int> first_subseq(0);
    vector <int> second_subseq(0);
    int first_pos=0, second_pos=0;

    for(int j=0;j<limit;j++) { //loop over all the possible positions
        overlap=0; //setting to zero the overlap counter
        //cycling on predator and resource genes that are overlapping
        //with two variables
        first_subseq.clear();
        second_subseq.clear();
        for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
            if (p<0)
                continue; //safe check
            else {
                if ( match(first_seq[p],second_seq[r]) ){ //if the base pair matches
                    overlap++;
                    first_subseq.push_back(first_seq[p]);
                    second_subseq.push_back(second_seq[r]);
                    if (overlap==1){
                        first_pos=p;
                        second_pos=r;
                    }
                }
                else {  //if there is a mismatch
                    if (overlap< MCO and overlap > Max and (first_pos!=first_pos_def or second_pos!=second_pos_def) ){
                        Max = overlap;
                    }
                    
                    overlap=0;
                    first_subseq.clear();
                    second_subseq.clear();
                    }
                }
        }//end for p,r
        
        if (overlap< MCO and overlap>Max and (first_pos!=first_pos_def or second_pos!=second_pos_def) ){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            Max = overlap; //update max
        }

    } //end for j

    MCO_2nd=Max;
    
    return;
}


//************************* not used *************************

void compute_MCO_TMO (vector<char> first_seq, vector<char> second_seq, int& MCO, int& MCO_2nd, int& MCO_3rd, int& TMO) { //compute the value of the max consecutive overlap between two strands

int Max = 0, Max_2nd = 0, Max_3rd=0, actual_TMO=0; //it will be Max overlaps
int overlap, tot_mix_overlap; //counters for overlap (consecutive) and total overlap (non consecutive)
int L=first_seq.size();
int l=second_seq.size();
int limit = L+l-1; //rightmost possible coordinate of a nucleotide
int p,r;
vector <int> first_subseq(0);
vector <int> second_subseq(0);
vector <int> first_subseq_def(0);
vector <int> second_subseq_def(0);
   int first_pos=0, second_pos=0;
int first_pos_def=0, second_pos_def=0;

for(int j=0;j<limit;j++) { //loop over all the possible positions
    overlap=0; //setting to zero the overlap counter
    tot_mix_overlap=0;
    //cycling on predator and resource genes that are overlapping
    //with two variables.
    first_subseq.clear();
    second_subseq.clear();
    for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
        if (p<0)
            continue; //safe check
        else {
            if ( match(first_seq[p],second_seq[r]) ){ //if the base pair matches
                overlap++;
                tot_mix_overlap++;
                first_subseq.push_back(first_seq[p]);
                second_subseq.push_back(second_seq[r]);
                if (overlap==1){
                    first_pos=p;
                    second_pos=r;
                }
            }
            else {  //if there is a mismatch
                if (overlap > Max){
                    if (Max>Max_2nd)
                        Max_2nd=Max; //update MCO_2nd
                    Max = overlap;
                    first_subseq_def=first_subseq;
                    second_subseq_def=second_subseq;
                    first_pos_def=first_pos;
                    second_pos_def=second_pos;
                    //TMO
                    actual_TMO=tot_mix_overlap;
                }
                //check 2nd largest MCO
                if (overlap<Max and overlap >= Max_2nd){
                    if (Max_2nd>Max_3rd)
                        Max_3rd=Max; //update MCO_3rd
                    Max_2nd = overlap; //update max
                }
                
                if (overlap<Max_2nd and overlap >= Max_3rd){
                    Max_3rd = overlap; //update max
                }
                
                overlap=0;
                first_subseq.clear();
                second_subseq.clear();
                }
            }
        //cout << p << "   " << overlap << endl;
    }//end for p,r
    
    if (overlap>Max){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
        //MCO
        if (Max>Max_2nd)
            Max_2nd=Max; //update MCO_2nd
        Max = overlap; //update max
        first_subseq_def=first_subseq;
        second_subseq_def=second_subseq;
        first_pos_def=first_pos;
        second_pos_def=second_pos;
        
        //TMO
        actual_TMO=tot_mix_overlap;
    }
    
    //check 2nd largest MCO
    if (overlap<Max and overlap >= Max_2nd){
        if (Max_2nd>Max_3rd)
            Max_3rd=Max; //update MCO_3rd
        Max_2nd = overlap; //update max
    }
    
    if (overlap<Max_2nd and overlap > Max_3rd){
        Max_3rd = overlap; //update max
    }

} //end for j

MCO=Max;
MCO_2nd=Max_2nd;
MCO_3rd=Max_3rd;
TMO=actual_TMO;

return;
}

