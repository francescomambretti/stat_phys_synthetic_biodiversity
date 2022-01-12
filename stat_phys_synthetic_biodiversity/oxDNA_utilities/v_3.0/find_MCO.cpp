// Written by Francesco Mambretti, 20/09/2021
// Finds the MCO between two strands and returns the subsequences involved
// Note: if there are two overlapping regions with the same MCO, the code only finds one of them

// Output: the MCO with the sequence file as it is and the MCO with the second sequence reverted

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

//functions
int compute_MCO(vector<char>,vector<char>);
bool match(char,char);

int main(int argc, char** argv){

	if (argc!=2)
	{
		cerr << "Please provide a file with the two sequences to be compared as an argument of the executable!" << endl;
		exit(EXIT_FAILURE);
	}	

	string filename=argv[1];
	ifstream in;
	in.open(filename);
	
	string first;
	string second;

	in >> first;
	in >> second;

	//cout<< first << endl;
	//cout << second << endl;

	in.close();

	//convert strings into vectors
	std::vector<char> first_seq;
    	std::copy(first.begin(), first.end(), std::back_inserter(first_seq)); 

    	for (const char &c: first_seq) {
    	    cout << c;
    	}
	cout<<endl;

	std::vector<char> second_seq;
    	std::copy(second.begin(), second.end(), std::back_inserter(second_seq)); 

    	for (const char &c: second_seq) {
    	    cout << c;
    	}
	cout<<endl;

	//check max cons overlap between first and second

	int MCO;// = compute_MCO(first_seq,second_seq);

	//cout << MCO << endl;

	std::reverse(second_seq.begin(), second_seq.end());

	MCO = compute_MCO(first_seq,second_seq);

	cout << MCO << endl;

	return 0;
}

//***********************************************************//
// FUNCTIONS 

int compute_MCO (vector<char> first_seq, vector<char> second_seq) { //compute the value of the max consecutive overlap between two strands

    	int Max = 0; //it will be Max overlaps
    	int overlap; //counter for overlap
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
        	//cycling on predator and resource genes that are overlapping
        	//with two variables.
		first_subseq.clear();
		second_subseq.clear();
        	for (p=L-1-j, r=0; p<L && r<l; p++,r++) { //p: track the left end of the resource in the first sequence reference frame (fixed); r: track the internal coordinate of the second sequence
			if (p<0)
				continue; //safe check
			else {
                		if ( match(first_seq[p],second_seq[r]) ){
                    			overlap++;
					first_subseq.push_back(first_seq[p]);
					second_subseq.push_back(second_seq[r]);
					if (overlap==1){
						first_pos=p;
						second_pos=r;
					}
				}
                		else {
                   			if (overlap > Max){
                        			Max = overlap;
						first_subseq_def=first_subseq;
						second_subseq_def=second_subseq;
						first_pos_def=first_pos;
						second_pos_def=second_pos;
					}
                    			overlap=0;
					first_subseq.clear();
					second_subseq.clear();
                		}
            		}
			//cout << p << "   " << overlap << endl;
		}//end for p,r
        
        	if (overlap>Max){ // useful when there is a group of matching bases at the last loop, which is not ended by any non-matching pair
            		Max = overlap; //update max
			first_subseq_def=first_subseq;
			second_subseq_def=second_subseq;
			first_pos_def=first_pos;
			second_pos_def=second_pos;
		}

	} //end for j

	cout << "Overlapping sequences are: " << endl;
	for (const char &c: first_subseq_def) 
		cout << c;
	cout << endl;

 	for (const char &d: second_subseq_def)
		cout << d;
	cout << endl;

	cout << "... and they start at position " << first_pos_def << " and " << second_pos_def << endl; 

	return Max;
}

//decide pairing rules
bool match (char a, char b){
	if ((a=='A' and b=='T') or (a=='T' and b=='A') or (a=='C' and b=='G')  or (a=='G' and b=='C') )
		return true;
	else
		return false;
}
