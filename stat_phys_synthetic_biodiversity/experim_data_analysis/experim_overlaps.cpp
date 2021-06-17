//
//  experim_overlaps.cpp
//
//  Created by Francesco Mambretti on 07/06/21.
//
//  Read a list of unknown length of (experimental) 50-long sequences of A,T,C,G.
//  For each of them, check their largest maximum consecutive overlap with the target strand, a_max. We compare here the predators with the reversed and complemented input target, so each base pair is regarded as 'matching' if the bases are identical (while in the real world A-T and C-G are the matching bases).
//  Once found, compute also the tot overlap with the target, in the position corresponding to the (last recorded) largest max consecutive overlap for each of them.
//  At the end, generates two files which can be processed by two external scripts:
//  1) histo.py for the relative species abundance histogram
//  2) compare.py which makes a 2D histogram of the largest a_max and of a_tot 
//
//  17/06/2021 - OK version

#include "experim_overlaps.hpp"

using namespace std;

int main(int argc, char** argv){
    
    if (argc!=4){
        cerr << "Error! Correct usage: " << argv[0] << " target_filepath    predators_filepath      output_folder" << endl;
        exit(EXIT_FAILURE);
    }
    
    else{
        cout << "The input target strand is reverted and complemented. Checking the overlap of each 50mer with it..." << endl;
    }
    
    int l=20; //target length
    int L=50; //predator length
    int a,b; //for loops indexes
    
    vector<int> target (l);  //allocate the target
    // map DNA letters onto numbers
    std::map <char, int> alphabet;
    alphabet['A'] = 0;
    alphabet['T'] = 1;
    alphabet['C'] = 2;
    alphabet['G'] = 3;
    
    char char_temp;
    
    ifstream in(argv[1]);
    
    for (a=0; a < l; ++a){
        in >> char_temp;
        target[a]=alphabet[char_temp];
        cout << a << '\t' << target[a] << endl;
    }
       
    in.close();
    
    target=rev_and_comp(target);
    
    vector<int> predator (L);
    
    //read population of strands of arbitrary length and check the largest max consecutive overlap + the total overlap in the position corresponding to the (last recorded) max consecutive overlap for each of them 
    int appo, appo2, pos;
    
    in.open(argv[2]);
    string folder=argv[3];
    ofstream outFile(folder+"/all_max_cons_overlaps.dat");
    ofstream outFile2(folder+"/all_tot_overlaps.dat");
    
    while (in >> char_temp) {
        predator[0]=alphabet[char_temp];

        for (b=1; b < L; ++b){
            in >> char_temp;
            predator[b]=alphabet[char_temp];
        }
        
        appo=check_max_cons_overlap(target,predator,l,L,pos);
        appo2=check_tot_overlap(target,predator,l,L,pos);
        outFile << appo << endl;
        outFile2 << appo2 << endl;
    }
    
    in.close();
    
    outFile.close();
    outFile2.close();
    
    return 0;
}


// FUNCTIONS ==================================================================================

//======================================
// check_max_cons_overlap
//======================================

int check_max_cons_overlap(vector<int> target, vector<int> predator, int l, int L, int & pos){
    
    int my_max = 0, temp_max;
    int p,q; //counter
    vector<int> overlaps (l,0);
    vector<int> overlaps_q (l,0);
    int i,j;
    
    // overlap with resource completely over predator
    for (i=0;i <= L-l;++i) {
        fill(overlaps.begin(), overlaps.end(), 0);
        p=0; //it counts consecutive overlaps
        for (j=i;j<l+i;++j) { // internal target index, always take l values
            if (match (target[j-i], predator[j])) //they match
                overlaps[p]++;
            else p++;
        }
        temp_max=*max_element(overlaps.begin(), overlaps.end());
        if (temp_max >= my_max){
            my_max = temp_max;
            pos = i; //save the position of the target left end
        }
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
        if (temp_max >= my_max){
            my_max = temp_max;
            pos = -(l-1-i); //save the position of the target left end
        }
        
        temp_max=*max_element(overlaps_q.begin(), overlaps_q.end());
        if (temp_max >= my_max){
            my_max = temp_max;
            pos = L-1-i;
        }
     }
    
    return my_max;
}

//======================================
// check_tot_overlap
//======================================

int check_tot_overlap(vector<int> target, vector<int> predator, int l, int L, int pos){
    
    int my_max = 0, tot=0, j, start;
    int limit_loop;
    // we only want to know the total overlap in correspondence of the given relative position 'pos' between target and predator

    if (0<= pos and pos <=L-l){
        start=pos;
        limit_loop=pos+l;
    }

    else if (pos < 0){
        start=0;
        limit_loop=l+pos;
    }

    else if (pos > L-l){
        start=pos;
        limit_loop=L;
    }


    tot=0;
    for (j=start; j < limit_loop; ++j){
        if (match (target[j-pos], predator[j])) //they match
            tot++;
    }

    return tot;
}

//=============
// match
//=============

bool match (int x, int y){
    
    //here, the method for pairing A-T and C-G is commented
    //if the string is reversed and complemented, a simple x==y is needed
    /*if (x<2){ // x=0 or 1
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
    }*/
    if (x==y)
        return true;
    else
        return false;
}

//=================
// rev_and_comp
//=================

vector<int> rev_and_comp(vector<int>& strand){
    reverse(strand.begin(),strand.end());
    int size=strand.size();
    for (int i=0; i < size; ++i) {
        if (strand[i]%2==0)
            strand[i]+=1;
        else
            strand[i]-=1;
    }
    return strand;
};
