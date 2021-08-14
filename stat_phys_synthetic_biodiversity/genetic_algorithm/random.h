//
//  random.h
//  
//
//  Created by Francesco Mambretti on 10/07/21.
//

#ifndef __random__
#define __random__

#include <fstream>
#include <cstdlib>
#include <iostream>
#include "global.h"

class Random {

private:
    int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

public:
    // constructors
    Random();
    // destructor
    ~Random();
    // methods
    void SetRandom(int * , int, int);
    void SaveSeed();
    double Rannyu(void);
    double Rannyu(double min, double max);
    double Gauss(double mean, double sigma);
};

#endif // __random__
