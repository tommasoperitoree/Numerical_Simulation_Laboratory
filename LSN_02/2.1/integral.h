#ifndef __Integral__
#define __Integral__

#include <cmath>
#include <iostream>
#include "functions.h"
#include "random.h"

using namespace std;

class Integral {

    /*
    This class contains all the methods used for computing the integral of a `Functions` object
    */
    
    public:
        Integral(double a, double b, Functions* f, Random* rnd);                        // Constructor
        double arithAverage (unsigned int N);                                           // Uniform sampling
        double impSampling (unsigned int N, Functions* p, Functions* inverse);               // Importance sampling

    private:
        double I_min;                                                                   // integral domain
        double I_max;
        double I_integral;                                                              // used to store the integral value
        double I_sign;                                                                  // '' the ingral sign
        Functions* I_f;                                                                 
        Random* I_rnd;
};


#endif //__Integral__