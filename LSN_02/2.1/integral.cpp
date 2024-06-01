#include "integral.h"

using namespace std;

Integral::Integral (double a, double b, Functions* f, Random* rnd){

    /*
    Integral generic class Constructor
    */

    I_min = min(a, b);                                                                  // Define the domain limits
    I_max = max(a, b);
    I_f = f;                                                                            // Data member representing the function
    I_rnd = rnd;                                                                        // Data member for the random seed
    if (a>b) I_sign = -1.;      
    else I_sign = 1.;                                                                   // Integral sign data member
    I_integral = 0;                                                                     // Value of the integral
}

double Integral::arithAverage(unsigned int N){

   /*
   Uniform sampling technique.
   Input unsigned int N: size of the partition used for the integral evaluation
   */

   I_integral = 0;                                                                     // Reinitialize to zero the value

   for (int i = 0; i < N; i++){

      I_integral += I_f -> Evaluate(I_rnd -> Rannyu(I_min, I_max));                   // Evaluate the function on a point within the range

   }

   I_integral = I_sign * (I_max - I_min) * I_integral / (double) N;                    // compute average, add the sign

   return I_integral;
}

double Integral::impSampling(unsigned int N, Functions *p, Functions *inv_p){

    /*
    Importance sampling technique.
    Inputs:
    - unsigned int N : size of the partition used for the integral evaluation
    - Functions *p : pointer to the function used for importance sampling
    - Functions *inv_p : 
    */

    I_integral = 0;                                                                // Reinitialize to zero the value

    for (int i = 0; i < N; i++){

        double x = inv_p->Evaluate(I_rnd->Rannyu());                       // draw a number distributed as p(x) using the inverse of p                     
        I_integral += I_f->Evaluate(x) / p->Evaluate(x) ;                  // Evaluate the integral

    }

    I_integral = I_sign * (I_max - I_min) * I_integral / (double)N;                     // compute the mean and add the sign

    return I_integral;
}