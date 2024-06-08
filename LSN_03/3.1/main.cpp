#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "random.h"

using namespace std ;

double Error (double average, double squared, int n) ;      // Function for statistical uncertainty estimation


int main () {

   Random rnd ;
   rnd.RandomImplementation();                                          // Initialization of Random generator

   int iterations = pow(10,6) ;                                         // total number of random numbers to be drawn
   int blocks = 100 ;                                                  	// number of blocks
	int blockSize = (int)(iterations/blocks) ;                           // Number of throws in each block

   cout << "Working with " << blocks <<" blocks, with "
		  << blockSize << " throws in each block" << endl ; 

   // Define variables of the pricing process
   double r = 0.1;                                                      // risk-free interest rate                                       
   double T = 1.;                                                       // delivery time
   double sigma = 0.25;                                                 // volatility
   double K = 100;                                                      // strike price
   double S_0 = 100;                                                    // initial asset price
    

   /*
   --------------------------------------------- 3.1.1 ----------------------------------------------
   
                           Direct Sampling with geometric Browmnian Motion GBM
   
   */

  	// Open output files
  	ofstream callFile("Data/call_direct.dat");
	ofstream putFile("Data/put_direct.dat");

	// Header files
	callFile << "# BLOCK:		CALL AVERAGE:			ERROR:" << endl ;	
	putFile << "# BLOCK:		PUT AVERAGE:			ERROR:" << endl ;

	// Define all the auxiliary variables needed for the computation
   double rand = 0.;
   double S_T = 0;
   double call = 0., put = 0.;

   double progressiveSumCall = 0., progressiveSquaredCall = 0.;
	double progressiveSumPut = 0., progressiveSquaredPut = 0.;
	double errorCall = 0., errorPut = 0.;


	for (int i = 0; i < blocks; i++){                                           	// Iterate over the blocks

		call = 0. ;
		put = 0. ;

		for (int j = 0; j < blockSize; j++){                                     	// Iterate over the throws in each block
			rand = rnd.Gauss(0, T) ;																// generate random with normal distribution N(0,T)
			S_T = S_0 * exp( (r - 0.5*pow(sigma,2)) * T + (sigma * rand) ) ;        // Compute the value of the asset at maturity using B&S model
			call += std::max(0., S_T - K) * exp (- r * T) ;                       	// Compute call price accordingly
			put += std::max(0., K - S_T)* exp (- r * T) ;                         	// Compute put price accordingly
		}

		call /= blockSize ;                                                      // Average over the number of throws                                  
		put /= blockSize ;

		// Call - block averaging
      progressiveSumCall = ((progressiveSumCall * i) + call) / (i + 1);          
		progressiveSquaredCall = ((progressiveSquaredCall * i) + pow(call, 2)) / (i + 1);
		errorCall = Error(progressiveSumCall, progressiveSquaredCall, i);
		callFile << setw(8) << i 
				  << setw(16) << progressiveSumCall 
				  << setw(16) << errorCall << endl;

      // Put - block averaging
      progressiveSumPut = ((progressiveSumPut * i) + put) / (i + 1);             
		progressiveSquaredPut = ((progressiveSquaredPut * i) + pow(put, 2)) / (i + 1);
		errorPut = Error(progressiveSumPut, progressiveSquaredPut, i);
		putFile << setw(8) << i 
				  << setw(16) << progressiveSumPut 
				  << setw(16) << errorPut << endl;

	}

	callFile.close() ;
	putFile.close() ;

   /*
   --------------------------------------------- 3.1.2 ----------------------------------------------
   
                           Indirect Sampling with geometric Browmnian Motion GBM
   
   */
   
	// Open output files
  	callFile.open("Data/call_indirect.dat");
	putFile.open("Data/put_indirect.dat");

	// Header files
	callFile << "# BLOCK:		CALL AVERAGE:			ERROR:" << endl ;	
	putFile << "# BLOCK:		PUT AVERAGE:			ERROR:" << endl ;

	// Reset all the auxiliary variables needed for the computation
   rand = 0.; S_T = 0;
   call = 0.; put = 0.;

   progressiveSumCall = 0.; progressiveSquaredCall = 0.;
	progressiveSumPut = 0.; progressiveSquaredPut = 0.;
	errorCall = 0.; errorPut = 0.;

	// Define new variables needed to quantize the time interval
	int steps = 100;
	double increment = T / (double) steps ;

	for (int i = 0; i < blocks; i++){                                           	// Iterate over the blocks

		call = 0. ;
		put = 0. ;

		for (int j = 0; j < blockSize; j++){                                     	// Iterate over the throws in each block

			S_T = S_0 ;
			
			for (int k = 0; k < steps; k++){

				rand = rnd.Gauss(0,1);
				S_T *= exp((r - 0.5*pow(sigma,2)) * increment + (sigma * rand * sqrt(increment)));     

			}

			call += std::max(0., S_T - K) * exp (- r * T) ;                       	// Compute call price accordingly
			put += std::max(0., K - S_T)* exp (- r * T) ;                         	// Compute put price accordingly
		}

		call /= blockSize ;                                                    // Average over the number of throws                                  
		put /= blockSize ;

		// Call - block averaging
      progressiveSumCall = ((progressiveSumCall * i) + call) / (i + 1);          
		progressiveSquaredCall = ((progressiveSquaredCall * i) + pow(call, 2)) / (i + 1);
		errorCall = Error(progressiveSumCall, progressiveSquaredCall, i);
		callFile << setw(8) << i 
				  << setw(16) << progressiveSumCall 
				  << setw(16) << errorCall << endl;

      // Put - block averaging
      progressiveSumPut = ((progressiveSumPut * i) + put) / (i + 1);             
		progressiveSquaredPut = ((progressiveSquaredPut * i) + pow(put, 2)) / (i + 1);
		errorPut = Error(progressiveSumPut, progressiveSquaredPut, i);
		putFile << setw(8) << i 
				  << setw(16) << progressiveSumPut 
				  << setw(16) << errorPut << endl;

	}


   rnd.SaveSeed();


   return 0;
}

double Error (double average, double squared, int n) {     // Function for statistical uncertainty estimation

   if (n==0) {
      return 0.;
   }
   else 
      return sqrt (fabs (squared - pow(average, 2) ) / n);
}