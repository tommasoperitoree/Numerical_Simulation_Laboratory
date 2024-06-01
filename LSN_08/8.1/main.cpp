#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "random.h"
#include "lib.h"

using namespace std ;


int main () {

   Random rand ;
   rand.RandomImplementation();                                          // Initialization of Random generator

	// set starting point
	double x = 0.;

	// Define parameters of simulation
	double mu = 0.0;
	double sigma = 1.0;

	// monte carlo variables
	double integral = 0.;
	double delta = 2.;

	// data blocking variables
   int iterations = pow(10,6) ;                                         // total number of random numbers to be drawn
   int blocks = 100 ;                                                  	// number of blocks
	int blockSize = (int)(iterations/blocks) ;                           // Number of throws in each block


	// Open output files
	string dirHamFile = "./HamiltonianExpValue_" + to_string(mu).substr(0,3) + "_" + to_string(sigma).substr(0,3) + ".dat" ;
  	ofstream hamiltonianFile(dirHamFile) ;
	string dirCoordFile = "./Coordinates_" + to_string(mu).substr(0,3) + "_" + to_string(sigma).substr(0,3) + ".dat" ;
	ofstream coordinatesFile(dirCoordFile) ;

	// Header files
	hamiltonianFile << "# BLOCK:		ACTUAL_HAM:	 	AVERAG_HAM:			ERROR:" << endl ;
	coordinatesFile << "# x: 		x: " << endl ;


	double accumulative = 0., accumulativeSquared = 0.;
   double progressiveSum = 0. ;
	double progressiveSumSquared= 0. ;
	double progressiveError = 0. ;

	int accepted = 0;
	int attempted = 0;

   Equilibrate(100, pow(10, 3), x, rand, delta, mu, sigma);


//Blocking average variables
   for(int j = 0; j < blocks; ++j){
      integral = 0;
      accepted = 0;
      attempted = 0;
      for (int i = 0; i < blockSize; ++i){
         Metropolis(x, rand, delta, accepted, attempted, mu, sigma);
         integral += evalHamiltonian(x,mu,sigma);
			coordinatesFile << setw(8) << j + i*blockSize 
							 << setw(16) << x << endl;
		}
        		
		accumulative += integral / double(blockSize);
		accumulativeSquared += pow(integral / double(blockSize), 2);
		
		progressiveSum = accumulative / double(j + 1);
		progressiveSumSquared = accumulativeSquared / double(j + 1);
		progressiveError = Error(progressiveSum, progressiveSumSquared, j);
        
		hamiltonianFile << setw(8) << j
				  << setw(16) << integral / double(blockSize) 
				  << setw(16) << progressiveSum
				  << setw(16) << progressiveError << endl;
   }

	
	hamiltonianFile.close() ;

	coordinatesFile.close() ;


   rand.SaveSeed();

   return 0;
}