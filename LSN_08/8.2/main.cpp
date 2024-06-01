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


	// define parameters of simulation
   double T = 2;   //temperature
   double beta = 1./T;
	double TMin = 0.01;

	double mu = 0.;
	double sigma = 1.;
	double muBest = mu, sigmaBest = sigma;
   
	// annealing
   double energy = 1.;
   double currentEnergy = energy;
	double bestEnergy = energy;

	double position = 0.;
	double integral = 0.;


	// monte carlo variables
	double delta = 2.;

	// data blocking variables
   int iterations = pow(10,6) ;                                         // total number of random numbers to be drawn
   int blocks = 1000 ;                                                  	// number of blocks
	int blockSize = int(iterations/blocks) ;                           // Number of throws in each block
	int step = 0;

	// Open output files
  	ofstream hamiltonianFile("./HamiltonianExpValue.dat") ;

	// Header files
	hamiltonianFile << "# STEP: \t\t\t\t T: \t\t\t\t MU: \t\t\t\t SIGMA: \t\t  ENERGY: \t\t ENERGY_BEST: \t\t ERROR:" << endl ;

	while (T >= TMin) {

		double muOld = mu;
     	double sigmaOld = sigma;

		// update parameters
		mu = abs(muOld + rand.Rannyu(-1.,1.) * 0.5 * 1./beta) ;
		sigma = abs(sigmaOld + rand.Rannyu(-1.,1.) * 0.25 * (1./beta)) ; 
		
		// trial parameters
		//mu = abs(muOld + rand.Rannyu(-1.,1.) * 1./beta) ;
		//sigma = abs(sigmaOld + rand.Rannyu(-1.,1.) * (1./beta)) ; 
		

		integral = 0;
		double accumulative = 0., accumulativeSquared = 0.;
   	double progressiveSum = 0. ;
		double progressiveSumSquared= 0. ;
		double progressiveError = 0. ;

		Equilibrate(100, pow(10, 3), position, rand, delta, mu, sigma);


      for(int j = 0; j < blocks; ++j){
      	
			integral = 0.;
      	int accepted = 0, attempted = 0;

      	for (int i = 0; i < blockSize; ++i){
      	   Metropolis(position, rand, delta, accepted, attempted, mu, sigma);
      	   integral += evalHamiltonian(position,mu,sigma);
			}
	
			accumulative += integral / double(blockSize);
			accumulativeSquared += pow(integral / double(blockSize), 2);

			progressiveSum = accumulative / double(j+1);
			progressiveSumSquared = accumulativeSquared / double(j+1);
			progressiveError = Error(progressiveSum, progressiveSumSquared, j);

		}

		currentEnergy = progressiveSum ;

		double acceptance = min(1., (exp(- beta * (currentEnergy - energy))));
	
		if (rand.Rannyu() < acceptance)
         energy = currentEnergy;
      else {
         mu = muOld;
         sigma = sigmaOld;
      }

		hamiltonianFile << setw(8) << step
			<< setw(16) << T
			<< setw(16) << mu
			<< setw(16) << sigma
			<< setw(16) << energy
			<< setw(16) << bestEnergy
			<< setw(16) << progressiveError << endl;
		
		step += 1;
		T *= 0.997 ;
      beta = 1. / T;
      if (bestEnergy > energy) {
         bestEnergy = energy;
         muBest = mu;
         sigmaBest = sigma;
		}

	}
	
	std::cout << "best parameters found :  mu_best = " << muBest << ", sigma_best = " << sigmaBest << endl; 
	std::cout << "with final energy : energy_best = " << bestEnergy << endl ;

	hamiltonianFile.close(); 


	// now 8.1 with best values of parameters

	// Open output files
  	hamiltonianFile.open("./GSHamiltonianExpValue.dat") ;
	ofstream coordinatesFile("./GSCoordinates.dat") ;

	// Header files
	hamiltonianFile << "# BLOCK:		ACTUAL_HAM:	 	AVERAG_HAM:			ERROR:" << endl ;	
	coordinatesFile << "# STEP: 		POSITION: " << endl ;


	double accumulative = 0., accumulativeSquared = 0.;
   double progressiveSum = 0. ;
	double progressiveSumSquared= 0. ;
	double progressiveError = 0. ;

	int accepted = 0;
	int attempted = 0;

	integral = 0.;
	position = 0.;

   Equilibrate(100, pow(10, 3), position, rand, delta, muBest, sigmaBest);


//Blocking average variables
   for(int j = 0; j < blocks; ++j){
      integral = 0;
      accepted = 0;
      attempted = 0;
      for (int i = 0; i < blockSize; ++i){
         Metropolis(position, rand, delta, accepted, attempted, muBest, sigmaBest);
         integral += evalHamiltonian(position,muBest,sigmaBest);
			coordinatesFile << setw(8) << j + i*blockSize 
							 << setw(16) << position << endl;
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
	
	coordinatesFile.close() ;
	hamiltonianFile.close(); 

   rand.SaveSeed();

   return 0;
}