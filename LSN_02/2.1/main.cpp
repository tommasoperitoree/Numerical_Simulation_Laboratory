#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// Include custom random number generator and function libraries
#include "random.h"
#include "functions.h"
#include "integral.h"
#include "lib.h"

using namespace std;

int main() {

	// Initialize the random number generator
	Random* rnd = new Random();
	rnd->RandomImplementation();

	// Number of random numbers to be drawn
	int iterations = pow(10, 5); // 100,000 random numbers
	int blocks = 100; // Number of blocks
	int blockSize = iterations / blocks; // Number of throws in each block

	// Variables for progressive statistics
	double x = 0.0;
	double progressiveSum = 0.0;
	double progressiveSumSquared = 0.0;
	double progressiveError = 0.0;

	/* ------------------------------- 2.1.1 -------------------------------
   
                  Sampling with uniform distribution in [0,1)

   */

	// Open an output file for results
	ofstream outputFile("2.1.1.dat");

	outputFile << "# BLOCK: \t\t AVERAGE: \t\t ERROR:" << endl ;

	// Initialize an integrand and integral object
	Cosine* cosine = new Cosine(M_PI / 2., M_PI / 2., 0.);
	Integral* integral = new Integral(0, 1, cosine, rnd);

	// Iterate over blocks blocks to estimate the integral
	for (int i = 0; i < blocks; i++) {
		// Calculate the integral using uniform sampling
		x = integral->arithAverage(blockSize);

		// Calculate progressive sums for averaging and squared values
		progressiveSum = ((progressiveSum * i) + x) / (i + 1);
		progressiveSumSquared = ((progressiveSumSquared * i) + pow(x, 2)) / (i + 1);

		// Calculate the error
		progressiveError = Error(progressiveSum, progressiveSumSquared, i);

		// Write the results to the output file
		outputFile << setw(8) << i 
				  << setw(16) << progressiveSum
				  << setw(16) << progressiveError << endl;
	}

	outputFile.close(); // Close the output file



	/* ------------------------------- 2.1.2 -------------------------------

					Sampling with first-order Taylor approximation in 1

	*/

	// Create pointers for different function classes
	Parabola* p = new Parabola(0, -2, 2);
	SquareRoot* inv_p = new SquareRoot(1, -1, 1);

	// Open a new output file for results
	outputFile.open("2.1.2.dat");

	outputFile << "# BLOCK: \t\t AVERAGE: \t\t ERROR:" << endl ;

	// Reset progressive variables
	x = 0.0;
	progressiveSum = 0.0;
	progressiveSumSquared = 0.0;
	progressiveError = 0.0;

	// Iterate over blocks blocks with non-uniform sampling
	for (int i = 0; i < blocks; i++) {
		// Calculate the integral using importance sampling
		x = integral->impSampling(blockSize, p, inv_p);

		// Calculate progressive sums for averaging and squared values
		progressiveSum = ((progressiveSum * i) + x) / (i + 1);
		progressiveSumSquared = ((progressiveSumSquared * i) + pow(x, 2)) / (i + 1);

		// Calculate the error
		progressiveError = Error(progressiveSum, progressiveSumSquared, i);

		// Write the results to the output file
		outputFile << setw(8) << i 
				  << setw(16) << progressiveSum
				  << setw(16) << progressiveError << endl;
	}

	outputFile.close(); // Close the output file


	// Save the random number generator seed
	rnd->SaveSeed();

	return 0;
}
