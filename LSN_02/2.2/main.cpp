#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "random.h"  // Custom random number generator
#include "lib.h"     // Custom library for additional functions

using namespace std;

int main() {

	Random rnd;  // Create a random number generator instance
	rnd.RandomImplementation();  // Initialize the RNG with a predefined method

	int iterations = pow(10, 4);  // Total number of iterations for random walks
	int blocks = 100;  // Number of blocks to divide iterations into
	int maxSteps = 100;  // Maximum steps for the random walks
	int blockSize = iterations / blocks;  // Number of iterations per block

	// Initialize a 2D vector to track 3D positions for discrete random walks
	vector <vector <double> > positionsDiscrete(iterations, vector<double>(3, 0.0));  // Start at origin

	// Variables for progressive averages and errors
	double progressiveSum = 0.0;
	double progressiveSumSquared = 0.0;
	double progressiveError = 0.0;

	double x;

	/*
	------------------------------- 2.2.1 -------------------------------
	
								Random Walk on a cubic lattice
	
	*/

	double lattice_const = 1.0;  // Lattice step length for discrete random walks

	// Output file for results and another for positions
	ofstream outputFile("2.2.1.dat");
	outputFile << "# STEP: \t\t AVERAGE: \t\t ERROR:" << endl;

	string positionsDiscreteFile = "2.2.1_positions.dat";  // File to store positions
	ofstream outPositionDiscrete(positionsDiscreteFile);
	outPositionDiscrete << "# STEP \t\t\t x \t\t\t y \t\t\t z" << endl;
	outPositionDiscrete.close();

	// Iterate over each step from 0 to maxSteps
	for (int n = 0; n < maxSteps; n++) {
		progressiveSum = 0.0;
		progressiveSumSquared = 0.0;

		// Divide iterations into blocks and calculate averages
		for (int b = 0; b < blocks; b++) {
			x = 0.0;

			// Loop over block size for each block
			for (int j = 0; j < blockSize; j++) {
					// Perform a step in a discrete random walk
					x += latticeWalk(positionsDiscrete[b * blockSize + j], rnd, lattice_const, 1);
					// Save the position for the first iteration of the block
					if (b * blockSize + j == 0) {
						savePosition(positionsDiscrete[b * blockSize + j], positionsDiscreteFile, n);
					}
			}

			// Compute progressive averages and errors
			progressiveSum += sqrt(x / blockSize);  // Average distance for each block
			progressiveSumSquared += x / blockSize;  // Squared average distance
		}

		// Calculate the progressive error for the current step
		progressiveError = Error(progressiveSum / blocks, progressiveSumSquared / blocks, blocks);

		// Write the results to the output file
		outputFile << setw(8) << n + 1
						<< setw(16) << progressiveSum / blocks
						<< setw(16) << progressiveError << endl;
	}

	outputFile.close();

	/*
	------------------------------- 2.2.2 -------------------------------
	
							Random Walk on the continuum
	
	*/

	double step_length = 1.0;  // Step length for continuous random walks
	progressiveSum = 0.0;
	progressiveSumSquared = 0.0;
	progressiveError = 0.0;

	// Initialize a 2D vector to track 3D positions for continuous random walks
	vector<vector<double>> positionsContinuum(iterations, vector<double>(3, 0.0));

	// Output file for results and another for positions
	outputFile.open("2.2.2.dat");
	outputFile << "# STEP: \t\t AVERAGE: \t\t ERROR:" << endl;

	string positionsContinuumFile = "2.2.2_positions.dat";  // File to store positions
	ofstream outPositionContinuum(positionsContinuumFile);
	outPositionContinuum << "# STEP \t\t\t x \t\t\t y \t\t\t z" << endl;
	outPositionContinuum.close();

	// Iterate over each step from 0 to maxSteps for continuous random walks
	for (int n = 0; n < maxSteps; n++) {
		progressiveSum = 0.0;
		progressiveSumSquared = 0.0;

		// Iterate over blocks for continuous random walks
		for (int b = 0; b < blocks; b++) {
			double x = 0.0;

			// Loop over block size for each block
			for (int j = 0; j < blockSize; j++) {
					// Perform a step in a continuous random walk
					x += continuumWalk(positionsContinuum[b * blockSize + j], rnd, step_length, 1);
					// Save the position for the first iteration of the block
					if (b * blockSize + j == 0) {
						savePosition(positionsContinuum[b * blockSize + j], positionsContinuumFile, n);
					}
			}

			// Compute progressive averages and errors
			progressiveSum += sqrt(x / blockSize);  // Average distance for each block
			progressiveSumSquared += x / blockSize;  // Squared average distance
		}

		// Calculate the progressive error for the current step
		progressiveError = Error(progressiveSum / blocks, progressiveSumSquared / blocks, blocks);

		// Write the results to the output file
		outputFile << setw(8) << n + 1
						<< setw(16) << progressiveSum / blocks
						<< setw(16) << progressiveError << endl;
	}

	// Additional steps for discrete and continuous random walks position 
	for (int i = maxSteps + 1; i < 5000; i++) {
		x += latticeWalk(positionsDiscrete[0], rnd, lattice_const, 1);
		savePosition(positionsDiscrete[0], positionsDiscreteFile, i);

		x += continuumWalk(positionsContinuum[0], rnd, step_length, 1);
		savePosition(positionsContinuum[0], positionsContinuumFile, i);
	}

	// Close the output file and save the RNG seed
	outputFile.close();

	rnd.SaveSeed();  // Save the random number generator seed

	return 0;
}
