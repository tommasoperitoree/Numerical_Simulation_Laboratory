#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// Custom random number generator header
#include "random.h"

using namespace std;

// Function to calculate statistical uncertainty
double Error(double average, double squared, int n) {
    if (n == 0) {
        return 0.0; // If there are no samples, return zero uncertainty
    } else {
        // Calculate the standard error
        return sqrt(fabs(squared - pow(average, 2)) / n);
    }
}

int main(int argc, char *argv[]) {
    // Parameters for Buffon's Needle experiment
    double l = 0.5; // Needle length
    double d = 2.0; // Distance between grid lines

    // Setup for block averaging
    int iterations = pow(10, 6); // Total number of needle throws
    int blocks = 100; // Number of blocks
    int blockSize = iterations / blocks; // Size of each block

    // Variables for tracking the progress and results
    double pi, pi2, progressivePi = 0, progressivePiSquared = 0, progressiveError;

    // Buffon's Needle experiment variables
    double end1, x, y, sin_theta, end2;
    int N_hits; // Counter for hits

    // Initialize random number generator
    Random rnd;
    rnd.RandomImplementation();

    // Create an output file to store results
    ofstream outFile("1.3.dat");
    outFile << "# BLOCKS\t\t\tAVE :\t\tERROR :" << endl;

    // File to check random theta distribution
    ofstream theta("theta_check.dat");

    // Loop through each block
    for (int i = 0; i < blocks; i++) {
        N_hits = 0; // Reset the hit counter for each block
        
        // Loop through each throw in the block
        for (int j = 0; j < blockSize; j++) {
            // Generate the position of one end of the needle
            end1 = rnd.Rannyu(0., d);
            
            // Generate a random angle (theta) for the needle
            do {
                x = rnd.Rannyu(); 
                y = rnd.Rannyu(-1., 1.); // Generate a random point in a circle
            } while (sqrt(pow(x, 2) + pow(y, 2)) > 1); // Ensure it's within the unit circle
            
            sin_theta = y / sqrt(pow(x, 2) + pow(y, 2)); // Calculate sin(theta)
            end2 = end1 + l * sin_theta; // Calculate the other end of the needle
            
            // Check if the needle crosses a line
            if (end2 >= d || end2 <= 0) {
                N_hits++; // Increment the hit counter
            }
        }

        // Calculate the estimated value of pi
        pi = 2 * l * blockSize / ((double)N_hits * d); // Formula for Buffon's Needle experiment
        pi2 = pow(pi, 2); // Square of the estimated pi
        
        // Update the progressive averages
        progressivePi += pi;
        progressivePiSquared += pi2;
        
        // Calculate the progressive error
        progressiveError = Error(progressivePi / (i + 1), progressivePiSquared / (i + 1), i);
        
        // Write the block results to the output file
        outFile << setw(10) << i 
                << setw(16) << progressivePi / (i + 1) 
                << setw(16) << progressiveError 
                << endl;
    }

    // Close the output file
    outFile.close();

    // Save the random number generator seed
    rnd.SaveSeed();

    return 0;
}
