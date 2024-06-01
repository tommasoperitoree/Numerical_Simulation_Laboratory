#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// Include the custom random number generator
#include "random.h"

using namespace std;

// Function to calculate statistical uncertainty
double Error(double average, double squared, int n) {
    if (n == 0) {
        return 0.0; // If there are no samples, return zero uncertainty
    } else {
        // Compute the standard error from the mean
        return sqrt(fabs(squared - pow(average, 2)) / n);
    }
}

int main(int argc, char *argv[]) {
    // Initialize the random number generator
    Random rnd;
    rnd.RandomImplementation();

    // Variables for averaging
    double ave_exp = 0.0, sum_exp = 0.0;  // Exponential distribution
    double ave_std = 0.0, sum_std = 0.0;  // Uniform distribution
    double ave_lrt = 0.0, sum_lrt = 0.0;  // Cauchy-Lorentz distribution

    // Parameters for the distributions
    double lambda = 1.0; // Parameter for exponential distribution
    double mean = 0.0, gamma = 1.0; // Parameters for Cauchy-Lorentz distribution

    // Number of samples to generate
    int realizations = pow(10, 5); // 100,000 samples
    vector<int> N = {1, 2, 10, 100}; // Different sample sizes for averaging

    // Create output file streams for each distribution
    ofstream outFile_exponential("1.2_exponential.dat");
    ofstream outFile_standard("1.2_standard.dat");
    ofstream outFile_lorentzian("1.2_lorentzian.dat");

    // Write headers to the output files
    outFile_exponential << "# iteration\t\t\tN = 1 :\t\t N = 2 :\t\t N = 10 :\t\t N = 100 :" << endl;
    outFile_standard << "# iteration\t\t\tN = 1 :\t\t N = 2 :\t\t N = 10 :\t\t N = 100 :" << endl;
    outFile_lorentzian << "# iteration\t\t\tN = 1 :\t\t N = 2 :\t\t N = 10 :\t\t N = 100 :" << endl;

    // Loop over the number of realizations
    for (int i = 0; i < realizations; i++) {
        sum_exp = 0.0; // Reset the sum for each distribution
        sum_std = 0.0;
        sum_lrt = 0.0;
        int sum_dim = 0; // Reset the total dimension of samples
        
        // Write the current iteration number
        outFile_exponential << setw(12) << i << setw(16);
        outFile_standard << setw(12) << i << setw(16);
        outFile_lorentzian << setw(12) << i << setw(16);

        // Loop over the different sample sizes
        for (int j = 0; j < (int)N.size(); ++j) {
            int D = N[j] - sum_dim; // Compute the number of new samples needed

            // Generate random numbers and update the sums
            for (int k = 0; k < D; ++k) { 
                sum_exp += rnd.Exponent(lambda); // Add exponential random numbers
                sum_lrt += rnd.CauchyLorentz(mean, gamma); // Add Cauchy-Lorentz random numbers
                sum_std += rnd.Rannyu(); // Add uniform random numbers

                sum_dim++; // Update the total dimension
            }

            // Compute the average for each distribution
            ave_exp = sum_exp / N[j];
            ave_std = sum_std / N[j];
            ave_lrt = sum_lrt / N[j];

            // Write the averages to the output files
            outFile_exponential << ave_exp;
            outFile_standard << ave_std;
            outFile_lorentzian << ave_lrt;

            // Add spacing between the averages
            if (j != 3) {
                outFile_exponential << setw(16);
                outFile_standard << setw(16);
                outFile_lorentzian << setw(16);
            }
        }

        // Move to the next line for the next iteration
        outFile_exponential << endl;
        outFile_standard << endl;
        outFile_lorentzian << endl;
    }

    // Close the output files
    outFile_exponential.close();
    outFile_standard.close();
    outFile_lorentzian.close();

    // Save the random number generator seed
    rnd.SaveSeed();

    return 0;
}
