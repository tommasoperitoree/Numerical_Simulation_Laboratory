#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "random.h" // Custom random number generator

using namespace std;

// Function to estimate statistical uncertainty
double Error(double average, double squared, int n) {
    if (n == 1) {
        return 0.0; // If there's only one sample, no uncertainty
    } else {
        // Calculate standard error
        return sqrt(fabs(squared - pow(average, 2)) / (n - 1));
    }
}

int main(int argc, char *argv[]) {
    // Total number of throws
    int iterations = pow(10, 6); 
    // Number of blocks for averaging
    int blocks = 100; 
    // Number of throws in each block
    int blockSize = iterations / blocks; 

    // Variables for progressive averages and errors
    double sum, average, averageSquared;
    double progressiveSum = 0, progressiveSquared = 0, progressiveError;

    // Initialize the random number generator
    Random rnd;
    rnd.RandomImplementation();

    /*
    --------------------------------------------- 1.1.1 ----------------------------------------------
    
    Estimating the average value of a random variable, generated uniformly in the range [0,1)
    
    */

    // Open a file to record results for the first subproblem
    ofstream outFile("1.1.1.dat");
    outFile << "#   BLOCKS :\t\tAVE :\t\t ERROR :" << endl; // Header for the output file

    // Loop through each block for averaging
    for (int i = 1; i <= blocks; i++) {
        sum = 0;
        // Loop through each throw in the block
        for (int j = 0; j < blockSize; j++) {
            sum += rnd.Rannyu(); // Generate a random number uniformly in [0, 1)
        }
        average = sum / blockSize; // Calculate the block average
        averageSquared = pow(average, 2); // Square of the block average
        progressiveSum += average; // Update the cumulative sum
        progressiveSquared += averageSquared; // Update the cumulative squared sum
        // Calculate the progressive error
        progressiveError = Error(progressiveSum / i, progressiveSquared / i, i);
        // Write the block results to the output file
        outFile << setw(12) << i
                << setw(16) << progressiveSum / i 
                << setw(16) << progressiveError << endl;
    }
    outFile.close(); // Close the output file

    /*
    --------------------------------------------- 1.1.2 ----------------------------------------------
    
    Estimating the variance of a uniform distribution and its uncertainty
    
    */

    // Open a file to record results for the second subproblem
    outFile.open("1.1.2.dat");
    outFile << "#   BLOCKS :\t\tAVE :\t\t ERROR :" << endl; // Header for the output file

    // Reset the cumulative sums and error
    progressiveSum = 0;
    progressiveSquared = 0;
    progressiveError = 0;

    // Loop through each block
    for (int i = 1; i <= blocks; i++) {
        sum = 0;
        // Loop through each throw in the block
        for (int j = 0; j < blockSize; j++) {
            // Calculate squared deviation from the expected mean
            sum += pow(rnd.Rannyu() - 0.5, 2);
        }
        average = sum / blockSize; // Calculate the block average of squared deviations
        averageSquared = pow(average, 2); // Square of the block average
        progressiveSum += average; // Update the cumulative sum
        progressiveSquared += averageSquared; // Update the cumulative squared sum
        // Calculate the progressive error
        progressiveError = Error(progressiveSum / i, progressiveSquared / i, i);
        // Write the block results to the output file
        outFile << setw(12) << i
                << setw(16) << progressiveSum / i 
                << setw(16) << progressiveError << endl;
    }
    outFile.close(); // Close the output file

    /*
    --------------------------------------------- 1.1.3 ----------------------------------------------
    
    Estimating the chi-squared test to evaluate the uniformity of the random number generator
    
    */

    // Open a file to record results for the third subproblem
    outFile.open("1.1.3.dat");
    outFile << "#   BLOCKS :\t\tCHI2 :" << endl; // Header for the output file

    double sum_chi = 0;
    int intervals = 100; // Number of subintervals for the chi-squared test
    int tests = 100; // Number of tests
    double expected = iterations / (double)intervals; // Expected count in each subinterval
    std::vector<int> n_i(intervals); // Vector to count occurrences in each subinterval

    // Loop through each test
    for (int i = 0; i < tests; i++) {
        sum_chi = 0;
        std::fill(n_i.begin(), n_i.end(), 0); // Reset the counts

        // Generate random numbers and count occurrences in each subinterval
        for (int j = 0; j < iterations; j++) {
            double rand_chi = rnd.Rannyu();
            int index = floor(intervals * rand_chi); // Determine the subinterval
            n_i[index]++; // Increment the count for the subinterval
        }

        // Calculate the chi-squared statistic
        for (int j = 0; j < intervals; j++) {
            sum_chi += pow(n_i[j] - expected, 2) / expected; // Chi-squared formula
        }

        // Write the chi-squared result to the output file
        outFile << setw(12) << i
                << setw(16) << sum_chi 
                << endl;
    }
    outFile.close(); // Close the output file

    // Additional tests to demonstrate chi-squared distribution
    outFile.open("1.1.3_extra.dat");
    outFile << "#   BLOCKS :\t\tCHI2 :" << endl; // Header for the output file

    tests = 3000; // Increased number of tests
    iterations = 1000; // Reduced iterations for each test
    expected = iterations / (double)intervals;

    // Loop through each additional test
    for (int i = 0; i < tests; i++) {
        sum_chi = 0;
        std::fill(n_i.begin(), n_i.end(), 0); // Reset the counts

        // Generate random numbers and count occurrences in each subinterval
        for (int j = 0; j < iterations; j++) {
            double rand_chi = rnd.Rannyu();
            int index = floor(intervals * rand_chi); // Determine the subinterval
            n_i[index]++; // Increment the count
        }

        // Calculate the chi-squared statistic
        for (int j = 0; j < intervals; j++) {
            sum_chi += pow(n_i[j] - expected, 2) / expected; // Chi-squared formula
        }

        // Write the chi-squared result to the output file
        outFile << setw(12) << i 
                << setw(16) << sum_chi 
                << endl;
    }
    outFile.close(); // Close the output file

    // Save the seed for the random number generator
    rnd.SaveSeed();

    return 0;
}
