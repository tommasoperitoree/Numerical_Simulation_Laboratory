/**
 * @file lib.cpp
 * @brief This file contains the implementation of various functions used in the simulation.
 */

#include <iostream>
#include <cmath>

#include "random.h"
#include "lib.h"

using namespace std;

/**
 * @brief Evaluates the potential energy of a particle at a given position.
 * @param x The position of the particle.
 * @param a The coefficient of the quadratic term.
 * @param b The coefficient of the quartic term.
 * @return The potential energy at the given position.
 */
double evalPotential(double x, double a, double b){
	
	return a * pow(x, 2) + b * pow(x, 4);
}

/**
 * @brief Evaluates the Gaussian function at a given position.
 * @param x The position at which to evaluate the Gaussian function.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The value of the Gaussian function at the given position.
 */
double evalGauss (double x, double mu, double sigma){	

	double sigma2 = pow(sigma,2) ;
	return exp(-(pow(x-mu,2)/(2*sigma2))) ;
}

/**
 * @brief Evaluates the squared wave function at a given position.
 * @param x The position at which to evaluate the wave function.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The squared value of the wave function at the given position.
 */
double evalWaveFunctionSquared(double x, double mu, double sigma){

	return pow( abs(evalGauss(x,mu,sigma) + evalGauss(x,-1*mu,sigma)), 2);
}

/**
 * @brief Evaluates the Hamiltonian of a particle at a given position.
 * @param x The position of the particle.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The Hamiltonian at the given position.
 */
double evalHamiltonian(double x, double mu, double sigma) {

	return (( -0.5 * evalWaveFunctionSecondDerivative(x, mu, sigma) ) / evalWaveFunction(x, mu, sigma))  + evalPotential(x);
}

/**
 * @brief Equilibrates the system using the Metropolis algorithm.
 * @param blocks The number of blocks to perform.
 * @param blockSize The size of each block.
 * @param position The current position of the particle.
 * @param rnd The random number generator.
 * @param delta The step size for the Metropolis algorithm.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 */
void Equilibrate(int blocks, int blockSize, double &position, Random &rnd, double delta, double mu, double sigma){

	int accepted = 0.;
	int attempted = 0.;
	for (int j = 0; j < blocks; j++){
		for (int i = 0; i < blockSize; i++){
			Metropolis(position, rnd, delta, accepted, attempted, mu, sigma);
		}
	}
}

/**
 * @brief Performs a Metropolis step to update the position of the particle.
 * @param position The current position of the particle.
 * @param rnd The random number generator.
 * @param delta The step size for the Metropolis algorithm.
 * @param accepted The number of accepted moves.
 * @param attempted The number of attempted moves.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 */
void Metropolis(double &position, Random &rnd, double delta, int &accepted, int &attempted, double mu, double sigma){

	double future_position = position + rnd.Rannyu(-1, 1) * delta;
	double alpha = min(1., (evalWaveFunctionSquared(future_position, mu, sigma) / evalWaveFunctionSquared(position, mu, sigma)));

	double p = rnd.Rannyu();
	if (p < alpha){
		position = future_position;
		accepted++;
	}
	attempted++;
}

/**
 * @brief Evaluates the second derivative of the wave function at a given position.
 * @param x The position at which to evaluate the second derivative.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The second derivative of the wave function at the given position.
 */
double evalWaveFunctionSecondDerivative(double x, double mu, double sigma){

	double minusExp = evalGauss(x,mu,sigma) ;
	double plusExp = evalGauss(x,-1*mu,sigma) ;
	double sigma2 = pow(sigma,2) ;

	return ( pow((x - mu)/sigma2, 2) - 1./sigma2 ) * minusExp + ( pow((x + mu)/sigma2, 2) -1./sigma2 )* plusExp;
}

/**
 * @brief Evaluates the wave function at a given position.
 * @param x The position at which to evaluate the wave function.
 * @param mu The mean of the Gaussian distribution.
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The value of the wave function at the given position.
 */
double evalWaveFunction(double x, double mu, double sigma){

	double sigma2 = pow(sigma, 2) ;
	return exp(- pow(x - mu, 2) / (2 *  sigma2)) + exp(- pow(x + mu, 2) / (2 *  sigma2));
}

/**
 * @brief Calculates the statistical uncertainty of a quantity.
 * @param average The average value of the quantity.
 * @param squared The squared average value of the quantity.
 * @param blockCount The number of blocks used to calculate the average.
 * @return The statistical uncertainty of the quantity.
 */
double Error (double average, double squared, int blockCount) {	 // Function for statistical uncertainty estimation

   if (blockCount==0) {
	  return 0.;
   }
   else 
	  return sqrt (fabs (squared - pow(average, 2) ) / double(blockCount));
}