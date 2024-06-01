#include <iostream>
#include <cmath>

#include "random.h"
#include "lib.h"

using namespace std;


double evalPotential(double x, double a, double b){
	
	return a * pow(x, 2) + b * pow(x, 4);
}

double evalGauss (double x, double mu, double sigma){

	double sigma2 = pow(sigma,2) ;
	return exp(-(pow(x-mu,2)/(2*sigma2))) ;
}

double evalWaveFunctionSquared(double x, double mu, double sigma){

	return pow( abs(evalGauss(x,mu,sigma) + evalGauss(x,-1*mu,sigma)), 2);
}

double evalHamiltonian(double x, double mu, double sigma) {

	return (( -0.5 * evalWaveFunctionSecondDerivative(x, mu, sigma) ) / evalWaveFunction(x, mu, sigma))  + evalPotential(x);
}

void Equilibrate(int blocks, int blockSize, double &position, Random &rnd, double delta, double mu, double sigma){

	int accepted = 0.;
	int attempted = 0.;
	for (int j = 0; j < blocks; j++){
		for (int i = 0; i < blockSize; i++){
			Metropolis(position, rnd, delta, accepted, attempted, mu, sigma);
		}
	}
}

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

double evalWaveFunctionSecondDerivative(double x, double mu, double sigma){

	double minusExp = evalGauss(x,mu,sigma) ;
	double plusExp = evalGauss(x,-1*mu,sigma) ;
	double sigma2 = pow(sigma,2) ;

	return ( pow((x - mu)/sigma2, 2) - 1./sigma2 ) * minusExp + ( pow((x + mu)/sigma2, 2) -1./sigma2 )* plusExp;
}

double evalWaveFunction(double x, double mu, double sigma){

	double sigma2 = pow(sigma, 2) ;
	return exp(- pow(x - mu, 2) / (2 *  sigma2)) + exp(- pow(x + mu, 2) / (2 *  sigma2));
}

double Error (double average, double squared, int blockCount) {	 // Function for statistical uncertainty estimation

   if (blockCount==0) {
	  return 0.;
   }
   else 
	  return sqrt (fabs (squared - pow(average, 2) ) / double(blockCount));
}