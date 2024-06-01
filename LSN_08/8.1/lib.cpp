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
/*
double evalGauss (double x, double sigma, double mu){
	double sigma2 = pow(sigma,2) ;
	return exp(-(pow(x-mu,2)/(2*sigma2))) ;
}

double evalWaveFunction (double x, double sigma, double mu){
	return evalGauss(x,sigma,mu) + evalGauss(x,sigma,-mu) ;
}

double evalWaveFunctionSquared (double x, double sigma, double mu){
   // return pow(abs(evalWaveFunction(x,sigma,mu)), 2);
   return  pow( abs(exp(- pow(x - mu, 2) / (2 *  pow(sigma, 2))) + exp(- pow(x + mu, 2) / (2 *  pow(sigma, 2)))), 2);
}

double evalWaveFunctionSecondDerivative(double x, double sigma, double mu){

	double minusExp = exp(-0.5 * ( pow(x - mu, 2) / ( pow(sigma, 2))));
	double plusExp = exp(-0.5 * ( pow(x + mu, 2) / ( pow(sigma, 2))));

	return ((-1 /  pow(sigma, 2)) * minusExp) + ((-1 /  pow(sigma, 2)) * plusExp) + (( pow(x - mu, 2) /  pow(sigma, 4)) * minusExp) + (( pow(x + mu, 2) /  pow(sigma, 4)) * plusExp);
}

double evalPotential (double x, double a, double b){
	return a*pow(x,4) + b*pow(x,2) ;
}

double evalKinetic (double x, double sigma, double mu) {
	double sigma2 = pow(sigma,2) ;
	// double expt = exp(-2*mu*x/sigma2);
	// return ( 1/sigma2 - pow((x-mu)/sigma2,2) + (1/sigma2 - pow((x+mu)/sigma2,2))*expt ) / (2+2*expt) ;
	double kinetic = evalGauss(x,sigma,mu)*(pow((x-mu/sigma2),2)-1./sigma2) +  evalGauss(x,sigma,-mu)*(pow((x+mu/sigma2),2)-1./sigma2) ;
	return - 0.5 * kinetic / evalWaveFunction(x,sigma,mu) ;
}

double evalHamiltonian (double x, double sigma, double mu, double a, double b){
	return evalKinetic(x,sigma,mu) + evalPotential(x,a,b) ;
}


void Metropolis(double &currentPosition, Random &rnd, double delta, int &accepted, int &attempted, double sigma, double mu){

	double futurePosition = currentPosition + rnd.Rannyu(-1, 1) * delta;
	double acceptance = min(1., (evalWaveFunctionSquared(futurePosition, mu, sigma) / evalWaveFunctionSquared(currentPosition, mu, sigma)));

	double p = rnd.Rannyu();
	if (p < acceptance){
		currentPosition = futurePosition;
		accepted++;
	}
	attempted++;

}

void Equilibrate(int nblocks, int  blockSize, double &position, Random &rnd, double delta, double mu, double sigma){
	int accepted = 0;
	int attempted = 0;
	for (int j = 0; j < nblocks; j++){
		for (int i = 0; i <  blockSize; i++){
			Metropolis(position, rnd, delta, accepted, attempted, mu, sigma);
		}
	}
}
*/