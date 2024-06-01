#pragma once

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <string>
#include "random.h"


double evalPotential(double x, double a = -2.5, double b = 1.);

double evalWaveFunctionSquared(double x, double mu, double sigma);

double evalWaveFunctionSecondDerivative(double x, double mu, double sigma);

double evalHamiltonian(double x, double mu, double sigma) ;

double evalWaveFunction(double x, double mu, double sigma);

double evalGauss (double x, double mu, double sigma) ;


void Equilibrate(int nblocks, int L, double &initialPosition, Random &rnd, double delta, double mu, double sigma);

void Metropolis(double &initialPosition, Random &rnd, double delta, int &accepted, int &attempted, double mu, double sigma);

//double Error(double AV, double AV2, int n);
double Error(double average, double averageSquared, int blockCount);

