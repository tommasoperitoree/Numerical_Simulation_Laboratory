#ifndef __Lib__
#define __Lib__

#include <iostream>
#include <cmath>
#include <vector>


#include "random.h"

using namespace std;

double Error (double average, double squared, int n) ;      // Function for statistical uncertainty estimation

double printVector (vector<double> vec) ;

double distanceSquared (vector<double> vec1, vector<double> vec2) ;

double modulusSquared (vector<double> vec) ;

void latticeStep (vector<double> &vec, Random &rnd, double step_length) ;

double latticeWalk (vector<double> &vec, Random &rnd, double step_length, int n) ;

void continuumStep (vector<double> &vec, Random &rnd, double step_length) ;

double continuumWalk (vector<double> &vec, Random &rnd, double step_length, int n) ;

void savePosition (vector<double> vec, string outpuFile, int step) ;


#endif //__Lib__