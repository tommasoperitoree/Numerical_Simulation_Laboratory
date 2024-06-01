#include <iostream>
#include <cmath>

#include "random.h"

using namespace std;

double Error (double average, double squared, int n) {     // Function for statistical uncertainty estimation
   if (n==0) {
      return 0.;
   }
   else 
      return sqrt (fabs (squared - pow(average, 2) ) / n);
}