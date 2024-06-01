#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "random.h"

using namespace std;

double Error (double average, double squared, int n) {     // Function for statistical uncertainty estimation

   if (n==0) {
      return 0.;
   }
   else 
      return sqrt (fabs (squared - pow(average, 2) ) / n);
}

void printVector (vector<double> vec){

   for (int i = 0; i < vec.size(); i++){
      cout << vec[i] << " " ;
   }
   cout << endl;
}

double distanceSquared (vector<double> vec1, vector<double> vec2){

   double sum=0.;
   if(vec1.size()!=vec2.size()) 
      cerr << "incompatible vectors sizes" << endl ;

   for(int i=0; i<vec1.size(); i++)
      sum += pow(vec1[i]-vec2[i],2) ;
   
   return sum;
}

double modulusSquared (vector<double> vec){

   vector<double> orig = {0.,0.,0.} ;
   return distanceSquared(vec,orig) ;

}

void latticeStep (vector<double> &vec, Random &rnd, double step_length){
   
   /*
   vector of a point in input and function picks random direction among x,y,z to progress/regress 
   */
   
   int deg = floor(rnd.Rannyu(0,6)) ;
   int dir = floor(deg/2.);
   vec[dir] += step_length*pow(-1.,deg); 
}

double latticeWalk (vector<double> &vec, Random &rnd, double step_length, int n){

   /*
   Random Walk on a lattice with N steps and output the modulus squared of the arrival point
   */

   for (int i = 0; i < n; i++){                     // steps from i=0 to N fixed in [0,100]
      latticeStep(vec, rnd, step_length);
   }
   return modulusSquared (vec);
}

void continuumStep (vector<double> &vec, Random &rnd, double step_length){
   
   /*
   vector of a point in input and function picks random direction uniformly in the solid angle to progress/regress 
   */
  
   double theta = rnd.Theta()/2. ;
   double phi = rnd.Theta() ;

   vec[0] += step_length * sin(theta) * cos(phi) ; 
   vec[1] += step_length * sin(theta) * sin(phi) ; 
   vec[2] += step_length * cos(theta) ; 

}

double continuumWalk (vector<double> &vec, Random &rnd, double step_length, int n){

   /*
   Random Walk on the continuum with N steps and output the modulus squared of the arrival point
   */

   for (int i = 0; i < n; i++){                     // steps from i=0 to N fixed in [0,100]
      continuumStep(vec, rnd, step_length);
   }
   return modulusSquared (vec);
}

void savePosition (vector<double> vec, string outpuFile, int step){
	
	ofstream out (outpuFile, ios::app);

	out << setw(8) << step ;

	for(int i=0; i<vec.size(); i++){
		out << setw(12) << vec[i] ;
	}
	out << endl;

}
