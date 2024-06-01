#pragma once

#include <armadillo>
#include <iostream>

#include "random.h"


using namespace std;
using namespace arma;

class City {

	private:
		const int _dimension = 2;															// Dimensionality of the system
		vec _location;

	public:
		City () ;
		void set_location (vec location, int dimension) ;
		vec return_location ();
};