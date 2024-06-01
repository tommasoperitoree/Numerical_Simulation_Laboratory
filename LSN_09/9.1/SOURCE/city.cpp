#include <armadillo>
#include <iostream>

#include "city.h"

City :: City (){
   _location.resize(_dimension);
   return ;
}

void City :: set_location (vec location, int dimension) {
	if (dimension == _dimension)
		_location = location ;
	else cerr << "ERROR : uncompatible location dimension" << endl ;	
	return ;
}

vec City :: return_location () {
	return _location ;
}
