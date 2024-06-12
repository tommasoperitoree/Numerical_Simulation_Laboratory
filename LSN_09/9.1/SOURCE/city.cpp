#include <armadillo>
#include <iostream>

#include "city.h"


/**
 * @brief Default constructor for the City class.
 * 
 * This function initializes the location vector of the city.
 * 
 */
City :: City (){
   _location.resize(_dimension);
   return ;
}

/**
 * @brief Sets the location of the city.
 * 
 * This function sets the location of the city using the provided vector and dimension.
 * sIf the dimension is not compatible with the city's dimension, an error message is printed.
 * 
 * @param location The vector representing the location of the city.
 * @param dimension The dimension of the location vector.
 */
void City :: set_location (vec location, int dimension) {
	if (dimension == _dimension)
		_location = location ;
	else cerr << "ERROR : uncompatible location dimension" << endl ;	
	return ;
}

/**
 * @brief Returns the location of the city.
 * 
 * This function returns the location of the city as a vector.
 * 
 * @return The location of the city as a vector.
 */
vec City :: return_location () {
	return _location ;
}
