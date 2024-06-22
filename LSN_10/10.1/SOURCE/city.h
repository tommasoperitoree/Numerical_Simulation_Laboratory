#pragma once

#include <armadillo>
#include <iostream>

#include "random.h"

using namespace std;
using namespace arma;

/**
 * @brief Represents a city in a two-dimensional system.
 * 
 * The City class stores the location of a city in a two-dimensional system.
 * It provides methods to set and retrieve the location of the city.
 */
class City {
	private:
		const int _dimension = 2;	// Dimensionality of the system
		arma::Row<double> _location;				// Location of the city

	public:
		/**
		 * @brief Default constructor for the City class.
		 */
		City();

		/**
		 * @brief Sets the location of the city.
		 * 
		 * @param location The location of the city as a vector.
		 * @param dimension The dimensionality of the system.
		 */
		void set_location(arma::Row<double> location, int dimension);

		/**
		 * @brief Returns the location of the city.
		 * 
		 * @return The location of the city as a vector.
		 */
		arma::Row<double> return_location();
		int get_dimension() {return _dimension; };
};
