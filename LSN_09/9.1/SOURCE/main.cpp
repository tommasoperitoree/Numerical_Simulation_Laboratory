#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "random.h"
#include "city.h"
#include "tsp.h"

using namespace std ;


int main () {

	TSP tsp ;
	tsp.initialize() ;

	
	for (int g=0; g<tsp.get_n_generations(); g++) {
		tsp.evolution() ;
		tsp.output_best_travel(g);
	}

	tsp.finalize() ;

   return 0;
}