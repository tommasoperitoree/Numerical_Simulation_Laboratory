#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "random.h"
#include "city.h"
#include "tsp.h"
#include "armadillo"

#include "mpi.h" // parallelization


using namespace std ;
using namespace arma ;

int main(int argc, char* argv[]) {

	int size, rank;

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	// int mpi_error_code;

	TSP tsp ;
	tsp.initialize(rank) ; // initialize the seed of the rnd as the rank
	Random rnd; 
	rnd.RandomImplementation();

	// variables needed all initialized before the loops
	arma::Mat<double> cities_position ;
	int n_cities = tsp.get_n_cities(), dimension = tsp.get_dimension();
	cities_position.resize(n_cities, dimension);	
	
	arma::Mat<int> migrators_mat;
	arma::Col<int> migrator ;
	migrator.resize(n_cities) ;
	migrators_mat.resize(n_cities,size);
	
	double i_loss ;
	arma::Col<double> loss_vec;
	// if ( rank == 0 ) 
	loss_vec.resize(size);

	// parameters
	int generations = tsp.get_n_generations() ;
	int migration_step = tsp.get_migration_step() ;
	int min_loss_rank ;

	// ------------ normalize cities positions for all ranks ------------ // 
	if (rank == 0) { // rank 0 has the initialisation of the cities
      cities_position = tsp.get_cities_position();
      // cities_position.print("Cities Position:");
   }
	// now they all get the cities position
	MPI_Bcast(cities_position.memptr(), cities_position.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);	

	if (rank != 0) { // other ranks set the cities position and can now initialize the starting population
		tsp.set_cities_position(cities_position) ;
		tsp.initialize_starting_population() ;
		tsp.order_by_loss() ;
	}

	// --------------------------- evolution --------------------------- //

	// cout << "Rank " << rank << " starting evolution" << endl;

	for (int gen=0; gen<generations; gen++) {
		tsp.evolution() ;
		// cout << "Rank " << rank << " finished evolution on gen " << gen << endl;
		// MPI_Barrier(MPI_COMM_WORLD);

		if ( (gen % migration_step) == 0 and gen != 0){
			// Synchronize all processes before checking for migration step
			if (rank == 0) cout << endl << "Starting migration on gen " << gen << endl << endl;
			migrator = tsp.get_best_travel() ;
			i_loss = tsp.get_best_loss() ;
			// all nodes send to rank 0 their best travel
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(migrator.memptr(), migrator.size(), MPI_INTEGER, migrators_mat.memptr(), migrator.size(), MPI_INTEGER, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(&i_loss, 1, MPI_DOUBLE, loss_vec.memptr(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
			// steve : MPI_Gather(&i_loss, 1, MPI_DOUBLE, rank == 0 ? loss_vec.memptr() : NULL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;


			if (rank == 0) { // have rank 0 save best overall travel and shuffle the matrix
				// cout << "matrix and loss vec for gen " << gen << endl; 
				//migrators_mat.print(); cout << endl;
				//loss_vec.print();
				min_loss_rank = arma::index_min(loss_vec);
				tsp.output_best_travel_migr(gen,migrators_mat.col(min_loss_rank),loss_vec(min_loss_rank));
				// finally shuffle
				migrators_mat = arma::shuffle(migrators_mat, 1) ;
			}

			// now they all get the reshuffled matrix
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(migrators_mat.memptr(), migrators_mat.size(), MPI_INTEGER, 0, MPI_COMM_WORLD);

			// each rank now has the whole matrix, the coloumn corresponding to its rank is its new given best travel
			migrator = migrators_mat.col(rank) ;
			// migrator.print("migrator for rank " + to_string(rank));
			tsp.set_best_travel(migrator) ; 
			// migrator.clear(); 
			// migrators_mat.clear();
			// cout << "Rank " << rank << " saved new best travel" << endl ;
			tsp.loss_evaluation(); // need to recalculate all losses of the population with the new element
			// cout << "Rank " << rank << " finished migration on gen " << gen << endl;
		}
		// migrator.clear(); 
		// migrators_mat.clear() ;

	}

	tsp.finalize() ;
	MPI_Finalize();

   return 0;
}