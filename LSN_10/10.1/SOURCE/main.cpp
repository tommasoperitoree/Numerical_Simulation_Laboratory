#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "random.h"
#include "city.h"
#include "tsp.h"
#include "armadillo"

#include "mpi.h"

using namespace std ;


int main(int argc, char* argv[]){

	int size, rank;

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;

	TSP tsp ;
	tsp.initialize() ;
	Random rnd; rnd.RandomImplementation();

	// variables needed for exchange of best travel among nodes
	int* migrator = new int (tsp.get_n_cities()); 

	int i_sender, i_receiver, i_tag;
	int n_cities = tsp.get_n_cities() ;
	int n_migrations = tsp.get_n_migrations() ;


	for (int gen=0; gen<n_migrations; gen++) {
		tsp.evolution() ;
		// tsp.output_best_travel(g);

		if(! gen % n_migrations ) {	// exchange of best travel among nodes once every n_migrations
			i_tag = gen + 1 ;
			// rank 0 chooses the nodes that will exchange the best travel			
			if (rank == 0){
				i_sender = int(rnd.Rannyu(0, size)) ;
					do{
						i_receiver = int(rnd.Rannyu(0, size)) ;
						cout << "i_sender = " << i_sender << ", i_receiver = " << i_receiver << endl ;
					} while (i_receiver == i_sender) ; // don't send to yourself
			}

			// rank 0 broadcasts the indices of the nodes that will exchange the best travel
			MPI_Bcast(&i_sender, 1, MPI_INTEGER, 0, MPI_COMM_WORLD) ;
			MPI_Bcast(&i_receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD) ;

			if (rank == i_sender){ // sender sends
				migrator = tsp.get_best_travel() ;
				MPI_Send(migrator, n_cities, MPI_INTEGER, i_receiver, i_tag, MPI_COMM_WORLD);
			}
         if(rank == i_receiver){ // receiver receives
            MPI_Recv(migrator, n_cities, MPI_INTEGER, i_sender, i_tag, MPI_COMM_WORLD, &stat);
				tsp.set_best_travel(migrator) ; // function to be written still
            }
			// need to forcefully recalculate losses as the population is modified for the receiver
		}
		// if I want to print the best travel at each generation I need to add a function to check amongst nodes which one has the best travel	

	}

	tsp.finalize() ;

	MPI_Finalize();

   return 0;
}