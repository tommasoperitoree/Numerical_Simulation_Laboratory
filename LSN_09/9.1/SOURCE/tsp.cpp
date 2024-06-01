#include <armadillo>
#include <iostream>
#include <set>

#include "tsp.h"
#include "random.h"
#include "mpi.h"

using namespace std;

/**
 * This file contains the implementation of the TSP (Traveling Salesman Problem) class.
 * It includes the necessary headers, defines the class methods, and initializes the TSP parameters.
 * The TSP class handles the initialization of cities, population, and loss evaluation.
 * It also includes mutation and crossover methods for solving the TSP.
 */

/**
 * Initializes the TSP (Traveling Salesman Problem) simulation.
 * 
 * This function reads input from a file and initializes the necessary variables and data structures
 * for the simulation. It sets the TSP type, norm order, number of cities, population size, and other
 * parameters based on the input file. It also initializes the random number generator and opens an
 * output file for writing simulation results.
 * 
 * After initializing the necessary variables, it calls other methods to initialize the cities' positions,
 * create the starting population, and order the population by loss.
 */
void TSP :: initialize() {
	// Initialize the random number generator
	_rand.RandomImplementation();
	std::ifstream input("../INPUT/input.dat");
	std::ofstream coutf;
	std::string property;
	while (!input.eof()) {
		input >> property;
		if (property == "TSP_TYPE") {
			input >> _tsp_type;
			if (_tsp_type > 1) {
				std::cerr << "PROBLEM: unknown simulation type" << std::endl;
				exit(EXIT_FAILURE);
			}
			if (_tsp_type == 0) {
				_output_path = "../OUTPUT/Circle/";
				coutf.open(_output_path + "output.dat");
				coutf << "GENETIC ALGORITHM WITH CITIES AROUND CIRCUMFERENCE" << std::endl;
			} else if (_tsp_type == 1) {
				_output_path = "../OUTPUT/Square/";
				coutf.open(_output_path + "output.dat");
				coutf << "GENETIC ALGORITHM WITH CITIES INSIDE A SQUARE" << std::endl;
			}
		} else if (property == "NORM_ORDER") {
			input >> _norm_order;
			coutf << "Norm order = " << _norm_order << std::endl;
		} else if (property == "N_CITIES") {
			input >> _n_cities;
			_cities.set_size(_n_cities);
			coutf << "Number of cities = " << _n_cities << std::endl;
		} else if (property == "N_INDIVIDUALS") {
			input >> _population_size;
			coutf << "Population size = " << _population_size << std::endl;
		} else if (property == "WEIGHT_POWER") {
			input >> _weight_power;
			coutf << "Weight power = " << _weight_power << std::endl;
		} else if (property == "N_GENERATIONS") {
			input >> _n_generations;
			coutf << "Number of generations = " << _n_generations << std::endl;
		} else if (property == "PROB_MUTATIONS") {
			input >> _prob_size;
			double prob;
			// cout << "prob size " << _prob_size << endl;
			for (int i = 0; i < _prob_size; i++) {
				input >> prob;
				//cout << "prob " << i << " = " << prob << endl;
				_prob_mutations = arma::join_vert(_prob_mutations, arma::Col<double>({prob}));
			}
			coutf << "Probabilities of mutations = " << _prob_mutations.t() << std::endl;
		} else if (property == "ENDINPUT") {
			coutf << "Reading input completed!" << std::endl;
			break;			
		} else {
			std::cerr << "PROBLEM: unknown input" << std::endl;
		}
	}
	input.close();
	// Ensure that _n_cities and _population_size have been initialized before using them
	if (_n_cities > 0 && _population_size > 0) {
		_population.resize(_n_cities, _population_size);
		_loss.resize(_count_loss, _population_size);
	} else {
		std::cerr << "ERROR: Invalid _n_cities or _population_size" << std::endl;
		return;
	}
	coutf << "System initialized!" << std::endl;
	coutf.close();
	// Ensure the following methods are defined and working correctly
	this->initialize_cities_position();
	this->initialize_starting_population();
	this->order_by_loss();
}

/**
 * Calculates the distance between two cities.
 * @param i_city The index of the first city.
 * @param i_travel The index of the travel (individual) in the population.
 * @return The distance between the two cities.
 */
double TSP :: distance(int i_city, int i_travel) {
	vec distance = _cities(_population.at(i_city, i_travel)).return_location() - 
						_cities(_population.at(boundary_condition(i_city + 1), i_travel)).return_location();
	return arma::norm(distance, _norm_order);
}

/**
 * Applies periodic boundary conditions to a city index.
 * @param i_city The city index.
 * @return The city index after applying periodic boundary conditions.
 */
int TSP :: boundary_condition(int i_city) {
	if (i_city >= _n_cities) i_city -= _n_cities;
	else if (i_city < 0) i_city += _n_cities;
	return i_city;
}

/**
 * Applies periodic boundary conditions to a city index, excluding zero.
 * @param i_city The city index.
 * @return The city index after applying periodic boundary conditions.
 */
int TSP :: boundary_condition_no_zero(int i_city) {
	if (i_city >= (_n_cities - 1)) i_city -= (_n_cities - 1);
	else if (i_city < 0) i_city += (_n_cities - 1);
	return i_city;
}


/**
 * Initializes the positions of the cities.
 */
void TSP :: initialize_cities_position() {
	vec position;
	position.resize(_dimension);

	for (int i_city = 0; i_city < _n_cities; i_city++) {
		if (_tsp_type == 0) { // cities position on a circumference of unit radius
			double theta = _rand.Theta();
			position(0) = cos(theta);
			position(1) = sin(theta);
			_cities(i_city).set_location(position, position.size());
		} else { // _tsp_type == 1 i.e. cities position inside a square of side 2
			position(0) = _rand.Rannyu(-1., 1.);
			position(1) = _rand.Rannyu(-1., 1.);
			_cities(i_city).set_location(position, position.size());
		}
	}
	ofstream coutf;
	coutf.open(_output_path + "output.dat", ios::app);
	coutf << "Cities position initialized!" << endl;
	cities_details_print();

	return;
}

/**
 * Initializes the starting population.
 */
void TSP :: initialize_starting_population() {
	arma::Col<int> cities_index = arma::regspace<arma::Col<int>>(1, _n_cities - 1);

	for (int i = 0; i < _population_size; i++) {
		arma::Col<int> missing_index = cities_index;
		arma::Col<int> temporary_path = {0};

		for (int j = _n_cities; j > 1; j--) {
			int selected_index = int(_rand.Rannyu() * (j - 1));
			int selection = missing_index(selected_index);

			temporary_path.insert_rows(temporary_path.n_elem, arma::Col<int>({selection}));
			missing_index.shed_row(selected_index);
		}
		check_constraints(temporary_path);
		_population.col(i) = temporary_path;
	}
	ofstream coutf;
	coutf.open(_output_path + "output.dat", ios::app);
	coutf << "Initial population created!" << endl;

	this->loss_evaluation();

	return;
}

/**
 * Checks if there are any duplicate cities within an individual.
 * @param travel The travel (individual) to check.
 */
void TSP :: check_constraints(arma::Col<int> travel) {
	set<int> num_set;
	for (int i_city = 0; i_city < travel.n_elem; i_city++) {
		if (!num_set.insert(travel(i_city)).second) {
			cerr << "ERROR: duplicate city within an individual" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

/**
 * Calculates the loss function of a defined travel.
 * @param i_travel The index of the travel (individual) in the population.
 * @return The loss value of the travel.
 */
double TSP :: loss_function(int i_travel) {
	double loss = 0.;
	for (int i_city = 0; i_city < _n_cities; i_city++)
		loss += distance(i_city, i_travel);

	return loss;
}

/**
 * Evaluates the loss function for each travel in the population.
 */
void TSP :: loss_evaluation() {
	for (int i_travel = 0; i_travel < _population_size; i_travel++)
		_loss.at(0, i_travel) = loss_function(i_travel);
	
	this->order_by_loss();
	this->weight_evaluation();

	return;
}

/**
 * Evaluates the weights for each travel in the population.
 */
void TSP :: weight_evaluation() {
	double loss_norm = 0; // normalization factor for the weights
	// double loss_avg = 0;
	// int evaluation_range = 5;
	for (int i_travel = 0; i_travel < _population_size; i_travel++)
		loss_norm += pow(_loss.at(0, i_travel), -1.*_weight_power); // sum of the inverse of the loss function
	for (int i_travel = 0; i_travel < _population_size; i_travel++){
		_loss.at(1, i_travel) = pow(_loss.at(0, i_travel), -1*_weight_power) / loss_norm;	// weight of the travel	
		// cout << i_travel << ": loss = " << _loss.at(0, i_travel) << " weight = " << _loss.at(1,i_travel) << endl;
		// if (i_travel < evaluation_range) loss_avg += _loss.at(0, i_travel);	
	}
	// loss_avg /= evaluation_range;
	// cout << "loss average = " << loss_avg << endl;
	return;
}


/**
 * Orders the population by loss value in ascending order.
 */
void TSP :: order_by_loss() {
	for (int i = 0; i < _population_size - 1; i++)
		for (int j = i + 1; j < _population_size; j++)
			if (_loss.at(0, i) > _loss.at(0, j))
				swap_travels(i, j);
	return;
}

/**
 * Swaps two travels in the population.
 * @param i_travel The index of the first travel.
 * @param j_travel The index of the second travel.
 */
void TSP :: swap_travels(int i_travel, int j_travel) {
	arma::Col<int> temp_travel = _population.col(i_travel);
	_population.col(i_travel) = _population.col(j_travel);
	_population.col(j_travel) = temp_travel;

	arma::Col<double> temp_loss = _loss.col(i_travel);
	_loss.col(i_travel) = _loss.col(j_travel);
	_loss.col(j_travel) = temp_loss;
}

/**
 * Prints the details of the cities to a file.
 */
void TSP :: cities_details_print() {
	ofstream coutcit(_output_path + "cities_details.dat");
	coutcit << "# CITY: \t\t POSITION X: \t POSITION Y: " << endl;
	for (int i = 0; i < _n_cities; i++) {
		vec pos = _cities(i).return_location();
		coutcit << setw(6) << i
				<< setw(18) << pos(0)
				<< setw(15) << pos(1) << endl;
	}
	return;
}

/**
 * Prints the details of a travel to a file.
 * @param i_travel The index of the travel.
 */
void TSP :: output_best_travel(int gen_count) {
	ofstream couttrv(_output_path + "/GEN_BEST/best_travel_gen" + to_string(gen_count) + ".dat");
	couttrv << "# CITY: \t\t POSITION X: \t POSITION Y: " << endl;
	
	arma::Col<int> travel = _population.col(0);
	for (int i = 0; i < _n_cities; i++) {
		int index = travel(i);
		vec pos = _cities(index).return_location();
		couttrv << setw(6) << index
				 << setw(18) << pos(0)
				 << setw(15) << pos(1) << endl;
	}
	return;
}

/**
 * Selects a travel from the population based on the weights.
 * @return The index of the selected travel.
 */
int TSP :: selection_operator() {
	double draw = _rand.Rannyu();
	double progr_weight = 0;
	for (int i_travel = 0; i_travel < _population_size; i_travel++) {
		progr_weight += _loss.at(1, i_travel); // cumulative weight of the travels
		if (draw < progr_weight) return i_travel;
	}
	return 0;
}

/**
 * Swaps two elements.
 * @param i The first element.
 * @param j The second element.
 */
void TSP :: swap_elements(int &i, int &j) {
	int temp = i;
	i = j;
	j = temp;
	return;
}

// methods to apply mutations to the population

/**
 * Performs a pair permutation on the selected row of the population matrix.
 *
 * This function selects a row from the population matrix, performs a pair permutation on it,
 * and updates the population matrix with the modified row. The pair permutation involves
 * selecting two adjacent cities in the row and swapping their positions.
 *
 * @param i_travel The index of the selected row in the population matrix.
 */
void TSP :: pair_permutation(int i_travel) {
	
	arma::Col<int> temp_travel = _population.col(i_travel);
	// temp_travel.print("selected row");
	temp_travel.shed_row(0);

	// Select a random index to start the pair permutation
	int i_start = int(_rand.Rannyu() * (_n_cities - 2)) + 1;
	// cout << "Swapping elements at index " << i_start << " and " << i_start + 1 << endl;

	// Swap the selected elements
	swap_elements(temp_travel(i_start), 
					  temp_travel(boundary_condition_no_zero(i_start + 1)));

	// Insert the modified row back into the population matrix
	arma::Col<int> zero = {0};
	temp_travel.insert_rows(0,zero) ;
	check_constraints(temp_travel);
	_population.col(i_travel) = temp_travel;

	// _population.col(i_travel).print("selected row after pair permutation");

	return;
}

/**
 * @brief Performs a block shift operation on the selected travel route.
 * 
 * This function shifts a block of cities in the selected travel route. The block is randomly determined
 * by selecting a starting city and the size of the block. The block is then shifted by a random amount.
 * The shifted travel route is stored in a new vector.
 * 
 * @param i_travel The index of the selected travel route.
 */
void TSP :: block_shift (int i_travel) {
	
	arma::Col<int> temp_travel = _population.col(i_travel);
	// temp_travel.print("selected row");
	temp_travel.shed_row(0);

	int block_start = int(_rand.Rannyu()*(_n_cities-2))+1;
		// cout << "block start " << block_start << endl;
	int block_size = int(_rand.Rannyu()*(_n_cities-1))+1;
		// cout << "block size " << block_size << endl; 
	int block_end = boundary_condition_no_zero(block_start+block_size-1);
	 	// cout << "block end " << block_end << endl;
	int block_shift = int(_rand.Rannyu()*(_n_cities-1))+1;
		// cout << "block shift " << block_shift << endl;

		for (int i=block_end; ; i--){
			swap_elements(temp_travel(i),
							  temp_travel(boundary_condition_no_zero(i+block_shift)));
			if (i==0) i=_n_cities-1;
			if (i==block_start) break;
		}

	arma::Col<int> zero = {0};
	temp_travel.insert_rows(0,zero) ;
	check_constraints(temp_travel);
	_population.col(i_travel) = temp_travel;
		
	// _population.col(i_travel).print("selected row after block shift");

}

/**
 * Performs a block permutation on the selected row of the population matrix.
 *
 * This function selects a row from the population matrix, performs a block permutation on it,
 * and updates the population matrix with the modified row. The block permutation involves
 * selecting two blocks of cities from the row and swapping their positions.
 *
 * @param i_travel The index of the selected row in the population matrix.
 */
void TSP :: block_permutation (int i_travel) {

	arma::Col<int> temp_travel = _population.col(i_travel);
		// temp_travel.print("selected row");
	temp_travel.shed_row(0);

	int blocks_size = int(_rand.Rannyu()*((_n_cities-4)/2.))+1;
		// cout << "blocks size = " << blocks_size << endl; 	
	int block_1_start = int(_rand.Rannyu()*(_n_cities-2));
		// cout << "block 1 start = " << block_1_start << endl;
	int block_2_start;
	bool no_overlap = false;

	do{
		block_2_start = int(_rand.Rannyu()*(_n_cities-2)); 
		no_overlap = (block_2_start < block_1_start and block_2_start + blocks_size < block_1_start) or 
					(block_1_start < block_2_start and block_1_start + blocks_size < block_2_start);
	} while (!no_overlap);
		// cout << "block 2 start = " << block_2_start << endl;

	for (int i=0; i<blocks_size; i++){
		swap_elements(temp_travel(boundary_condition_no_zero(block_1_start+i)),
						  temp_travel(boundary_condition_no_zero(block_2_start+i)));
	}

	arma::Col<int> zero = {0};
	temp_travel.insert_rows(0,zero) ; // temp_travel.print("new travel");

	check_constraints(temp_travel);
	_population.col(i_travel) = temp_travel;

	return ;
}

/**
 * Reverses a block of cities in the travel order.
 *
 * @param i_travel The index of the travel order to modify.
 */
void TSP :: order_inversion (int i_travel) {

	arma::Col<int> temp_travel = _population.col(i_travel);
		// temp_travel.print("selected row");
	temp_travel.shed_row(0);

	int block_start = int(_rand.Rannyu()*(_n_cities-2));
		// cout << "block start " << block_start << endl;
	int block_size = int(_rand.Rannyu()*(_n_cities-2))+1;	
		// cout << "block size " << block_size << endl;

	for(int i=0; i<block_size/2; i++){
		swap_elements(temp_travel(boundary_condition_no_zero(block_start+i)),
						  temp_travel(boundary_condition_no_zero(block_start+block_size-i)));
	}
		
	arma::Col<int> zero = {0};
	temp_travel.insert_rows(0,zero) ; // temp_travel.print("new travel");
	check_constraints(temp_travel);
	_population.col(i_travel) = temp_travel;


	return ;
}

/**
 * @brief Performs mutations on a given travel path.
 * 
 * This function applies different mutation operations to a given travel path based on the probabilities defined in _prob_mutations.
 * The mutation operations include pair permutation, block shift, block permutation, and order inversion.
 * 
 * @param i_travel The index of the travel path to operate mutations on.
 */
void TSP :: operate_mutations (int i_travel) {
	
	for (int i=1; i<_prob_mutations.size(); i++){
		double prob = _rand.Rannyu();
		if (prob < _prob_mutations(i)){
			switch (i){
				case 0:
					break;
				case 1:
					pair_permutation(i_travel);
					break;
				case 2:
					block_shift(i_travel);
					break;
				case 3:
					block_permutation(i_travel);
					break;
				case 4:
					order_inversion(i_travel);
					break;
				default:
					cerr << "ERROR: unknown mutation operation" << endl;
					break;
			}
		}
	}

	return ;

}

/**
 * @brief Performs the crossover operation for the TSP problem.
 * 
 * This function is responsible for performing the crossover operation, which is a genetic operator used in evolutionary algorithms for solving the Traveling Salesman Problem (TSP). The crossover operation combines genetic material from two parent individuals to create new offspring individuals.
 * 
 * In the context of the TSP, the crossover operation involves selecting a subset of cities from one parent and arranging them in the same order as they appear in the other parent. This process helps to create diverse and potentially better solutions for the TSP.
 * 
 * This function does not return any value.
 */
void TSP :: cross_over (int i_travel, int j_travel) {	
	
	arma::Col<int> parent_1 = _population.col(i_travel), parent_2 = _population.col(j_travel);
	parent_1.shed_row(0); parent_2.shed_row(0);
	int crossover_index = _rand.Rannyu()*(_n_cities-1)+1;
		// cout << "crossover index " << crossover_index << endl;
	arma::Col<int> child_1 = parent_1.subvec(0,crossover_index-1), child_2 = parent_2.subvec(0,crossover_index-1);
	
	set<int> child_1_set, child_2_set;
	for(int i=0; i<child_1.size(); i++){ // filling the sets with the cities already in the children
		child_1_set.insert(child_1(i));
		child_2_set.insert(child_2(i));
	}
	/* parent_1.print("parent 1");
	parent_2.print("parent 2");
	child_1.print("child 1");
	child_2.print("child 2"); */

	int count_1 = crossover_index, count_2 = crossover_index, i_city = 0;
	do{	// fill children with the remaining cities reading from order of parents
			if (child_1_set.insert(parent_2(i_city)).second){	// checks if the city is not already in the child
				child_1 = arma::join_cols(child_1, arma::Col<int>({parent_2(i_city)}));
				count_1++;
					/* cout << "inserting city " << parent_2(i_city) << " index " << i_city << " in child 1" << endl;
					child_1.print("progress child 1");
					cout << "count 1 progress " << count_1 << endl; */
			}
			if (child_2_set.insert(parent_1(i_city)).second){
				child_2 = arma::join_cols(child_2, arma::Col<int>({parent_1(i_city)}));
				count_2++;
					/* cout << "inserting city " << parent_1(i_city) << " index " << i_city << " in child 2" << endl;
					child_2.print("progress child 2");
					cout << "count 2 progress " << count_2 << endl; */
			}
		i_city++; // cout << "i_city " << i_city << endl;
	} while (count_1 < _n_cities-1 or count_2 < _n_cities-1);

	arma::Col<int> zero = {0};
	child_1.insert_rows(0,zero); child_2.insert_rows(0,zero);
	// child_1.print("child 1"); child_2.print("child 2");

	check_constraints(child_1);
	_new_generation = arma::join_horiz(_new_generation,child_1);
	// _new_generation.print("new generation");
	operate_mutations(_new_generation.n_cols-1); // operating mutations on the new child 1

	check_constraints(child_2);
	_new_generation = arma::join_horiz(_new_generation,child_2);
	// _new_generation.print("new generation");
	operate_mutations(_new_generation.n_cols-1); // operating mutations on the new child 2

	return ;
}

// to be deleted once code is done 
void TSP :: mutation_check_debug () {
	for (int i=0; i<10001; i++){
		int travel_pick = _rand.Rannyu()*_population_size; cout << endl;
		int second_pick = _rand.Rannyu()*_population_size; cout << endl;
		cout << "mutation " << i << endl; cout << endl;
		cross_over(travel_pick,second_pick);
	}
	return ;
}

void TSP :: evolution () {

	_new_generation.resize(_n_cities,0);
	// _new_generation.print("new generation");
	do {
		int parent_1 = selection_operator();
		int parent_2 = selection_operator();
		cross_over(parent_1,parent_2);
		// cout << "new generation size " << _new_generation.n_cols << endl;
	} while (_new_generation.n_cols < _population_size);

	_population = _new_generation;
	_new_generation.resize(0,0);
	_new_generation.clear();
	// _population.print("new population");
	loss_evaluation();

	std::ofstream coutf;
	coutf.open(_output_path + "output.dat", ios::app);
	_evolution_count++;
	coutf << "Evolution " << _evolution_count << " completed!" << std::endl;


	return ;
}

void TSP :: finalize () {
	_rand.SaveSeed();
	ofstream coutf;
	coutf.open(_output_path + "output.dat",ios::app);
	coutf << endl << "Simulation completed!" << endl;
	coutf.close();
}