#pragma once

#include <armadillo>
#include <iostream>

#include "random.h"
#include "city.h"
#include "random.h"


using namespace std;
using namespace arma;

/**
 * @file tsp.h
 * @brief This file contains the declaration of the TSP class, which represents the Traveling Salesman Problem.
 */

/**
 * @class TSP
 * @brief The TSP class represents the Traveling Salesman Problem.
 */
class TSP {

	private:
		const int _dimension = 2; /* The dimension of the problem. */
		const int _count_loss = 2; /* The number of loss parameters (actual loss, weight). */

		int _n_cities; /* The number of cities. */
		int _migration_step; /* The number of evolutions before migration, i.e. exchange of best travel among cores. */
		double _norm_order; /* The order of the norm used to calculate distances. */
		double _weight_power; /* The power used when evaluating weight. Higher means that fitter indivuals are more dominant*/
		int _population_size; /* The size of the population. */
		int _tsp_type; /* The type of the TSP problem. */
		int _n_generations; /* The number of generations. */
		int _evolution_count=0; /* The number of evolution steps. */
		string _output_path; /* The path to the output file. */

		arma::Col<double> _prob_mutations; /* The probabilities of the mutations. */
		int _prob_size; /* The size of the probabilities vector. */
			
		arma::Mat<double> _loss; /* The loss matrix. The first row represents the loss function of each travel, and the second row represents the weight, calculated in weight_evaluation. */
		arma::Mat<int> _population; /* The population matrix. Each column represents an individual travel, and each row represents a step of the travel. */
		arma::Mat<int> _new_generation; /* The new generation matrix. */
		arma::field<City> _cities; /* The field of cities. */
		Random _rand; /* The random number generator. */
		int _rank; /* The rank of the core. */

	public:

		void initialize (int rank);
		void initialize_cities_position ();
		void initialize_starting_population ();
		double distance (int i_city, int i_travel);
		double loss_function (int i_travel);
		// double loss_function (arma::Col<int> travel); // maybe unused
		void loss_evaluation ();
		int boundary_condition (int i_city);
		int boundary_condition_no_zero (int i_city);
		void check_constraints (arma::Col<int> travel);
		void order_by_loss ();
		void swap_travels (int i_travel, int j_travel);
		void cities_details_print ();
		void output_best_travel (int gen_count);
		void output_best_travel_migr (int gen_count, arma::Col<int> migrator, double loss);
		arma::Col<int> get_best_travel();
		void set_best_travel(arma::Col<int> migrator);
		void fitness_evaluation ();
		int selection_operator ();
		void swap_elements (int &i, int &j);
		void pair_permutation (int i_travel);
		void block_shift (int i_travel);
		void block_permutation (int i_travel);
		void order_inversion (int i_travel);
		void operate_mutations (int i_travel) ;
		void cross_over (int i_travel, int j_travel);
		void evolution ();
		void finalize ();
		int get_n_generations () { return _n_generations;};
		int get_migration_step () { return _migration_step;};
		int get_n_cities () { return _n_cities;};
		double get_best_loss () { return _loss(0,0);};
		int get_dimension () {return _dimension; };
		void set_cities_position (arma::Mat<double> cities_position);
		arma::Mat<double> get_cities_position ();
		
};