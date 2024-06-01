/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "system.h"

using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
	if(_sim_type == 0) this->Verlet();	// Perform a MD step
	else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
	_nattempts += _npart; //update number of attempts performed on the system
	return;
}

void System :: Verlet(){
	double xnew, ynew, znew;
	for(int i=0; i<_npart; i++){ //Force acting on particle i
		_fx(i) = this->Force(i,0);
		_fy(i) = this->Force(i,1);
		_fz(i) = this->Force(i,2);
	}
	for(int i=0; i<_npart; i++){ //Verlet integration scheme
		xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
		ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
		znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
		_particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
		_particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
		_particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
		_particle(i).acceptmove(); // xold = xnew
		_particle(i).setposition(0, xnew);
		_particle(i).setposition(1, ynew);
		_particle(i).setposition(2, znew);
	}
	_naccepted += _npart;
	return;
}

double System :: Force(int i, int dim){
	double f=0.0, dr;
	vec distance;
	distance.resize(_ndim);
	for (int j=0; j<_npart; j++){
		if(i != j){
			distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
			distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
			distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
			dr = sqrt(arma::dot(distance,distance) ); // dot is dot-product, hence this calculates the modulus of distance
			if(dr < _r_cut){
				f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
			}
		}
	}
	return f;
}

void System :: move(int i){ // Propose a MC move for particle i
	if(_sim_type == 3) { //Gibbs sampler for Ising
		double prob_i_up = 1. / (1. + exp (-2. * _beta * (_J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin()) + _H))); // conditional distribution of heat bath
		if (_rnd.Rannyu() < prob_i_up) _particle(i).setspin(+1);
		else _particle(i).setspin(-1);
		_naccepted++;
	} else {					 // M(RT)^2
		if(_sim_type == 1){			 // LJ system
			vec shift(_ndim);			 // Will store the proposed translation
			for(int j=0; j<_ndim; j++){
				shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
			}
			_particle(i).translate(shift, _side);	//Call the function Particle::translate
			if(this->metro(i)){ //Metropolis acceptance evaluation
				_particle(i).acceptmove();
				_naccepted++;
			} else _particle(i).moveback(); //If translation is rejected, restore the old configuration
		} else {									// Ising 1D (_sym_type = 2)
			if(this->metro(i)){		 //Metropolis acceptance evaluation for a spin flip involving spin i
				_particle(i).flip();	//If accepted, the spin i is flipped
				_naccepted++;
			}
		}
	}
	return;
}

bool System :: metro(int i){ // Metropolis algorithm
	bool decision = false;
	double delta_E, acceptance;
	if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
	else delta_E = 2.0 * _particle(i).getspin() * 
								 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
	acceptance = exp(-_beta*delta_E);
	if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
	return decision;
}

double System :: Boltzmann(int i, bool xnew){
	double energy_i=0.0;
	double dx, dy, dz, dr;
	for (int j=0; j<_npart; j++){
		if(j != i){
			dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
			dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
			dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
			if(dr < _r_cut){
				energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}
	return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
	return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
	if(i >= _npart) i = i - _npart;
	else if(i < 0)	i = i + _npart;
	return i;
} 

void System :: initialize(string input_dir){ // Initialize the System object according to the content of the input files in the ../_INPUT/ directory

	int p1, p2; // Read from ../_INPUT/Primes a pair of numbers to be used to initialize the RNG
	ifstream Primes("../_INPUT/Primes");
	Primes >> p1 >> p2 ;
	Primes.close();
	int seed[4]; // Read the seed of the RNG
	ifstream Seed("../_INPUT/seed.in");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	_rnd.SetRandom(seed,p1,p2);

	_input_dir = input_dir;

	ifstream input(_input_dir); // Start reading _input_dir
	string property;
	double delta;
	bool init = false;

	ofstream coutf;
	ofstream couta;
	
	while ( !input.eof() ){
		input >> property;
		if(property == "SIMULATION_TYPE"){ 
			input >> _sim_type;
			if(_sim_type == 0){	// _sim_type format : SIMULATION_TYPE 0 #phase #equilibration
				int phase;
				input >> _int_phase;
				input >> _equilibration;
				if (_int_phase == 1) _phase = "Gas";
				else if (_int_phase == 2) _phase = "Liquid";
				else if (_int_phase == 3) _phase = "Solid" ;
			}
			if(_sim_type > 1){
				input >> _J;
				input >> _H;
				_equilibration = false;
			}
			if(_sim_type > 3){
				cerr << "PROBLEM: unknown simulation type" << endl;
				exit(EXIT_FAILURE);
			}
		} else if( property == "TEMP" ){
			input >> _temp;
			_beta = 1.0/_temp;
			string temp_string = to_string (_temp) ;
			temp_string = temp_string.substr(0, 4);
			if(_sim_type==0 and _equilibration){
				_output_dir = "../" + _phase + "/EQUILIBRATION/" + temp_string + "_" ;
			} else if (_sim_type==0) {
				_output_dir = "../" + _phase +"/OUTPUT/";
			} else if (_sim_type == 2){
				_output_dir = "../_OUTPUT/METROPOLIS/" + to_string(_temp).substr(0, 4) + "_" + to_string(_H).substr(0,4) + "_" ;
			} else if (_sim_type == 3){
				_output_dir = "../_OUTPUT/GIBBS/" + to_string(_temp).substr(0, 4) + "_" + to_string(_H).substr(0,4) + "_" ;
			}
			coutf.open(_output_dir + "output.dat");
			if(_sim_type == 0)		coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"	<< endl;
			else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"			<< endl;
			else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
			else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
			couta.open(_output_dir + "acceptance.dat"); // Set the heading line in file _OUTPUT/acceptance.dat
			couta << "#	 N_BLOCK:	ACCEPTANCE:" << endl;
			couta.close();
			coutf << "TEMPERATURE = " << _temp << endl;
		} else if( property == "RESTART" ){
			input >> _restart;
		} else if( property == "NPART" ){
			input >> _npart;
			_fx.resize(_npart);
			_fy.resize(_npart);
			_fz.resize(_npart);
			_particle.set_size(_npart);
			for(int i=0; i<_npart; i++){ 
				_particle(i).initialize();
				if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
			}
			coutf << "NPART = " << _npart << endl;
		} else if( property == "RHO" ){
			input >> _rho;
			_volume = _npart/_rho;
			_side.resize(_ndim);
			_halfside.resize(_ndim);
			double side = pow(_volume, 1.0/3.0);
			for(int i=0; i<_ndim; i++) _side(i) = side;
			_halfside=0.5*_side;
			coutf << "SIDE = ";
			for(int i=0; i<_ndim; i++){
				coutf << setw(12) << _side[i];
			}
			coutf << endl;
		} else if( property == "R_CUT" ){
			input >> _r_cut;
			coutf << "R_CUT = " << _r_cut << endl;
		} else if( property == "DELTA" ){
			input >> delta;
			coutf << "DELTA = " << delta << endl;
			_delta = delta;
		} else if( property == "NBLOCKS" ){
			input >> _nblocks;
			coutf << "NBLOCKS = " << _nblocks << endl;
		} else if( property == "NSTEPS" ){
			input >> _nsteps;
			coutf << "NSTEPS = " << _nsteps << endl;
		} else if( property == "ENDINPUT" ){
			coutf << "Reading input completed!" << endl;
			break;
		} else cerr << "PROBLEM: unknown input" << endl;
	}
	input.close();
	this->read_configuration();
	this->initialize_velocities();
	this->initialize_spin();

	coutf << "System initialized!" << endl;
	coutf.close();
	return;
}

void System :: initialize_spin(){
	double spin;
	ifstream cins;
	if(!_restart and _sim_type > 1){
		for(int i=0; i<_npart; i++){
			spin = _rnd.Rannyu();
			if (spin < 0.5) _particle(i).setspin(-1);
			else _particle(i).setspin(1);
		}
	}
	else if(_restart and _sim_type > 1){
		cins.open("../_INPUT/CONFIG/config.spin");
		for(int i=0; i<_npart; i++){
			cins >> spin;
			_particle(i).setspin(spin);
		}
		cins.close();
	}
	return;
}

void System :: initialize_velocities(){
	if(_restart and _sim_type==0){
		ifstream cinf;
		cinf.open("../_INPUT/CONFIG/velocities.in");
		if(cinf.is_open()){
			double vx, vy, vz;
			for(int i=0; i<_npart; i++){
				cinf >> vx >> vy >> vz;
				_particle(i).setvelocity(0,vx);
				_particle(i).setvelocity(1,vy);
				_particle(i).setvelocity(2,vz);
			}
		} else cerr << "PROBLEM: Unable to open INPUT file velocities.in"<< endl;
		cinf.close();
	} else {
		vec vx(_npart), vy(_npart), vz(_npart);
		vec sumv(_ndim);
		sumv.zeros();
		for (int i=0; i<_npart; i++){
			vx(i) = _rnd.Gauss(0.,sqrt(_temp));	// velocities in reduced units inizialized as sqrt(temperature)	
			vy(i) = _rnd.Gauss(0.,sqrt(_temp));
			vz(i) = _rnd.Gauss(0.,sqrt(_temp));
			sumv(0) += vx(i);
			sumv(1) += vy(i);
			sumv(2) += vz(i);
		}
		for (int idim=0; idim<_ndim; idim++) sumv(idim) = sumv(idim)/double(_npart);

		double sumv2 = 0.0, scalef;
		for (int i=0; i<_npart; i++){
			vx(i) = vx(i) - sumv(0);
			vy(i) = vy(i) - sumv(1);
			vz(i) = vz(i) - sumv(2);
			sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
		}
		sumv2 /= double(_npart);
		scalef = sqrt(3.0 * _temp / sumv2);	 // velocity scale factor 
		for (int i=0; i<_npart; i++){
			_particle(i).setvelocity(0, vx(i)*scalef);
			_particle(i).setvelocity(1, vy(i)*scalef);
			_particle(i).setvelocity(2, vz(i)*scalef);
		}
	}
	if(_sim_type == 0){
	double xold, yold, zold;
		for (int i=0; i<_npart; i++){
			xold = this->pbc( _particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
			yold = this->pbc( _particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
			zold = this->pbc( _particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);
			_particle(i).setpositold(0, xold);
			_particle(i).setpositold(1, yold);
			_particle(i).setpositold(2, zold);
		}
	}
	return;
}

void System :: initialize_properties(){ // Initialize data members used for measurement of properties

	string property;
	int index_property = 0;
	_nprop = 0;

	_measure_penergy 	= false; //Defining which properties will be measured
	_measure_kenergy 	= false;
	_measure_tenergy 	= false;
	_measure_pressure = false;
	_measure_gofr 		= false;
	_measure_magnet 	= false;
	_measure_cv 		= false;
	_measure_chi 		= false;

	ifstream input("../_INPUT/properties.dat");
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "POTENTIAL_ENERGY" ){ 
				if(!_equilibration){
					ofstream coutp(_output_dir + "potential_energy.dat");
					coutp << "#		 BLOCK:	ACTUAL_PE:		 PE_AVE:			ERROR:" << endl;
					coutp.close();
				}
					_nprop++;
					_index_penergy = index_property;
					_measure_penergy = true;
					index_property++;
					_vtail = 0.0; // TO BE FIXED IN EXERCISE 7

			} else if( property == "KINETIC_ENERGY" ){
				if(!_equilibration){
					ofstream coutk(_output_dir + "kinetic_energy.dat");
					coutk << "#		 BLOCK:	 ACTUAL_KE:		KE_AVE:			ERROR:" << endl;
					coutk.close();
				}
					_nprop++;
					_measure_kenergy = true;
					_index_kenergy = index_property;
					index_property++;
			} else if( property == "TOTAL_ENERGY" ){
				if(!_equilibration and _sim_type!=2 and _sim_type!=3){
					ofstream coutt(_output_dir + "total_energy.dat");
					coutt << "#		 BLOCK:	 ACTUAL_TE:		TE_AVE:			ERROR:" << endl;
					coutt.close();
				}
				if(_sim_type == 2){
					ifstream check ("../_OUTPUT/METROPOLIS/total_energy_" + to_string(_H).substr(0,4) + ".dat");				 
					if (check.peek() == std::ifstream::traits_type::eof()){
						ofstream coutpr("../_OUTPUT/METROPOLIS/total_energy_" + to_string(_H).substr(0,4) + ".dat");
						coutpr << "TEMPERATURE:	 	FINAL_TE:	ERROR:" << endl;
						coutpr.close();
					}
				} else if(_sim_type == 3){
					ifstream check ("../_OUTPUT/GIBBS/total_energy_" + to_string(_H).substr(0,4) + ".dat");				 
					if (check.peek() == std::ifstream::traits_type::eof()){
						ofstream coutpr("../_OUTPUT/GIBBS/total_energy_" + to_string(_H).substr(0,4) + ".dat");
						coutpr << "TEMPERATURE:	 	FINAL_TE:	ERROR:" << endl;
						coutpr.close();
					}
				}
					_nprop++;
					_measure_tenergy = true;
					_index_tenergy = index_property;
					index_property++;

			} else if( property == "TEMPERATURE"){
				ofstream coutte(_output_dir + "temperature.dat");
				coutte << "#		 BLOCK:	 ACTUAL_T:		 T_AVE:			 ERROR:" << endl;
				coutte.close();
				_nprop++;
				_measure_temp = true;
				_index_temp = index_property;
				index_property++;
			} else if( property == "PRESSURE" ){
				if(!_equilibration and _sim_type!=2 and _sim_type!=3){
					ofstream coutpr(_output_dir + "pressure.dat");
					coutpr << "#		 BLOCK:	 ACTUAL_P:		 P_AVE:			 ERROR:" << endl;
					coutpr.close();
				}
					_nprop++;
					_measure_pressure = true;
					_index_pressure = index_property;
					index_property++;
					_ptail = 0.0; // TO BE FIXED IN EXERCISE 7
				 
			} else if( property == "GOFR" ){
				if(!_equilibration and _sim_type!=2 and _sim_type!=3){
					ofstream coutgr(_output_dir + "gofr.dat");
					coutgr << "# DISTANCE:		AVE_GOFR:	 ERROR:" << endl;
					coutgr.close();
				}
					input>>_n_bins;
					_nprop+=_n_bins;
					_bin_size = (_halfside.min() )/(double)_n_bins;
					_measure_gofr = true;
					_index_gofr = index_property;
					index_property+= _n_bins;
				 
			} else if( property == "MAGNETIZATION" ){
				if(_H!=0){
					if(!_equilibration and _sim_type!=2 and _sim_type!=3){
						ofstream coutpr(_output_dir + "magnetization.dat");
						coutpr << "#		 BLOCK:	 ACTUAL_M:		 M_AVE:			 ERROR:" << endl;
						coutpr.close();
					}
					if(_sim_type == 2){
						ifstream check ("../_OUTPUT/METROPOLIS/magnetization_" + to_string(_H).substr(0,4) + ".dat");				 
						if (check.peek() == std::ifstream::traits_type::eof()){
							ofstream coutpr("../_OUTPUT/METROPOLIS/magnetization_" + to_string(_H).substr(0,4) + ".dat",ios::app);
							coutpr << "TEMPERATURE:	 	FINAL_M:		ERROR:" << endl;
							coutpr.close();
						}
					} else if(_sim_type == 3){
						ifstream check ("../_OUTPUT/GIBBS/magnetization_" + to_string(_H).substr(0,4) + ".dat");				 
						if (check.peek() == std::ifstream::traits_type::eof()){
							ofstream coutpr("../_OUTPUT/GIBBS/magnetization_" + to_string(_H).substr(0,4) + ".dat",ios::app);
							coutpr << "TEMPERATURE:	 	FINAL_M:		ERROR:" << endl;
							coutpr.close();
						}
					}
						_nprop++;
						_measure_magnet = true;
						_index_magnet = index_property;
						index_property++;
				}
			} else if( property == "SPECIFIC_HEAT" ){
				if(_H==0){
					if(!_equilibration and _sim_type!=2 and _sim_type!=3){
						ofstream coutpr(_output_dir + "specific_heat.dat");
						coutpr << "#		 BLOCK:	 ACTUAL_CV:		CV_AVE:			ERROR:" << endl;
						coutpr.close();
					}
					if(_sim_type == 2){
						ifstream check ("../_OUTPUT/METROPOLIS/specific_heat_" + to_string(_H).substr(0,4) + ".dat");				 
						if (check.peek() == std::ifstream::traits_type::eof()){
							ofstream coutpr("../_OUTPUT/METROPOLIS/specific_heat_" + to_string(_H).substr(0,4) + ".dat",ios::app);
							coutpr << "TEMPERATURE:	 	FINAL_M:		ERROR:" << endl;
							coutpr.close();
						}
					} else if(_sim_type == 3){
						ifstream check ("../_OUTPUT/GIBBS/specific_heat_" + to_string(_H).substr(0,4) + ".dat");				 
						if (check.peek() == std::ifstream::traits_type::eof()){
							ofstream coutpr("../_OUTPUT/GIBBS/specific_heat_" + to_string(_H).substr(0,4) + ".dat",ios::app);
							coutpr << "TEMPERATURE:	 	FINAL_M:		ERROR:" << endl;
							coutpr.close();
						}
					}
						_nprop++;
						_measure_cv = true;
						_index_cv = index_property;
						index_property++;
				}
			} else if( property == "SUSCEPTIBILITY" ){
				if(_H==0){
					if(!_equilibration and _sim_type!=2 and _sim_type!=3){
						ofstream coutpr(_output_dir + "susceptibility.dat");
						coutpr << "#		 BLOCK:	 ACTUAL_X:		 X_AVE:			 ERROR:" << endl;
						coutpr.close();
					}
					if(_sim_type == 2){
						ifstream check ("../_OUTPUT/METROPOLIS/susceptibility_" + to_string(_H).substr(0,4) + ".dat");				 
						if (check.peek() == std::ifstream::traits_type::eof()){
							ofstream coutpr("../_OUTPUT/METROPOLIS/susceptibility_" + to_string(_H).substr(0,4) + ".dat",ios::app);
							coutpr << "TEMPERATURE:	 	FINAL_M:		ERROR:" << endl;
							coutpr.close();
						}
					} else if(_sim_type == 3){
						ifstream check ("../_OUTPUT/GIBBS/susceptibility_" + to_string(_H).substr(0,4) + ".dat");				 
						if (check.peek() == std::ifstream::traits_type::eof()){
							ofstream coutpr("../_OUTPUT/GIBBS/susceptibility_" + to_string(_H).substr(0,4) + ".dat",ios::app);
							coutpr << "TEMPERATURE:	 	FINAL_M:		ERROR:" << endl;
							coutpr.close();
						}
					}
						_nprop++;
						_measure_chi = true;
						_index_chi = index_property;
						index_property++;
				}
			} else if( property == "ENDPROPERTIES"){
				ofstream coutf;
				coutf.open(_output_dir + "output.dat",ios::app);
				coutf << "Reading properties completed!" << endl;
				coutf.close();
				break;
			} else cerr << "PROBLEM: unknown property" << endl;
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open properties.dat" << endl;

	// according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
	_measurement.resize(_nprop);
	_average.resize(_nprop);
	_block_av.resize(_nprop);
	_global_av.resize(_nprop);
	_global_av2.resize(_nprop);
	_average.zeros();
	_global_av.zeros();
	_global_av2.zeros();
	_nattempts = 0;
	_naccepted = 0;
	return;
}

void System :: save_final_block (int blk) {
	// for Gibbs simulation, save final block data in a cumulative data file
	double error;
	if (_sim_type > 1 and blk == _nblocks) {
		if (_measure_tenergy) {
			ofstream coutte;
			if(_sim_type==2) coutte.open("../_OUTPUT/METROPOLIS/total_energy_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if(_sim_type==3) coutte.open("../_OUTPUT/GIBBS/total_energy_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if (coutte.is_open()){
				error = this->error(_global_av(_index_tenergy),_global_av2(_index_tenergy), this->get_nbl()) ;
				// cout << _temp << " " <<  _global_av(_index_tenergy) << " " << error << endl;
				coutte << setw(12) << _temp
					<< setw(12) << _global_av(_index_tenergy)/double(blk)
					<< setw(12) << error << endl;
			} else cerr << " Error: unable to open file" << endl;
			coutte.close();
		}
		if (_measure_chi and _H==0) {
			ofstream coutchi ;
			if(_sim_type==2) coutchi.open("../_OUTPUT/METROPOLIS/susceptibility_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if(_sim_type==3) coutchi.open("../_OUTPUT/GIBBS/susceptibility_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if (coutchi.is_open()){
				error = this->error(_global_av(_index_chi),_global_av2(_index_chi), this->get_nbl()) ;
				// cout << _temp << " " <<  _global_av(_index_chi) << " " << error << endl;
				coutchi << setw(12) << _temp
						<< setw(12) << _global_av(_index_chi)/double(blk)
						<< setw(12) << error << endl;
			} else cerr << " Error: unable to open file" << endl;
			coutchi.close();
		}
		if (_measure_cv and _H==0) {
			ofstream coutcv;
			if(_sim_type==2) coutcv.open("../_OUTPUT/METROPOLIS/specific_heat_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if(_sim_type==3) coutcv.open("../_OUTPUT/GIBBS/specific_heat_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if (coutcv.is_open()){
				error = this->error(_global_av(_index_cv),_global_av2(_index_cv), this->get_nbl()) ;
				// cout << _temp << " " <<  _global_av(_index_cv) << " " << error << endl;
				coutcv << setw(12) << _temp
						<< setw(12) << _global_av(_index_cv)/double(blk)
						<< setw(12) << error << endl;
			} else cerr << " Error: unable to open file" << endl;
			coutcv.close();
		}
		if (_measure_magnet and _H!=0) {
			ofstream coutmag;
			if(_sim_type==2) coutmag.open("../_OUTPUT/METROPOLIS/magnetization_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if(_sim_type==3) coutmag.open("../_OUTPUT/GIBBS/magnetization_" + to_string(_H).substr(0,4) + ".dat",ios::app);
			if (coutmag.is_open()){
				error = this->error(_global_av(_index_magnet),_global_av2(_index_magnet), this->get_nbl()) ;
				// cout << _temp << " " <<  _global_av(_index_magnet) << " " << error << endl;
				coutmag << setw(12) << _temp
						<< setw(12) << _global_av(_index_magnet)/double(blk)
						<< setw(12) << error << endl;
			} else cerr << " Error: unable to open file" << endl;
			coutmag.close();
		}
	}
}

void System :: finalize(){
	ofstream coutf;
	coutf.open(_output_dir + "output.dat",ios::app);
	if (!_equilibration){
		coutf << "Writing final configurations" << endl;
		this->write_configuration();
	}
	_rnd.SaveSeed();
	coutf << "Simulation completed!" << endl;
	coutf.close();
	return;
}

// Write current configuration as a .xyz file in directory ../CONFIG/
void System :: write_configuration(){
	ofstream coutf;
	if(!_equilibration and _sim_type < 2){
		coutf.open(_output_dir + "/CONFIG/config.xyz");
		if(coutf.is_open()){
			coutf << _npart << endl;
			coutf << "#Comment!" << endl;
			for(int i=0; i<_npart; i++){
				coutf << "LJ" << "	" 
							<< setw(12) << _particle(i).getposition(0,true)/_side(0)					// x
							<< setw(12) << _particle(i).getposition(1,true)/_side(1)					// y
							<< setw(12) << _particle(i).getposition(2,true)/_side(2) << endl; // z
			}
		} else cerr << "PROBLEM: Unable to open config.xyz" << endl;
		coutf.close();
		this->write_velocities();
	} else {
		coutf.open("../_INPUT/CONFIG/config.spin");
		for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
		coutf.close();
	}
	return;
}

// Write configuration nconf as a .xyz file in directory ../_phase/CONFIG/
void System :: write_XYZ(int nconf){
	ofstream coutf;
	coutf.open(_output_dir + "/CONFIG/config_" + to_string(nconf) + ".xyz");
	if(coutf.is_open()){
		coutf << _npart << endl;
		coutf << "#Comment!" << endl;
		for(int i=0; i<_npart; i++){
			coutf << "LJ" << "	" 
						<< setw(12) << _particle(i).getposition(0,true)					// x
						<< setw(12) << _particle(i).getposition(1,true)					// y
						<< setw(12) << _particle(i).getposition(2,true) << endl; // z
		}
	} else cerr << "PROBLEM: Unable to open config.xyz" << endl;
	coutf.close();
	return;
}

void System :: write_velocities(){
	ofstream coutf;
	coutf.open(_output_dir + "/CONFIG/velocities.out");
	if(coutf.is_open()){
		for(int i=0; i<_npart; i++){
			coutf << setw(12) << _particle(i).getvelocity(0)					// vx
						<< setw(12) << _particle(i).getvelocity(1)					// vy
						<< setw(12) << _particle(i).getvelocity(2) << endl; // vz
		}
	} else cerr << "PROBLEM: Unable to open velocities.dat" << endl;
	coutf.close();
	return;
}

// Read configuration from a .xyz file in directory
void System :: read_configuration(){
	ifstream cinf;
	string comment;
	string particle;
	double x, y, z;
	int ncoord;
	int spin;
	if(_sim_type < 2) cinf.open("../_INPUT/CONFIG/config.xyz");
	else cinf.open("../_INPUT/CONFIG/config.ising");
	if(cinf.is_open()){
		cinf >> ncoord;
		if (ncoord != _npart){
			cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
			exit(EXIT_FAILURE);
		}
		cinf >> comment;
		for(int i=0; i<_npart; i++){
			cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
			_particle(i).setposition(0, this->pbc(_side(0)*x, 0));
			_particle(i).setposition(1, this->pbc(_side(1)*y, 1));
			_particle(i).setposition(2, this->pbc(_side(2)*z, 2));
			_particle(i).acceptmove(); // _x_old = _x_new
		}
	} else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
	cinf.close();	
	return;
}

void System :: block_reset(int blk){ // Reset block accumulators to zero
	ofstream coutf;
	if(blk>0){
		coutf.open(_output_dir + "output.dat",ios::app);
		coutf << "Block completed: " << blk << endl;
		coutf.close();
	}
	_block_av.zeros();
	return;
}

void System :: measure(){ // Measure properties
	_measurement.zeros();
	// POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
	int bin;
	vec distance;
	distance.resize(_ndim);
	double penergy_temp = 0., dr; // temporary accumulator for potential energy
	double kenergy_temp = 0.; // temporary accumulator for kinetic energy
	double tenergy_temp = 0.; // temporary accumulator for total (kin.+pot.) energy
	double pressure_temp = 0.;
	double cv_temp = 0.;
	double magnetization = 0.;
	double virial = 0.;
	double spin_sum = 0.;

	if (_measure_penergy or _measure_pressure or _measure_gofr) {
		for (int i=0; i<_npart-1; i++){
			for (int j=i+1; j<_npart; j++){
				distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
				distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
				distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
				dr = sqrt( arma::dot(distance,distance) );
				// GOFR ... TO BE FIXED IN EXERCISE 7
				if(dr < _r_cut){
					if(_measure_penergy)	penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6) ; // POTENTIAL ENERGY
					if(_measure_pressure) pressure_temp += 1.0/pow(dr,12) - 1.0/2.0 * 1.0/pow(dr,6) ; // PRESSURE
				}
			}
		}
	}
	// POTENTIAL ENERGY //////////////////////////////////////////////////////////
	if (_measure_penergy){
		penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
		_measurement(_index_penergy) = penergy_temp;
	}
	// KINETIC ENERGY ////////////////////////////////////////////////////////////
	if (_measure_kenergy){
		for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * arma::dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
		kenergy_temp /= double(_npart);
		_measurement(_index_kenergy) = kenergy_temp;
	}
	// TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
	if (_measure_tenergy){
		if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp; // Molecular Dynamics
		else { 	// Ising
			double s_i, s_j;
			for (int i=0; i<_npart; i++){
				 s_i = double(_particle(i).getspin());
				 s_j = double(_particle(this->pbc(i+1)).getspin());
				 tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
			}
			tenergy_temp /= double(_npart);		// energy per particle
			_measurement(_index_tenergy) = tenergy_temp;
		}
	}
	// TEMPERATURE ///////////////////////////////////////////////////////////////
	if (_measure_temp and _measure_kenergy) 
		_measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
	// PRESSURE //////////////////////////////////////////////////////////////////
	if (_measure_pressure and _measure_temp and _measure_kenergy){
		pressure_temp = 48./3. * 1./_volume * pressure_temp ;
		_measurement(_index_pressure) = _rho * _measurement(_index_temp) + pressure_temp ;
	}
	// TOTAL SPIN ////////////////////////////////////////////////////////////////
	if (_sim_type > 1)
		for(int i = 0; i<_npart; i++) 
			spin_sum += double(_particle(i).getspin());	
	// MAGNETIZATION /////////////////////////////////////////////////////////////
	if (_measure_magnet)
		_measurement(_index_magnet) = spin_sum ;
	// SPECIFIC HEAT /////////////////////////////////////////////////////////////
	if(_measure_cv and _measure_tenergy)
		_measurement(_index_cv) = pow(_measurement(_index_tenergy)*_npart,2) ; // total energy squared
	// SUSCEPTIBILITY //////////////////////////////////////////////////////////// 
	if(_measure_chi)
		_measurement(_index_chi) = _beta * pow(spin_sum,2);
	
	_block_av += _measurement; //Update block accumulators

	return;
}

void System :: averages(int blk){

	ofstream coutf;
	double average, sum_average, sum_ave2;

	_average = _block_av / double(_nsteps);
	
	if (_measure_cv){
		_average(_index_cv) -= pow((_average(_index_tenergy)*_npart),2) ;
		_average(_index_cv) *= pow(_beta,2) ;	// re multiply indexcv
	}
	
	_global_av	+= _average;
	_global_av2 += _average % _average; // % -> element-wise multiplication

	// POTENTIAL ENERGY //////////////////////////////////////////////////////////
	if (_measure_penergy and !_equilibration){
		coutf.open(_output_dir + "potential_energy.dat",ios::app);
		average	= _average(_index_penergy);
		sum_average = _global_av(_index_penergy);
		sum_ave2 = _global_av2(_index_penergy);
		coutf << setw(12) << blk 
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	// KINETIC ENERGY ////////////////////////////////////////////////////////////
	if (_measure_kenergy and !_equilibration){
		coutf.open(_output_dir + "kinetic_energy.dat",ios::app);
		average	= _average(_index_kenergy);
		sum_average = _global_av(_index_kenergy);
		sum_ave2 = _global_av2(_index_kenergy);
		coutf << setw(12) << blk
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	// TOTAL ENERGY //////////////////////////////////////////////////////////////
	if (_measure_tenergy and !_equilibration and _sim_type!=2 and _sim_type!=3){
		coutf.open(_output_dir + "total_energy.dat",ios::app);
		average	= _average(_index_tenergy);
		sum_average = _global_av(_index_tenergy);
		sum_ave2 = _global_av2(_index_tenergy);
		coutf << setw(12) << blk
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	// TEMPERATURE ///////////////////////////////////////////////////////////////
	if (_measure_temp){
		coutf.open(_output_dir + "temperature.dat",ios::app);
		average	= _average(_index_temp);
		sum_average = _global_av(_index_temp);
		sum_ave2 = _global_av2(_index_temp);
		coutf << setw(12) << blk
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	// PRESSURE //////////////////////////////////////////////////////////////////
	if (_measure_pressure and !_equilibration and _sim_type!=2 and _sim_type!=3){
		coutf.open(_output_dir + "pressure.dat",ios::app);
		average	= _average(_index_pressure);
		sum_average = _global_av(_index_pressure);
		sum_ave2 = _global_av2(_index_pressure);
		coutf << setw(12) << blk
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	
	// GOFR //////////////////////////////////////////////////////////////////////
	// TO BE FIXED IN EXERCISE 7
	
	// MAGNETIZATION /////////////////////////////////////////////////////////////
	if (_measure_magnet and !_equilibration and _sim_type!=2 and _sim_type!=3){
		coutf.open(_output_dir + "magnetization.dat",ios::app);
		average	= _average(_index_magnet);
		sum_average = _global_av(_index_magnet);
		sum_ave2 = _global_av2(_index_magnet);
		coutf << setw(12) << blk 
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	// SPECIFIC HEAT /////////////////////////////////////////////////////////////
	if (_measure_cv and !_equilibration and _sim_type!=2 and _sim_type!=3){
		coutf.open(_output_dir + "specific_heat.dat",ios::app);
		average	= _average(_index_cv);
		sum_average = _global_av(_index_cv);
		sum_ave2 = _global_av2(_index_cv);
		coutf << setw(12) << blk
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}
	// SUSCEPTIBILITY ////////////////////////////////////////////////////////////
	if (_measure_chi and !_equilibration and _sim_type!=2 and _sim_type!=3){
		coutf.open(_output_dir + "susceptibility.dat",ios::app);
		average	= _average(_index_chi);
		sum_average = _global_av(_index_chi);
		sum_ave2 = _global_av2(_index_chi);
		coutf << setw(12) << blk 
					<< setw(12) << average
					<< setw(12) << sum_average/double(blk)
					<< setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
		coutf.close();
	}

	// ACCEPTANCE ////////////////////////////////////////////////////////////////
	double fraction;
	coutf.open(_output_dir + "acceptance.dat",ios::app);
	if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
	else fraction = 0.0; 
	coutf << setw(12) << blk << setw(12) << fraction << endl;
	coutf.close();

	return;
}

double System :: error(double acc, double acc2, int blk){
	if(blk <= 1) return 0.0;
	else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
	return _nblocks;
}

int System :: get_nsteps(){
	return _nsteps;
}

int System :: get_phase(){
	return _int_phase;
}

bool System :: get_equilib () {
	return _equilibration;
}


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/