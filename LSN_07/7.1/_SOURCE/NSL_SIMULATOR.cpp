/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){


	if (argc != 2){
		cerr << "Give path to input file to be executed" << endl;
	}
	string dir = argv[1];

   int nconf = 1;
   vector<int> equil_blocks_MD_NVE = {400,200,200};     	// Molecular Dynamics, depending on phase, number of blocks necessary for equilibration

   System SYS;
   
   SYS.initialize(dir);
	SYS.initialize_properties();
	SYS.block_reset(0);

   if (!SYS.get_equilib() and SYS.get_sim_type()==0){										// cycle through equilibration steps for different _sim_type
      for (int i=0; i<equil_blocks_MD_NVE[SYS.get_phase()-1]; i++)
      	for(int j=0; j < 1000; j++) //loop over steps in a block
         	SYS.step();
   }

	for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
		// cout << "looping blocks" << endl;
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
				// cout << "looping step " << j << endl;
         SYS.step();
				// cout << "step done " << j << endl;
         SYS.measure();
				// cout << "measure done " << j << endl;
         if(j%10 == 0){
            // SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
            nconf++;
         }
      }
      SYS.averages(i+1);
      // SYS.save_final_block(i+1);
		SYS.block_reset(i+1);
   }
   SYS.finalize();

   return 0;
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