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
   vector<int> equil_blocks = {400,200,200};     // depending on phase, number of blocks necessary for equilibration

   System SYS;
   
   SYS.initialize(dir);
   SYS.initialize_properties();
   SYS.block_reset(0);

   if (!SYS.get_equilib())
	   for (int i=0; i<equil_blocks[SYS.get_phase()-1]; i++)
         for(int j=0; j < SYS.get_nsteps(); j++) //loop over steps in a block
            SYS.step();

   for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
         SYS.step();
         SYS.measure();
         if(j%10 == 0){
            // SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
            nconf++;
         }
      }
      SYS.averages(i+1);
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