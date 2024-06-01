## How to use my code

### Equilibration

The equilibration process is taken care of by the bash script `./_INPUT/equilibration_routine.sh` which requires a phase to be chosen for equilibration. Thus to execute:
- verify that the script is executable by going into `./_INPUT` and running `chmod +x equilibration_routine.sh`
- execute with a choice of phase with `./equilibration_routine.sh <phase>` where inputs can be {Gas, Liquid, Solid}. The compiling of the code is taken care of automatically by the bash script

### Simulation

The compilation of the simulation itself is slightly different. Go into `./_SOURCE` and
- to compile hit `make`
- to remove `*.o` and `*.exe` files, hit `make clean`
- to remove all data produced from previous simulations, hit `make remove`
- finally to execute the simulation, you have to provide the input file, hit `./simulator.exe ../<phase>/INPUT/input.in` where again phases are {Gas, Liquid, Solid} 
