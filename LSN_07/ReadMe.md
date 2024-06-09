## How to use my code - LSN Exercise 07

All of the subsections of the exercise have been tackled in the same directory, namely `./7.1`.
For each state phase that the simulations can be run on, there is a dedicated directory, namely `./7.1/StatePhase` which contains both the input files and the results of the simulation. 

The source codes of the simulation are in the `./7.1/_SOURCE` directory. Head there and follow the guidelines below
- To compile hit `make`
- To remove `*.o` and `*.exe` files, hit `make clean`
- To execute the code, hit `./simulator.exe` followed by the path of the input file. For example, `./simulator.exe ../Gas/INPUT/input_MC_NVT.in`