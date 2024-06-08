## How to use my code - LSN Exercise 06

### Execution

The execution process is taken care of by the bash script `LSN_06/6.1/_INPUT/execution_routine.sh`. To use it, you have to:
- verify that the script is executable by going into `./_INPUT` and running `chmod +x execution_routine.sh`
- modify the input file `./_INPUT/input_orig.in` to execute Metropolis or Gibbs sampling by modiyfing the first data under `SIMULATION_TYPE`. The options are 2 for Metropolis, 3 for Gibbs.
- execute by running `./execution_routine.sh`
- deletion of old files and compiling before execution is taken care of by the script
