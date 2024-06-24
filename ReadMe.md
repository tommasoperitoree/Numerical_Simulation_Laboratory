# Numerical Simulation Laboratory
## Tommaso Peritore

This is a collection of exercises written during Laboratorio di Simulazione Numerica (Numerical Simulation Laboratory) (A.A. 2023-2024) offered at the Physics Department of Università degli Studi di Milano Statale. 

## Prerequisites
In order to execute the commands described in this file you will need the following tools: make, gcc, mpi library. 

On a Debian-based system you can install them with:
```bash
sudo apt-get install make gcc opnmpi-bin
````
or on Arm-based systems with homebrew:
```bash
brew install make gcc openmpi
```

## Compilation and Execution
To compile the exercises specific instructions have been given in each subfolder dedicated to the single exercise.

In general, each subfolder contains a Makefile that can be used to compile the code. Some variations of the following instructions are used:
```bash
make
make run
```
In any case, it is necessary to enter the subfolder of the exercise to execute the commands.	

## Jupyter Notebooks
In each subdirectory, there is a Jupyter Notebook that contains the results of the exercises. The results are presented in a clear and readable way, with the help of graphs and tables.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.