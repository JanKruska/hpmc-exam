# hpmc-exam
Source Code for my end of semester project in the High Performance Matrix Computations course offered at RWTH Aachen in WS21/22. 

# Structure
## Directories
* `data` contains any produced data outputs, i.e. diagnostics produced from jobs.
* `jobs` contains SLURM jobs to be run on an HPC cluster, tested on teh RWTH CLAIX_18 partition.
* `plots` contains code to produce plots from the produced data.
* `report` contains LaTeX source for the report that was handed in.

## Source
* `main.c` The main program containing the procedure to solve the given problem, as well as surrounding measurement code
* `utils.c` Some helper functions for creating matrices and calculating some statistical measurements
* `dgemm.c` Determines MKL version, and performs a simple dgemm on two 1000x1000 matrices to get a reference for it's GFLOPS performed.

# Compiling
The code was compiled using gcc version 10.1.0 (GCC) and requires a working Intel MKL version to be installed properly.
To check if MKL is installed correctly check whether the ${MKLROOT} environment variable is set.
Furthermore the code requires OpenMP.