========================================================================
    CONSOLE APPLICATION : Wave2D Project Overview
========================================================================
This Euler2D Project implements the Swept2D Rule decomposition for the 2D
Wave Equation.  Use the included compile script to compile the code.
The code is executed as follows:

mpirun -np NUM_PROC ./wave2d.exe NUM_SPACIAL_POINTS NUM_X_PARTS NUM_Y_PARTS NUM_SWEPT_CYCLES

Where:
NUM_PROC: number of MPI processes to run
NUM_SPACIAL_POINTS: the square root of the number of grid points in each MPI process (square partition side length)
NUM_X_PARTS: number of squares in the X direction
NUM_Y_PARTS:  number of squares in the Y direction
NUM_SWEPT_CYCLES: number of substeps to run

A typical run is as follow:

mpirun -np 4 ./wave2d.exe 8 2 2 100
