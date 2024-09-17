MPI Initialization:
I initialize MPI with MPI_Init and obtain the rank (ID of the process) and the total number of processes using MPI_Comm_rank and MPI_Comm_size.
The grid is split among processes, and each process handles its own portion of the grid.

Grid Initialization:
Each process initializes its part of the grid. The central point of the grid (i.e., the source of fluid concentration) is placed in the middle of the grid in the appropriate process.

Parallel Computation with MPI:
Each process updates its portion of the grid independently, using a simple diffusion equation.
Boundary rows are exchanged between processes using MPI_Send and MPI_Recv to ensure neighboring rows are updated correctly.

Grid Updates and Communication:
The simulation runs in a loop, with each process calculating the new values for its portion of the grid and exchanging boundary rows with its neighboring processes to ensure consistent values at the borders.

Timing the Simulation:
I use MPI_Wtime() to time the simulation and calculate the total runtime.

Final Results:
After completing the simulation, the root process (rank == 0) outputs the total runtime. Optionally, we could gather the results from all processes and display the final grid.
