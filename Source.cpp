#include <iostream>
#include <vector>
#include <mpi.h>

#define SIZE 1000				// Grid size (1000x1000)
#define STEPS 100				// Number of simulation steps
#define DIFFUSION_COEFF 0.1		// Diffusion coefficient
#define DT 0.01					// Time step

using namespace std;

typedef vector<vector<double>> Matrix;

// Initialize the grid with some initial conditions (e.g., a concentration source in the middle)
void initialize(Matrix& grid, int rank, int rows_per_process) {
	int global_mid = SIZE / 2;
	int local_mid = global_mid - rank * rows_per_process;

	if (local_mid >= 0 && local_mid < rows_per_process) {
		grid[local_mid][SIZE / 2] = 100.0;	// Initial concentration at the center
	}

}

// // Display the grid (for small grid sizes only, not practical for large grids)
void display(const Matrix& grid) {
	for (int i = 0; i < grid.size(); ++i) {
		for (int j = 0; j < grid[i].size(); ++j) {
			cout << grid[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// Perform one step of the similation (calculate the new velocities)
void step(const Matrix& current, Matrix& next, int rank, int size) {
	int rows_per_process = current.size();

	for (int i = 1; i < rows_per_process - 1; ++i) {
		for (int j = 1; j < SIZE - 1; ++j) {
			next[i][j] = current[i][j] + DIFFUSION_COEFF * DT * (
				current[i + 1][j] + current[i - 1][j] +
				current[i][j + 1] + current[i][j - 1] - 4 * current[i][j]
				);
		}
	}

	// Exchange boundary rows with neighboring processes
	if (rank > 0) {
		// Send top boundary row to process above
		MPI_Send(current[1].data(), SIZE, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		// Receive bottom boundary row from the process above
		MPI_Recv(current[0].data(), SIZE, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if (rank < size - 1) {
		// Send bottom boundary row to process below
		MPI_Send(current[rows_per_process - 2].data(), SIZE, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
		// Send top boundary row from the process below
		MPI_Recv(current[rows_per_process - 1].data(), SIZE, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

}

// Run the simulation for a certain number of steps
void run_simulation(Matrix& grid, int rank, int size) {
	int rows_per_process = grid.size();
	Matrix next_grid = grid;

	for (int step_count = 0; step_count < STEPS; ++step_count) {
		step(grid, next_grid, rank, size);
		grid.swap(next_grid);	// Swap grids for the next iteration
		if (rank == 0 && step_count % 10 == 0) {
			cout << "Completed step: " << step_count << endl;
		}
	}
}


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Divide the grid into chunks for each process
	int rows_per_process = SIZE / size + 2;		// Add 2 extra rows for the boundary exchange
	Matrix grid(rows_per_process, vector<double>(SIZE, 0.0));

	// Initialize the grid
	initialize(grid, rank, rows_per_process);

	// Run the simulation
	double start_time = MPI_Wtime();
	run_simulation(grid, rank, size);
	double end_time = MPI_Wtime();

	// Display final results (optional for smaill grids)
	// If (rank == 0) display(grid);

	if (rank == 0) {
		cout << "Simulation completed in " << (end_time - start_time) << " seconds." << endl;
	}

	MPI_Finalize();
	return 0;

}




