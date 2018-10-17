/*
 * mpi2norm_timer.cpp
 * Description: Timing the MPI 2-norm computation.
 * @author Johnny Sellers
 * @version 0.1 05/30/2017
 */
#include <iostream>
#include <string>
#include <cmath>
#include <mpi.h>
#include "Vector.hpp"

double mpiTwoNorm(const Vector& v);

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (argc != 2) {
		if (0 == world_rank) std::cout << "Usage: " << argv[0] << " [ vectorSize ]" << std::endl;
		return -1;
	}

	size_t vectorSize = std::stol(argv[1]);
	Vector x(vectorSize);
	double global, mpi, seq_time;

	if (0 == world_rank) {
		Vector v(vectorSize*world_size);
		randomize(v);

		double seq_start, seq_finish;
		seq_start = MPI_Wtime();
		global = twoNorm(v);
		seq_finish = MPI_Wtime();
		seq_time = seq_finish - seq_start;

		MPI_Scatter(&v(world_rank*vectorSize), vectorSize, MPI_DOUBLE, &x(0), vectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Scatter(NULL, vectorSize, MPI_DOUBLE, &x(0), vectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	double local_start, local_finish, local_elapsed, para_time;
	MPI_Barrier(MPI_COMM_WORLD);
	local_start = MPI_Wtime();
	mpi = mpiTwoNorm(x);
	local_finish = MPI_Wtime();
	local_elapsed = local_finish - local_start;

	MPI_Reduce(&local_elapsed, &para_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (0 == world_rank) {
		std::cout << "#\tgvSize\tseq_time\tpara_time\tspeedup\tdiff" << std::endl;
		std::cout << world_size << "\t" << vectorSize*world_size << "\t" << seq_time << "\t" << para_time << "\t" << seq_time/para_time << "\t" << std::abs(global-mpi) << std::endl;
	}

  MPI_Finalize();

  return 0;
}

double mpiTwoNorm(const Vector& v) {
	double sum = 0.0, mpi;
  for (int i = 0; i < v.numRows(); ++i)
    sum += v(i)*v(i);
	MPI_Allreduce(&sum, &mpi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return std::sqrt(mpi);
}
