/**
		mpi2norm_timer.cpp, AMATH 583
		Description: mpi2norm_timer.cpp for PS8

		@author Johnny Sellers
		@version 0.1 05/30/2017
*/
#include <iostream>
#include <string>
#include <cmath>
#include <mpi.h>
#include "Vector.hpp"

double mpiTwoNorm(const Vector& v);

int main(int argc, char* argv[]) {
  MPI::Init();

  int myrank = MPI::COMM_WORLD.Get_rank();
  int mysize = MPI::COMM_WORLD.Get_size();

	if (argc != 2) {
		if (0 == myrank) std::cout << "Usage: " << argv[0] << " [ vectorSize ]" << std::endl;
		return -1;
	}

	size_t vectorSize = std::stol(argv[1]);
	Vector x(vectorSize);
	double global, mpi, seq_time;

	if (0 == myrank) {
		Vector v(vectorSize*mysize);
		randomize(v);

		double seq_start, seq_finish;
		seq_start = MPI_Wtime();
		global = twoNorm(v);
		seq_finish = MPI_Wtime();
		seq_time = seq_finish - seq_start;

		MPI::COMM_WORLD.Scatter(&v(myrank*vectorSize), vectorSize, MPI::DOUBLE, &x(0), vectorSize, MPI::DOUBLE, 0);
	}
	else {
		MPI::COMM_WORLD.Scatter(NULL, vectorSize, MPI::DOUBLE, &x(0), vectorSize, MPI::DOUBLE, 0);
	}

	double local_start, local_finish, local_elapsed, para_time;
	MPI::COMM_WORLD.Barrier();
	local_start = MPI_Wtime();
	mpi = mpiTwoNorm(x);
	local_finish = MPI_Wtime();
	local_elapsed = local_finish - local_start;

	MPI::COMM_WORLD.Reduce(&local_elapsed, &para_time, 1, MPI::DOUBLE, MPI::MAX, 0);

	if (0 == myrank) {
		std::cout << "#\tgvSize\tseq_time\tpara_time\tspeedup\tdiff" << std::endl;
		std::cout << mysize << "\t" << vectorSize*mysize << "\t" << seq_time << "\t" << para_time << "\t" << seq_time/para_time << "\t" << std::abs(global-mpi) << std::endl;
	}

  MPI::Finalize();

  return 0;
}

double mpiTwoNorm(const Vector& v) {
	double sum = 0.0, mpi;
  for (int i = 0; i < v.numRows(); ++i)
    sum += v(i)*v(i);
	MPI::COMM_WORLD.Allreduce(&sum, &mpi, 1, MPI::DOUBLE, MPI::SUM);
	return std::sqrt(mpi);
}
