/*
 * mpi2norm_driver.cpp, AMATH 583
 * Description: mpi2norm_driver.cpp for PS8
 * @author Johnny Sellers
 * @version 0.1 05/30/2017
 */
#include <iostream>
#include <string>
#include <cmath>
#include <mpi.h>
#include "Vector.hpp"
// #include "../../inc/Vector.hpp"

double mpiTwoNorm(const Vector& v);

int main(int argc, char* argv[]) {
  MPI_Init(NULL, NULL);

  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (argc != 2) {
		if (0 == world_rank) std::cout << "Usage: " << argv[0] << " [ vectorSize ]" << std::endl;
		return -1;
	}

	size_t vectorSize = std::stol(argv[1]);
	Vector x(vectorSize);
	double global, mpi;

	if (0 == world_rank) {
		Vector v(vectorSize*world_size);
		randomize(v);
		global = twoNorm(v);

		MPI_Scatter(&v(world_rank*vectorSize), vectorSize, MPI_DOUBLE, &x(0), vectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Scatter(NULL, vectorSize, MPI_DOUBLE, &x(0), vectorSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	mpi = mpiTwoNorm(x);

	if (0 == world_rank) {
		std::cout << "#\tglobal\tmpi\tdiff" << std::endl;
		std::cout << world_size << "\t" << global << "\t" << mpi << "\t" << std::abs(global-mpi) << std::endl;
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
