/**
		mpi2norm_driver.cpp, AMATH 583
		Description: mpi2norm_driver.cpp for PS8

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
	double global, mpi;

	if (0 == myrank) {
		Vector v(vectorSize*mysize);
		randomize(v);
		global = twoNorm(v);

		MPI::COMM_WORLD.Scatter(&v(myrank*vectorSize), vectorSize, MPI::DOUBLE, &x(0), vectorSize, MPI::DOUBLE, 0);
	}
	else {
		MPI::COMM_WORLD.Scatter(NULL, vectorSize, MPI::DOUBLE, &x(0), vectorSize, MPI::DOUBLE, 0);
	}

	mpi = mpiTwoNorm(x);

	if (0 == myrank) {
		std::cout << "#\tglobal\tmpi\tdiff" << std::endl;
		std::cout << mysize << "\t" << global << "\t" << mpi << "\t" << std::abs(global-mpi) << std::endl;
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
