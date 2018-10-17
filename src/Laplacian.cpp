/*
 * Laplacian.cpp
 * Description: 5-point Laplacian stencil class.
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#include <mpi.h>
#include <iostream>
#include <cmath>
#include "Laplacian.hpp"

Grid operator*(const Laplacian& A, const Grid& x) {
  Grid y(x);
  for (size_t i = 1; i < x.numX() - 1; ++i) {
   for (size_t j = 1; j < x.numY() - 1; ++j) {
     y(i, j) = x(i, j) - (x(i - 1, j) + x(i + 1, j) + x(i, j - 1) + x(i, j + 1)) / 4.0;
   }
  }
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Status status1, status2;

	int src  = (world_rank == 0 ? MPI_PROC_NULL : world_rank-1);
	int dest = (world_rank == world_size-1 ? MPI_PROC_NULL : world_rank+1);
	int M = x.numX()-1, N = x.numY()+2;

	MPI_Sendrecv(&y(M,0), N, MPI_DOUBLE, dest, 23, &y(1,0), N, MPI_DOUBLE,
               src, 23, MPI_COMM_WORLD, &status1);
	MPI_Sendrecv(&y(1,0), N, MPI_DOUBLE, src, 34, &y(M,0), N, MPI_DOUBLE,
               dest, 34, MPI_COMM_WORLD, &status2);

 return y;
}

double mpiDot(const Grid& X, const Grid& Y) {
  double sum = 0.0, rho;
  for (size_t i = 0; i < X.numX(); ++i) {
    for (size_t j = 0; j < X.numY(); ++j) {
      sum += X(i, j) * Y(i, j);
    }
  }
	MPI_Allreduce(&sum, &rho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return rho;
}

size_t ir(const Laplacian& A, Grid& x, const Grid& b, size_t max_iter, double tol) {
  for (size_t iter = 0; iter < max_iter; ++iter) {
    Grid   r   = b - A * x;
    double rho = mpiDot(r, r);
		std::cout << "||r|| = " << std::sqrt(rho) << std::endl;
    if (std::sqrt(rho) < tol) return iter;
    x += r;
  }
  return max_iter;
}

size_t cg(const Laplacian& A, Grid& x, const Grid& b, size_t max_iter, double tol) {
	Grid   r   = b - A * x, p(b);
  double rho = mpiDot(r, r), rho_1 = 0.0;
  for (size_t iter = 0; iter < max_iter; ++iter) {
    std::cout << "||r|| = " << std::sqrt(rho) << std::endl;

    if (iter == 0) {
      p = r;
    } else {
      double beta = (rho / rho_1);
      p           = r + beta * p;
    }

    Grid   q     = A * p;
    double alpha = rho / mpiDot(p, q);

    x += alpha * p;

    rho_1 = rho;
    r -= alpha * q;
    rho = mpiDot(r, r);

    if (rho < tol) return iter;
  }
  return max_iter;
}
