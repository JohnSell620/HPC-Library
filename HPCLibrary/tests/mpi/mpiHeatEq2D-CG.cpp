/*
 * mpiHeatEq2D.cpp
 * Description: Using CG method to solve 2D steady-state heat equation.
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#include <iostream>
#include <mpi.h>
#include "Grid.hpp"
#include "Laplacian.hpp"

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

  size_t xsize              = 128, ysize = 128;
  size_t max_iters          = xsize;
  double tol                = 1.E-4;

  if (argc >= 2) xsize = ysize = std::stol(argv[1]);
  if (argc >= 3) max_iters     = std::stol(argv[2]);
  if (argc >= 4) tol           = std::stol(argv[3]);

  Grid X0(xsize, ysize), X1(xsize, ysize);

  for (size_t i = 0; i < ysize+2; ++i) {
    X1(0, i) = X0(0, i) = 1.0;
  }

  Laplacian A;
  cg(A, X1, X0, max_iters, tol);

	MPI_Finalize();
  return 0;
}
