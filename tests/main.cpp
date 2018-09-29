/**
		main.cpp
		Description: Example using static library libhpc.h.

		@author Johnny Sellers
		@version 0.1 07/03/2018
*/
#include <iostream>
#include "libhpc.h"

#ifndef _THREADING
#define _THREADING
#endif

int main(int argc, char *argv[]) {

  Matrix A(20,20);
  randomizeMatrix(A);
  writeMatrix(A, std::cout);

  Matrix R = A;
  qr(A,R);
  writeMatrix(R, std::cout);

  CSCMatrix B(16,16);
  Vector x(16), y(16);
  B.piscretize(4,4);
  randomize(x);
  randomize(y);

  B.matvec(x, y);
  B.streamMatrix(std::cout);

  Vector const& const_x = x;
  size_t partitions = 2;
  double ansp = partitionedTwoNorm(const_x, 2);
  double ansr = recursiveTwoNorm(const_x,2);
  std::cout << "2-Norm of x: " << ansp << std::endl;
  std::cout << "2-Norm of x: " << ansr << std::endl;
  writeVector(x, std::cout);

  return 0;
}
