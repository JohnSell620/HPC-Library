/**
		main.cpp
		Description: Example using static library libhpc.h.

		@author Johnny Sellers
		@version 0.1 07/03/2018
*/
#include "libhpc.h"

int main(int argc, char *argv[]) {

  Matrix A(20,20);
  randomizeMatrix(A);
  writeMatrix(A, std::cout);

  Matrix R(20,20);
  qr(A,R);
  writeMatrix(R, std::cout);

  CSCMatrix B(16,16);
  Vector x(16), y(16);
  B.piscretize(4,4);
  randomize(x);
  randomize(y);

  B.matvec(x, y);
  B.streamMatrix(std::cout);

  return 0;
}
