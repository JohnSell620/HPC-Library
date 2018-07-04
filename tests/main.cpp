/**
		main.cpp
		Description: Example using static library libhpc.h.

		@author Johnny Sellers
		@version 0.1 07/03/2018
*/
#include <iostream>
#include <string>
#include <fstream>
#include "Matrix.hpp"

int main(int argc, char *argv[]) {

  Matrix A(20,20);
  randomizeMatrix(A);
  writeMatrix(A, std::cout);

  Matrix R(20,20);
  qr(A,R);
  writeMatrix(R, std::cout);


  return 0;
}
