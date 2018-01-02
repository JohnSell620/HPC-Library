/**
		CSR.cpp
		Description: Implements CSR.hpp function prototypes

		@author Johnny Sellers
		@version 0.1 05/10/2017
*/
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>

#include "CSR.hpp"

#ifdef _OPENMP

void ompMatvec(const CSRMatrix& A, const Vector& x, Vector& y) {
  A.ompMatvec(x, y);
}

#endif

Vector operator*(const CSRMatrix& A, const Vector& x) {
  assert(A.numCols() == x.numRows());

  Vector y(A.numRows(), 0.0);
  matvec(A, x, y);

  return y;
}

void matvec(const CSRMatrix& A, const Vector& x, Vector& y) {
  A.matvec(x, y);
}

void writeMatrix(const CSRMatrix& A, const std::string& filename) {
  std::ofstream outputFile(filename);
  streamMatrix(A, outputFile);
  outputFile.close();
}

void streamMatrix(const CSRMatrix&A) {
  A.streamMatrix(std::cout);
}

void streamMatrix(const CSRMatrix&A, std::ostream& outputFile) {
  A.streamMatrix(outputFile);
}
