/*
 * CSR.cpp
 * Description: Implements CSR.hpp function prototypes
 *
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
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

void CSRMatrix::piscretize(size_type xpoints, size_type ypoints) {
  assert(this->numRows() == this->numCols());
  assert(xpoints*ypoints == this->numRows());

  this->clear();
	tuple_list elements;

  for (size_type j = 0; j < xpoints; j++) {
    for (size_type k = 0; k < ypoints; k++) {
      size_type jrow = j*ypoints + k;

      if (j != 0) {
				size_type jcol = (j-1)*ypoints + k;
				std::tuple<size_type, size_type, double> tuple = std::make_tuple(jrow, jcol, -1.0);
				elements.push_back(tuple);
      }
      if (k != 0) {
				size_type jcol = j*ypoints + (k-1);
				std::tuple<size_type, size_type, double> tuple = std::make_tuple(jrow, jcol, -1.0);
				elements.push_back(tuple);
      }

			std::tuple<size_type, size_type, double> tuple = std::make_tuple(jrow, jrow, 4);
			elements.push_back(tuple);

      if (k != ypoints-1) {
				size_type jcol = j*ypoints + (k+1);
				std::tuple<size_type, size_type, double> tuple = std::make_tuple(jrow, jcol, -1.0);
				elements.push_back(tuple);
      }
      if (j != xpoints-1) {
				size_type jcol = (j+1)*ypoints + k;
				std::tuple<size_type, size_type, double> tuple = std::make_tuple(jrow, jcol, -1.0);
				elements.push_back(tuple);
      }
    }
  }
	elements.sort();
	this->openForPushBack();

	size_type colIndice, rowIndice; double value;
	for (tuple_list::iterator it = elements.begin(); it != elements.end(); ++it) {
		std::tuple<size_type, size_type, double> temp = *it;
		rowIndice = std::get<0>(temp);
		colIndice = std::get<1>(temp);
		value = std::get<2>(temp);
		if (value != 0)
			this->push_back(rowIndice, colIndice, value);
	}
	this->closeForPushBack();
}
