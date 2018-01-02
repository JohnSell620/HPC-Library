/**
		AOS.hpp
		Description: Implements the COO Matrix as array of structs.

		@author Johnny Sellers
		@version 0.1 05/30/2017
*/
#ifndef __AOS_HPP
#define __AOS_HPP

#include <cassert>
#include <tuple>
#include <vector>
#include <iostream>
#include <fstream>

#include "Vector.hpp"

class AOSMatrix {
private:
  typedef std::vector<double>::size_type size_type;
  typedef std::tuple<size_type, size_type, double> element;

public:
  AOSMatrix(size_type M, size_type N) : iRows(M), jCols(N) {}

  void push_back(size_type i, size_type j, double val) {
    assert(i < iRows && i >= 0);
    assert(j < jCols && j >= 0);

    arrayData.push_back(element(i, j, val));
  }

  void clear() {
    arrayData.clear();
  }

  void reserve(size_type n) {
    assert(n >= 0);

    arrayData.reserve(n);
  }

  size_type numRows()     const { return iRows; }
  size_type numCols()     const { return jCols; }
  size_type numNonzeros() const { return arrayData.size(); }

  void matvec(const Vector& x, Vector& y) const {
    for (size_type k = 0; k < arrayData.size(); ++k) {
      y(std::get<1>(arrayData[k])) += std::get<2>(arrayData[k]) * x(std::get<0>(arrayData[k]));
    }
  }

  void streamMatrix(std::ostream& outputFile) const {

		// Write header
    outputFile << "AOSMATRIX" << std::endl;
    outputFile << iRows << " " << jCols << std::endl;

    // Write data
    for (size_type i = 0; i < arrayData.size(); ++i) {
      outputFile << std::get<0>(arrayData[i]) << " ";
      outputFile << std::get<1>(arrayData[i]) << " ";
      outputFile << std::get<2>(arrayData[i]) << " ";
      outputFile << std::endl;
    }

    // Write tailer
    outputFile << "END" << std::endl;
  }


private:
  size_type iRows, jCols;
  std::vector<element> arrayData;
};

// Function prototypes
Vector operator*(const AOSMatrix& A, const Vector& x);
void matvec(const AOSMatrix& A, const Vector& x, Vector& y);
void piscetize(AOSMatrix& A, size_t xpoints, size_t ypoints);
void writeMatrix(const AOSMatrix& A, const std::string& filename);
void streamMatrix(const AOSMatrix&A);
void streamMatrix(const AOSMatrix&A, std::ostream& outputFile);

#endif // __AOS_HPP
