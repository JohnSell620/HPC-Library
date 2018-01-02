/**
		CSR.hpp
		Description: Compressed-sparse-row matrix implementation.

		@author Johnny Sellers
		@version 0.1 05/10/2017
*/
#ifndef CSR_HPP
#define CSR_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <tuple>
#include <algorithm>
#include "Vector.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::list<std::tuple<int, int, double> > tuple_list;
typedef unsigned long size_type;

// Prototype for wrapper function to read environment variables
std::string getenv_to_string(const char *in);
std::string getenv(const std::string& in);

class CSRMatrix {
public:
	CSRMatrix(size_type M, size_type N) : is_open(false), iRows(M), jCols(N) {}

	CSRMatrix() : iRows(0), jCols(0) {}

  ~CSRMatrix();

	void piscetize(size_type xpoints, size_type ypoints);

	void openForPushBack() { is_open = true; }

	void closeForPushBack() { is_open = false;
		// Compress rowIndices

		std::vector<size_type> ptr(1, 0); ptr[0] = 0;
		size_type cum = 1, i = 1, j = 0;
		do {
			if (rowIndices[i-1] == rowIndices[i]) {
				while (rowIndices[i-1] == rowIndices[i]) {
					++cum;
					++i;
				}
				ptr.push_back(cum + ptr[j]);
				++j;
			}
			else if (rowIndices[i-1] != rowIndices[i] && rowIndices[i] != rowIndices[i+1]) {
				cum = 1;
				ptr.push_back(cum + ptr[j]);
				++i;
				++j;
			}
			else {
				cum = 1;
				++i;
			}
		} while (i < rowIndices.size());
		if (rowIndices[rowIndices.size()-1] != rowIndices[rowIndices.size()-2]) {
			ptr.push_back(1 + ptr[ptr.size()-1]);
		}

		rowIndices.clear();
		rowIndices = ptr;
	}

	void push_back(size_type i, size_type j, double val) { assert(is_open);
		rowIndices.push_back(i);
    colIndices.push_back(j);
    arrayData.push_back(val);
  }

	void clear() {
    rowIndices.clear();
    colIndices.clear();
    arrayData.clear();
  }

	void reserve(size_type n) {
    assert(n >= 0);

    colIndices.reserve(n);
    arrayData.reserve(n);
  }

	size_type numRows()     const { return iRows; }
	size_type numCols()     const { return jCols; }
	size_type numNonzeros() const { return arrayData.size(); }

	void setDims(size_type i, size_type j) { iRows = i; jCols = j; }
	size_type ptrSize()			const { return colIndices.size(); }
	size_type ptrVal(size_type i)	const { return colIndices[i]; }

	void matvec(const Vector& x, Vector& y) const {
		for (size_type i = 0; i < iRows; ++i) {
			for (size_type j = rowIndices[i]; j < rowIndices[i+1]; ++j) {
				y(i) += arrayData[j] * x(colIndices[j]);
			}
		}
	}

private:
	size_type iRows, jCols;
	bool is_open;
	std::vector<size_type> rowIndices, colIndices;
	std::vector<double> arrayData;
};

void CSRMatrix::piscetize(size_type xpoints, size_type ypoints) {
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

// Function prototypes
void ompMatvec(const CSRMatrix& A, const Vector& x, Vector& y);
Vector operator*(const CSRMatrix& A, const Vector& x);
void matvec(const CSRMatrix& A, const Vector& x, Vector& y);
void writeMatrix(const CSRMatrix& A, const std::string& filename);
void streamMatrix(const CSRMatrix&A);
void streamMatrix(const CSRMatrix&A, std::ostream& outputFile);

#endif // CSR_HPP
