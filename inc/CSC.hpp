/**
		CSC.hpp
		Description: Compressed-sparse-column implementation.

		@author Johnny Sellers
		@version 0.1 04/28/2017
*/
#ifndef CSC_HPP
#define CSC_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <list>
#include <tuple>
#include <algorithm>
#include "Vector.hpp"

typedef std::list<std::tuple<int, int, double> > tuple_list;

class CSCMatrix {
public:
	CSCMatrix(int M, int N) : iRows(M), jCols(N) {}
	CSCMatrix() : iRows(0), jCols(0) {}

	void readMatrix(std::string fileName);

	void piscetize(int xpoints, int ypoints);

	void openForPushBack() { is_open = true; }

	void closeForPushBack() { is_open = false;

		// Compress colIndices
		std::vector<int> ptr(1, 0); ptr[0] = 0;
		int cum = 1, i = 1, j = 0;
		do {
			if (colIndices[i-1] == colIndices[i]) {
				while (colIndices[i-1] == colIndices[i]) {
					++cum;
					++i;
				}
				ptr.push_back(cum + ptr[j]);
				++j;
			}
			else if (colIndices[i-1] != colIndices[i] && colIndices[i] != colIndices[i+1]) {
				cum = 1;
				ptr.push_back(cum + ptr[j]);
				++i;
				++j;
			}
			else {
				cum = 1;
				++i;
			}
		} while (i < colIndices.size());
		if (colIndices[colIndices.size()-1] != colIndices[colIndices.size()-2]) {
			ptr.push_back(1 + ptr[ptr.size()-1]);
		}

		colIndices.clear();
		colIndices = ptr;
	}

	void push_back(int i, int j, double val) { assert(is_open);
		rowIndices.push_back(i);
    colIndices.push_back(j);
    arrayData.push_back(val);
  }

	void clear() {
    rowIndices.clear();
    colIndices.clear();
    arrayData.clear();
  }

	void reserve(int n) {
    assert(n >= 0);

    rowIndices.reserve(n);
    arrayData.reserve(n);
  }

	int numRows()     const { return iRows; }
	int numCols()     const { return jCols; }
	int numNonzeros() const { return arrayData.size(); }

	void setDims(int i, int j) { iRows = i; jCols = j; }
	int ptrSize()			const { return colIndices.size(); }
	int ptrVal(int i)	const { return colIndices[i]; }

	void matvec(const Vector& x, Vector& y) const {
		for (int i = 0; i < jCols; ++i) {
			for (int j = colIndices[i]; j < colIndices[i+1]; ++j) {
				y(rowIndices[j]) += arrayData[j] * x(rowIndices[j]);
			}
		}
	}

	void streamMatrix(std::ostream& os) const {
    assert(arrayData.size() == rowIndices.size());

		// Decompress colIndices
		std::vector<int> decomp(arrayData.size());

		if (!is_open) {
			int i = 0, j = 0, k = 0;
			while (i < colIndices.size()) {
				while (k < colIndices[i+1] - colIndices[i]) {
					decomp[j] = i;
					++j;
					++k;
				}
				++i;
				k = 0;
			}
		}
		else { decomp = colIndices; }

		// Write header
    os << "CSCMATRIX" << std::endl;
    os << jCols << std::endl;
		os << iRows << std::endl;
		os << arrayData.size() << std::endl;

		// Write data
    for (int i = 0; i < arrayData.size(); ++i) {
      os << rowIndices[i] << " ";
      os << decomp[i] << " ";
      os << arrayData[i] << " ";
      os << std::endl;
    }
		// Write trailer
    os << "END" << std::endl;
  }

	void streamMatrix(std::string file) const {
		assert(arrayData.size() == rowIndices.size());

		// Decompress colIndices
		std::vector<int> decomp(arrayData.size());

		if (!is_open) {
			int i = 0, j = 0, k = 0;
			while (i < colIndices.size()) {
				while (k < colIndices[i+1] - colIndices[i]) {
					decomp[j] = i;
					++j;
					++k;
				}
				++i;
				k = 0;
			}
		}
		else { decomp = colIndices; }

		std::ofstream ofs;
		ofs.open (file, std::ios::out);

		// Write header
		ofs << "CSCMatrix" << std::endl;
		ofs << jCols << std::endl;
		ofs << iRows << std::endl;
		ofs << arrayData.size() << std::endl;

		// Write data
		for (int i = 0; i < arrayData.size(); ++i) {
			ofs << rowIndices[i] << " ";
			ofs << decomp[i] << " ";
			ofs << arrayData[i] << std::endl;
		}
		// Write trailer
		ofs << "END" << std::endl;

		ofs.close();
	}

private:
	int iRows, jCols;
	bool is_open;
	std::vector<int> rowIndices, colIndices;
	std::vector<double> arrayData;
};

#endif // CSC_HPP
