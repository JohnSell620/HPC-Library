/**
		CSCMatrix.hpp
		Description: Compressed-sparse-column implementation.

		@author Johnny Sellers
		@version 0.1 04/28/2017
*/
#ifndef CSCMATRIX_HPP
#define CSCMATRIX_HPP

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
    os << "AMATH 583 COOMATRIX" << std::endl;
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
    os << "THIS IS THE END" << std::endl;
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
		ofs << "AMATH 583 COOMatrix" << std::endl;
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
		ofs << "THIS IS THE END" << std::endl;

		ofs.close();
	}

private:
	int iRows, jCols;
	bool is_open;
	std::vector<int> rowIndices, colIndices;
	std::vector<double> arrayData;
};


void CSCMatrix::readMatrix(std::string fileName) {
	std::string line;
	std::ifstream inputFile (fileName);

	int iRows, jCols, nonzeros;
	std::string string_input;

	std::getline(inputFile, string_input);
	if (string_input.compare("AMATH 583 COOMATRIX") != 0) {
		std::cout << "Error: incorrect header. Correct header: AMATH 583 COOMATRIX" << std::endl;
		std::exit(-2);
	}

	std::getline(inputFile, string_input);
	jCols = std::stoi(string_input);
	if (jCols <= 0) {
		std::cout << "Error: columns <= 0?" << std::endl;
		std::exit(-2);
	}

	std::getline(inputFile, string_input);
	iRows = std::stoi(string_input);
	if (iRows <= 0) {
		std::cout << "Error: rows <= 0?" << std::endl;
		std::exit(-2);
	}

	std::getline(inputFile, string_input);
	nonzeros = std::stoi(string_input);
	if (nonzeros < 0) {
		std::cout << "Error: nonzeros < 0?" << std::endl;
		std::exit(-2);
	}

	this->setDims(iRows, jCols);

	int rowIndice, colIndice;
	double number, value;
	tuple_list elements;
	std::vector<double> element;

	for (int i = 0; i < nonzeros; ++i) {
		std::getline(inputFile, string_input);
		std::stringstream iss(string_input);

		while(iss >> number) {
			element.push_back(number);
		}
		rowIndice = element[0]; colIndice = element[1];
		value = element[2];

		std::tuple<int, int, double> tuple = std::make_tuple(colIndice, rowIndice, value);

		elements.push_back(tuple);

		element.clear();
	}

	// Sort the matrix elements by column indices
	elements.sort();
	this->openForPushBack();

	for (tuple_list::iterator it = elements.begin(); it != elements.end(); ++it) {
		std::tuple<int, int, double> temp = *it;
		colIndice = std::get<0>(temp);
		rowIndice = std::get<1>(temp);
		value = std::get<2>(temp);
		if (value != 0)
			this->push_back(rowIndice, colIndice, value);
	}

	this->closeForPushBack();

	assert(nonzeros == this->numNonzeros());

	std::getline(inputFile, string_input);
	if (string_input.compare("THIS IS THE END") != 0) {
		std::cout << "Error: incorrect trailer. Correct trailer: THIS IS THE END" << std::endl;
		std::exit(-2);
	}

	inputFile.close();
}

void CSCMatrix::piscetize(int xpoints, int ypoints) {
  assert(this->numRows() == this->numCols());
  assert(xpoints*ypoints == this->numRows());

  this->clear();
	tuple_list elements;

  for (int j = 0; j < xpoints; j++) {
    for (int k = 0; k < ypoints; k++) {
      int jrow = j*ypoints + k;

      if (j != 0) {
				int jcol = (j-1)*ypoints + k;
				std::tuple<int, int, double> tuple = std::make_tuple(jcol, jrow, -1.0);
				elements.push_back(tuple);
      }
      if (k != 0) {
				int jcol = j*ypoints + (k-1);
				std::tuple<int, int, double> tuple = std::make_tuple(jcol, jrow, -1.0);
				elements.push_back(tuple);
      }

			std::tuple<int, int, double> tuple = std::make_tuple(jrow, jrow, 4);
			elements.push_back(tuple);

      if (k != ypoints-1) {
				int jcol = j*ypoints + (k+1);
				std::tuple<int, int, double> tuple = std::make_tuple(jcol, jrow, -1.0);
				elements.push_back(tuple);
      }
      if (j != xpoints-1) {
				int jcol = (j+1)*ypoints + k;
				std::tuple<int, int, double> tuple = std::make_tuple(jcol, jrow, -1.0);
				elements.push_back(tuple);
      }
    }
  }
	elements.sort();
	this->openForPushBack();

	int colIndice, rowIndice; double value;
	for (tuple_list::iterator it = elements.begin(); it != elements.end(); ++it) {
		std::tuple<int, int, double> temp = *it;
		colIndice = std::get<0>(temp);
		rowIndice = std::get<1>(temp);
		value = std::get<2>(temp);
		if (value != 0)
			this->push_back(rowIndice, colIndice, value);
	}
	this->closeForPushBack();
}

#endif // CSCMATRIX_HPP
