/*
 * CSC.cpp
 * Description: Implements CSC.hpp function readMatrix
 * and piscretize methods.
 *
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include "CSC.hpp"


void CSCMatrix::readMatrix(std::string fileName) {
	std::string line;
	std::ifstream inputFile (fileName);

	int iRows, jCols, nonzeros;
	std::string string_input;

	std::getline(inputFile, string_input);
	if (string_input.compare("CSCMATRIX") != 0) {
		std::cout << "Error: incorrect header. Correct header: CSCMATRIX\n";
		std::exit(-2);
	}

	std::getline(inputFile, string_input);
	jCols = std::stoi(string_input);
	if (jCols <= 0) {
		std::cout << "Error: columns <= 0?\n";
		std::exit(-2);
	}

	std::getline(inputFile, string_input);
	iRows = std::stoi(string_input);
	if (iRows <= 0) {
		std::cout << "Error: rows <= 0?\n";
		std::exit(-2);
	}

	std::getline(inputFile, string_input);
	nonzeros = std::stoi(string_input);
	if (nonzeros < 0) {
		std::cout << "Error: nonzeros < 0?\n";
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
	if (string_input.compare("END") != 0) {
		std::cout << "Error: incorrect trailer. Correct trailer: END\n";
		std::exit(-2);
	}

	inputFile.close();
}

void CSCMatrix::piscretize(int xpoints, int ypoints) {
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
