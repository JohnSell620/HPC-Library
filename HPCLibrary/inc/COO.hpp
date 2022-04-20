/*
* COO.hpp
* Description: Implements the COO Matrix as array of structs.
* @author Johnny Sellers
* @version 0.1 05/30/2017
*/
#ifndef COO_HPP
#define COO_HPP

#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>

#include "Vector.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

// Wrapper function to read environment variables (and avoid pointers)
std::string getenv_to_string(const char *in);
std::string getenv(const std::string& in);


class COOMatrix {
private:
    typedef std::vector<double>::size_type size_type;

public:
    COOMatrix(size_type M, size_type N) : iRows(M), jCols(N) {}
    COOMatrix() : iRows(0), jCols(0) {}

    ~COOMatrix() {}

    void push_back(size_type i, size_type j, double val) {
        assert(i < iRows && i >= 0);
        assert(j < jCols && j >= 0);

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

        rowIndices.reserve(n);
        colIndices.reserve(n);
        arrayData.reserve(n);
    }

    size_type numRows()     const { return iRows; }
    size_type numCols()     const { return jCols; }
    size_type numNonzeros() const { return arrayData.size(); }

    void matvec(const Vector& x, Vector& y) const {
        for (size_type k = 0; k < arrayData.size(); ++k) {
            y(rowIndices[k]) += arrayData[k] * x(colIndices[k]);
        }
    }

    void transposeMatvec(const Vector& x, Vector& y) const {
        for (size_type k = 0; k < arrayData.size(); ++k) {
            y(colIndices[k]) += arrayData[k] * x(rowIndices[k]);
        }
    }

    #ifdef _OPENMP

    void ompMatvec(const Vector& x, Vector& y) const {
        std::string env = getenv("OMP_NUM_THREADS");
        int thread_count;
        if (env.compare("") == 0) thread_count = 1;
        else thread_count = stoi(env);

        #pragma omp parallel for num_threads(thread_count)
        for (size_type k = 0; k < arrayData.size(); ++k)
            y(rowIndices[k]) += arrayData[k] * x(colIndices[k]);
    }


    #endif

    void streamMatrix(std::ostream& outputFile) const {
        assert(arrayData.size() == rowIndices.size() && arrayData.size() == colIndices.size());

        // Write header
        outputFile << "COOMATRIX" << std::endl;
        outputFile << iRows << " " << jCols << std::endl;

        // Write data
        for (size_type i = 0; i < arrayData.size(); ++i) {
            outputFile << rowIndices[i] << " ";
            outputFile << colIndices[i] << " ";
            outputFile << arrayData[i] << " ";
            outputFile << std::endl;
        }

        // Write tailer
        outputFile << "END" << std::endl;
    }

private:
    size_type iRows, jCols;
    std::vector<size_type> rowIndices, colIndices;
    std::vector<double> arrayData;
};

// Function prototypes
void ompMatvec(const COOMatrix& A, const Vector& x, Vector& y);
Vector operator*(const COOMatrix& A, const Vector& x);
void matvec(const COOMatrix& A, const Vector& x, Vector& y);
void piscretize(COOMatrix& A, int xpoints, int ypoints);
void writeMatrix(const COOMatrix& A, const std::string& filename);
void streamMatrix(const COOMatrix&A);
void streamMatrix(const COOMatrix&A, std::ostream& outputFile);

#endif // COO_HPP
