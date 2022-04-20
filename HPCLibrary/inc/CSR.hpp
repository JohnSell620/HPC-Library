/*
* CSR.hpp
* Description: Compressed-sparse-row matrix implementation.
* @author Johnny Sellers
* @version 0.1 05/10/2017
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

    ~CSRMatrix() {}

    void piscretize(size_type xpoints, size_type ypoints);

    void openForPushBack() { is_open = true; }

    void closeForPushBack() {
        is_open = false;

        // Compress rowIndices

        std::vector<size_type> ptr(1, 0); ptr[0] = 0;
        size_type col = 1, i = 1, j = 0;
        do {
            if (rowIndices[i-1] == rowIndices[i]) {
                while (rowIndices[i-1] == rowIndices[i]) {
                    ++col;
                    ++i;
                }
                ptr.push_back(col + ptr[j]);
                ++j;
            } else if (rowIndices[i-1] != rowIndices[i] && rowIndices[i] != rowIndices[i+1]) {
                col = 1;
                ptr.push_back(col + ptr[j]);
                ++i;
                ++j;
            } else {
                col = 1;
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

    void streamMatrix(std::ostream& os) const {
        assert(arrayData.size() == rowIndices.size());

        // Decompress colIndices
        std::vector<size_type> decomp(arrayData.size());

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
        } else {
            decomp = colIndices;
        }

        // Write header
        os << "CSRMATRIX" << std::endl;
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
    std::vector<size_type> decomp(arrayData.size());

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
    ofs << "CSRMatrix" << std::endl;
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
    size_type iRows, jCols;
    bool is_open;
    std::vector<size_type> rowIndices, colIndices;
    std::vector<double> arrayData;
};

// Function prototypes
void ompMatvec(const CSRMatrix& A, const Vector& x, Vector& y);
Vector operator*(const CSRMatrix& A, const Vector& x);
void matvec(const CSRMatrix& A, const Vector& x, Vector& y);
void writeMatrix(const CSRMatrix& A, const std::string& filename);
void streamMatrix(const CSRMatrix&A);
void streamMatrix(const CSRMatrix&A, std::ostream& outputFile);

#endif // CSR_HPP
