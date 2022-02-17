#if !defined(MATRIX_HPP)
#define MATRIX_HPP

#include <initializer_list>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <climits>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <cmath>
using std::vector;

#include "Vector.hpp"

template<typename T>
class Matrix {
protected:
  int iRows, jCols;
  vector<T> arrayData;

public:
  Matrix(int M, int N) : iRows(M), jCols(N), arrayData(iRows*jCols) {}
  Matrix() { iRows = 0; jCols = 0; }

  ~Matrix() {}

  int numRows() const { return iRows; }
  int numCols() const { return jCols; }
  int size() const { return arrayData.size(); }

  T &operator()(int i, int j);
  const T &operator()(int i, int j) const;

  void setValue(int i, int j, T value);
  void setValue(int k, T value);

  Matrix operator*(const Matrix& B);
  std::vector<T> operator*(const vector<T>& x);
  Matrix operator+(const Matrix& B);
  Matrix operator-(const Matrix& B);
  void operator==(const Matrix& B);
  bool isEqual(const Matrix& B);
  Vector matvec(const Vector& x);

  Matrix partition(int start_row, int end_row, int start_col, int end_col) const;
  static Matrix combine(vector<vector<Matrix>>& matrices);
  static Matrix squareMultiply(const Matrix& A, const Matrix& B);
  static Matrix strassenMultiply(const Matrix& A, const Matrix& B);
  double norm(char type);
  Matrix outerProduct(const vector<T>& x, const vector<T>& y);
  vector<Matrix> qr();      // Modified Gram-Schmidt QR Factorization
  Matrix qrHouseholder();   // Householder QR Factorization

  static Matrix chainMultiply(std::initializer_list<Matrix> matrices);
  static vector<vector<int>> matrixChainOrder(vector<int>& p);
  static Matrix chainMultiply(
    vector<Matrix>& matrices,
    vector<vector<int>>& s,
    int i,
    int j);
  static void printOptimalParens(vector<vector<int>>& s, int i, int j);

  void randomizeMatrix();
  void zeroizeMatrix();
  void overwriteMatrix(std::string fileName);
  static Matrix readMatrix(std::string fileName);
  void writeMatrix(std::string file);
  void writeMatrix(std::ostream& os);
  void prettyPrint(std::string file);
  void prettyPrint(std::ostream& os);

private:
  double twoNorm();
  double oneNorm();
  double infinityNorm();
  double norm(const vector<T>& v);

};

#endif // MATRIX_HPP
