/**
		Matrix.hpp
		Description: Matrix class.

		@author Johnny Sellers
		@version 0.1 05/10/2017
*/
#if !defined(MATRIX_HPP)
#define MATRIX_HPP

#include <vector>
#include <string>
#include "Vector.hpp"


class Matrix {
public:
  Matrix(int M, int N) : iRows(M), jCols(N), arrayData(iRows*jCols) {}
  Matrix() { iRows = 0; jCols = 0; }

  ~Matrix() {}

  double &operator()(int i, int j) { return arrayData[i*jCols + j]; }

  const double &operator()(int i, int j) const {
    return arrayData[i*jCols + j];
  }

  int numRows() const { return iRows; }
  int numCols() const { return jCols; }
  int size() const { return arrayData.size(); }

	void setValue(int k, double value) {
		arrayData[k] = value;
	}

  void matvec(const Vector& x, Vector& y) {
		for (int i = 0; i < iRows; ++i) {
			double t0 = 0; double t1 = 0;
			for (int j = 0; j < jCols; ++j) {
				t0 += arrayData[i*jCols + j] * x(j);
				t1 += arrayData[(i+1)*jCols + j] * x(j);
			}
			y(i) = t0;
			y(i+1) = t1;
		}
	}

  // void operator=(const Matrix& A);

protected:
  int iRows, jCols;
  std::vector<double> arrayData;
};


Matrix operator*(const Matrix& A, const Matrix &B);
Matrix operator+(const Matrix& A, const Matrix &B);
Matrix operator-(const Matrix& A, const Matrix &B);
void multiply(const Matrix& A, const Matrix &B, Matrix& C);
void hoistedMultiply(const Matrix& A, const Matrix &B, Matrix& C);
void tiledMultiply2x2(const Matrix& A, const Matrix&B, Matrix&C);
void hoistedTiledMultiply2x2(const Matrix& A, const Matrix&B, Matrix&C);
void blockedTiledMultiply2x2(const Matrix& A, const Matrix&B, Matrix&C);
void tiledMultiply2x4(const Matrix& A, const Matrix&B, Matrix&C);
void tiledMultiply4x2(const Matrix& A, const Matrix&B, Matrix&C);
void tiledMultiply4x4(const Matrix& A, const Matrix&B, Matrix&C);
void copyBlockedTiledMultiply2x2(const Matrix& A, const Matrix&B, Matrix&C);
void hoistedBlockedTiledMultiply2x2(const Matrix& A, const Matrix&B, Matrix&C);
void hoistedCopyBlockedTiledMultiply2x2(const Matrix& A, const Matrix&B, Matrix&C);
void hoistedCopyBlockedTiledMultiply4x4(const Matrix& A, const Matrix&B, Matrix&C);
double oneNorm(const Matrix& A);
double infinityNorm(const Matrix& A);
double frobeniusNorm(const Matrix& A);
void zeroizeMatrix(Matrix& C);
void randomizeMatrix(Matrix &A);
void matvec(const Matrix& A, const Vector& x, Vector& y);
Matrix outProd(const Vector&x, const Vector& y);
Matrix qr(const Matrix& A, Matrix& R);
Matrix readMatrix(std::string fileName);
void writeMatrix(const Matrix& A, std::string file);
void writeMatrix(const Matrix& A, std::ostream& os);

#endif // MATRIX_HPP
