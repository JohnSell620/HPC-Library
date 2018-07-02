/**
		Matrix.cpp
		Description: Matrix class.

		@author Johnny Sellers
		@version 0.1 05/10/2017
*/
#include <cassert>
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include "../inc/Matrix.hpp"


Matrix operator*(const Matrix& A, const Matrix& B) {
  assert(A.numCols() == B.numRows());

  Matrix C(A.numRows(), B.numCols());
  zeroizeMatrix(C);

  if (A.numCols() % 4 == 0)
    hoistedBlockedTiledMultiply2x2(A, B, C);
  else if (A.numCols() % 2 == 0)
    hoistedCopyBlockedTiledMultiply2x2(A, B, C);
  else
    multiply(A, B, C);

  return C;
}


Matrix operator+(const Matrix& A, const Matrix& B) {
  assert(A.numCols() == B.numCols() && A.numRows() == B.numRows());

  Matrix C(A.numRows(), A.numCols());

  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < A.numCols(); ++j) {
      C(i,j) = A(i,j) + B(i,j);
    }
  }
  return C;
}


Matrix operator-(const Matrix& A, const Matrix& B) {
  assert(A.numCols() == B.numCols() && A.numRows() == B.numRows());

  Matrix C(A.numRows(), A.numCols());

  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < A.numCols(); ++j) {
      C(i,j) = A(i,j) - B(i,j);
    }
  }
  return C;
}


void multiply(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < B.numCols(); ++j) {
      for (int k = 0; k < A.numCols(); ++k) {
				C(i,j) += A(i,k) * B(k,j);
      }
    }
  }
}


void hoistedMultiply(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < B.numCols(); ++j) {
      double t = C(i,j);
      for (int k = 0; k < A.numCols(); ++k) {
				t += A(i,k) * B(k,j);
      }
      C(i,j) = t;
    }
  }
}


void tiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); i += 2) {
    for (int j = 0; j < B.numCols(); j += 2) {
      for (int k = 0; k < A.numCols(); ++k) {
				C(i  , j  ) += A(i  , k) * B(k, j  );
				C(i  , j+1) += A(i  , k) * B(k, j+1);
				C(i+1, j  ) += A(i+1, k) * B(k, j  );
				C(i+1, j+1) += A(i+1, k) * B(k, j+1);
      }
    }
  }
}


void tiledMultiply2x4(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); i += 2) {
    for (int j = 0; j < B.numCols(); j += 4) {
      for (int k = 0; k < A.numCols(); ++k) {
				C(i  , j  ) += A(i  , k) * B(k, j  );
				C(i  , j+1) += A(i  , k) * B(k, j+1);
				C(i  , j+2) += A(i  , k) * B(k, j+2);
				C(i  , j+3) += A(i  , k) * B(k, j+3);
				C(i+1, j  ) += A(i+1, k) * B(k, j  );
				C(i+1, j+1) += A(i+1, k) * B(k, j+1);
				C(i+1, j+2) += A(i+1, k) * B(k, j+2);
				C(i+1, j+3) += A(i+1, k) * B(k, j+3);
      }
    }
  }
}


void tiledMultiply4x2(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); i += 4) {
    for (int j = 0; j < B.numCols(); j += 2) {
      for (int k = 0; k < A.numCols(); ++k) {
				C(i  , j  ) += A(i  , k) * B(k, j  );
				C(i  , j+1) += A(i  , k) * B(k, j+1);
				C(i+1, j  ) += A(i+1, k) * B(k, j  );
				C(i+1, j+1) += A(i+1, k) * B(k, j+1);
				C(i+2, j  ) += A(i+2, k) * B(k, j  );
				C(i+2, j+1) += A(i+2, k) * B(k, j+1);
				C(i+3, j  ) += A(i+3, k) * B(k, j  );
				C(i+3, j+1) += A(i+3, k) * B(k, j+1);
      }
    }
  }
}


void tiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); i += 4) {
    for (int j = 0; j < B.numCols(); j += 4) {
      for (int k = 0; k < A.numCols(); ++k) {
				C(i  , j  ) += A(i  , k) * B(k, j  );
				C(i  , j+1) += A(i  , k) * B(k, j+1);
				C(i  , j+2) += A(i  , k) * B(k, j+2);
				C(i  , j+3) += A(i  , k) * B(k, j+3);
				C(i+1, j  ) += A(i+1, k) * B(k, j  );
				C(i+1, j+1) += A(i+1, k) * B(k, j+1);
				C(i+1, j+2) += A(i+1, k) * B(k, j+2);
				C(i+1, j+3) += A(i+1, k) * B(k, j+3);
				C(i+2, j  ) += A(i+2, k) * B(k, j  );
				C(i+2, j+1) += A(i+2, k) * B(k, j+1);
				C(i+2, j+2) += A(i+2, k) * B(k, j+2);
				C(i+2, j+3) += A(i+2, k) * B(k, j+3);
				C(i+3, j  ) += A(i+3, k) * B(k, j  );
				C(i+3, j+1) += A(i+3, k) * B(k, j+1);
				C(i+3, j+2) += A(i+3, k) * B(k, j+2);
				C(i+3, j+3) += A(i+3, k) * B(k, j+3);
      }
    }
  }
}


void hoistedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  for (int i = 0; i < A.numRows(); i += 2) {
    for (int j = 0; j < B.numCols(); j += 2) {
      double t00 = C(i,j);       double t01 = C(i,j+1);
      double t10 = C(i+1,j);     double t11 = C(i+1,j+1);

      for (int k = 0; k < A.numCols(); ++k) {
				t00 += A(i  , k) * B(k, j  );
				t01 += A(i  , k) * B(k, j+1);
				t10 += A(i+1, k) * B(k, j  );
				t11 += A(i+1, k) * B(k, j+1);
      }
      C(i,  j) = t00;  C(i,  j+1) = t01;
      C(i+1,j) = t10;  C(i+1,j+1) = t11;
    }
  }
}


void blockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const int blocksize = std::min(A.numRows(), 32);

  for (int ii = 0; ii < A.numRows(); ii += blocksize) {
    for (int jj = 0; jj < B.numCols(); jj += blocksize) {
      for (int kk = 0; kk < A.numCols(); kk += blocksize) {

				for (int i = ii; i < ii+blocksize; i += 2) {
				  for (int j = jj; j < jj+blocksize; j += 2) {
				    for (int k = kk; k < kk+blocksize; ++k) {
				      C(i  , j  ) += A(i  , k) * B(k, j  );
				      C(i  , j+1) += A(i  , k) * B(k, j+1);
				      C(i+1, j  ) += A(i+1, k) * B(k, j  );
				      C(i+1, j+1) += A(i+1, k) * B(k, j+1);
				    }
				  }
				}
      }
    }
  }
}


void hoistedBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const int blocksize = std::min(A.numRows(), 32);

  for (int ii = 0; ii < A.numRows(); ii += blocksize) {
    for (int jj = 0; jj < B.numCols(); jj += blocksize) {
      for (int kk = 0; kk < A.numCols(); kk += blocksize) {

				for (int i = ii; i < ii+blocksize; i += 2) {
				  for (int j = jj; j < jj+blocksize; j += 2) {
				    double t00 = C(i,j);       double t01 = C(i,j+1);
				    double t10 = C(i+1,j);     double t11 = C(i+1,j+1);

				    for (int k = kk; k < kk+blocksize; ++k) {
				      t00 += A(i  , k) * B(k, j  );
				      t01 += A(i  , k) * B(k, j+1);
				      t10 += A(i+1, k) * B(k, j  );
				      t11 += A(i+1, k) * B(k, j+1);
				    }
				    C(i,  j) = t00;  C(i,  j+1) = t01;
				    C(i+1,j) = t10;  C(i+1,j+1) = t11;
				  }
				}
      }
    }
  }
}


void copyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const int blocksize = std::min(A.numRows(), 32);

  for (int ii = 0; ii < A.numRows(); ii += blocksize) {
    for (int jj = 0; jj < B.numCols(); jj += blocksize) {
      for (int kk = 0; kk < A.numCols(); kk += blocksize) {

				Matrix BB(blocksize, blocksize);
				for (int j = jj, jb = 0; j < jj+blocksize; ++j, ++jb) {
				  for (int k = kk, kb = 0; k < kk+blocksize; ++k, ++kb) {
				    BB(jb, kb) = B(k, j);
				  }
				}

				for (int i = ii; i < ii+blocksize; i += 2) {
				  for (int j = jj, jb = 0; j < jj+blocksize; j += 2, jb += 2) {
				    for (int k = kk, kb = 0; k < kk+blocksize; ++k, ++kb) {

				      C(i  , j  ) += A(i  , k) * BB(jb,   kb);
				      C(i  , j+1) += A(i  , k) * BB(jb+1, kb);
				      C(i+1, j  ) += A(i+1, k) * BB(jb,   kb);
				      C(i+1, j+1) += A(i+1, k) * BB(jb+1, kb);
				    }
				  }
				}
      }
    }
  }
}


void hoistedCopyBlockedTiledMultiply2x2(const Matrix& A, const Matrix& B, Matrix& C) {
  const int blocksize = std::min(A.numRows(), 32);

  for (int ii = 0; ii < A.numRows(); ii += blocksize) {
    for (int jj = 0; jj < B.numCols(); jj += blocksize) {
      for (int kk = 0; kk < A.numCols(); kk += blocksize) {

				Matrix BB(blocksize, blocksize);
				for (int j = jj, jb = 0; j < jj+blocksize; ++j, ++jb) {
				  for (int k = kk, kb = 0; k < kk+blocksize; ++k, ++kb) {
				    BB(jb, kb) = B(k, j);
				  }
				}

				for (int i = ii; i < ii+blocksize; i += 2) {
				  for (int j = jj, jb = 0; j < jj+blocksize; j += 2, jb += 2) {

				    double t00 = C(i,j);       double t01 = C(i,j+1);
				    double t10 = C(i+1,j);     double t11 = C(i+1,j+1);
				    for (int k = kk, kb = 0; k < kk+blocksize; ++k, ++kb) {
				      t00 += A(i  , k) * BB(jb  , kb);
				      t01 += A(i  , k) * BB(jb+1, kb);
				      t10 += A(i+1, k) * BB(jb  , kb);
				      t11 += A(i+1, k) * BB(jb+1, kb);
				    }
				    C(i,  j) = t00;  C(i,  j+1) = t01;
				    C(i+1,j) = t10;  C(i+1,j+1) = t11;
				  }
				}
      }
    }
  }
}


void hoistedCopyBlockedTiledMultiply4x4(const Matrix& A, const Matrix& B, Matrix& C) {
  const int blocksize = std::min(A.numRows(), 32);

  for (int ii = 0; ii < A.numRows(); ii += blocksize) {
    for (int jj = 0; jj < B.numCols(); jj += blocksize) {
      for (int kk = 0; kk < A.numCols(); kk += blocksize) {

				Matrix BB(blocksize, blocksize);
				for (int j = jj, jb = 0; j < jj+blocksize; ++j, ++jb) {
				  for (int k = kk, kb = 0; k < kk+blocksize; ++k, ++kb) {
				    BB(jb, kb) = B(k, j);
				  }
				}

				for (int i = ii; i < ii+blocksize; i += 4) {
				  for (int j = jj, jb = 0; j < jj+blocksize; j += 4, jb += 4) {

				    double t00 = C(i,  j);       double t01 = C(i,  j+1);
				    double t10 = C(i+1,j);       double t11 = C(i+1,j+1);
				    double t20 = C(i+2,j);       double t21 = C(i+2,j+1);
				    double t30 = C(i+3,j);       double t31 = C(i+3,j+1);
				    double t02 = C(i,  j+2);     double t03 = C(i,  j+3);
				    double t12 = C(i+1,j+2);     double t13 = C(i+1,j+3);
				    double t22 = C(i+2,j+2);     double t23 = C(i+2,j+3);
				    double t32 = C(i+3,j+2);     double t33 = C(i+3,j+3);
				    for (int k = kk, kb = 0; k < kk+blocksize; ++k, ++kb) {
				      t00 += A(i  , k) * BB(jb  , kb);
				      t10 += A(i+1, k) * BB(jb  , kb);
				      t20 += A(i+2, k) * BB(jb  , kb);
				      t30 += A(i+3, k) * BB(jb  , kb);

				      t01 += A(i  , k) * BB(jb+1, kb);
				      t11 += A(i+1, k) * BB(jb+1, kb);
				      t21 += A(i+2, k) * BB(jb+1, kb);
				      t31 += A(i+3, k) * BB(jb+1, kb);

				      t02 += A(i  , k) * BB(jb+2, kb);
				      t12 += A(i+1, k) * BB(jb+2, kb);
				      t22 += A(i+2, k) * BB(jb+2, kb);
				      t32 += A(i+3, k) * BB(jb+2, kb);

				      t03 += A(i  , k) * BB(jb+3, kb);
				      t13 += A(i+1, k) * BB(jb+3, kb);
				      t23 += A(i+2, k) * BB(jb+3, kb);
				      t33 += A(i+3, k) * BB(jb+3, kb);
				    }
				    C(i,  j  ) = t00;  C(i,  j+1) = t01;
				    C(i+1,j  ) = t10;  C(i+1,j+1) = t11;
				    C(i+2,j  ) = t20;  C(i+2,j+1) = t21;
				    C(i+3,j  ) = t30;  C(i+3,j+1) = t31;
				    C(i,  j+2) = t02;  C(i,  j+3) = t03;
				    C(i+1,j+2) = t12;  C(i+1,j+3) = t13;
				    C(i+2,j+2) = t22;  C(i+2,j+3) = t23;
				    C(i+3,j+2) = t32;  C(i+3,j+3) = t33;
				  }
				}
      }
    }
  }
}


double frobeniusNorm(const Matrix& A) {
  double sum = 0.0;
  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < A.numCols(); ++j) {
      sum += A(i, j)*A(i, j);
    }
  }

  return sqrt(sum);
}


double oneNorm(const Matrix& A) {
  std::vector<double> v(A.numCols(), 0.0);

  for (int j = 0; j < A.numCols(); ++j) {
    for (int i = 0; i < A.numRows(); ++i) {
      v[j] += std::abs(A(i, j));
    }
  }
  return *std::max_element(v.begin(), v.end());
}


double infinityNorm(const Matrix& A) {
  std::vector<double> v(A.numRows(), 0.0);

  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < A.numCols(); ++j) {
      v[i] += std::abs(A(i, j));
    }
  }
  return *std::max_element(v.begin(), v.end());
}


void randomizeMatrix(Matrix &A) {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(2.0, 32.0);
  static auto dice = std::bind(distribution, generator);

  for (int i = 0; i < A.numRows(); ++i) {
    for (int j = 0; j < A.numCols(); ++j) {
      A(i, j) = dice();
    }
  }
}


void zeroizeMatrix(Matrix& C) {
  for (int i = 0; i < C.numRows(); ++i) {
    for (int j = 0; j < C.numCols(); ++j) {
      C(i, j) = 0.0;
    }
  }
}

void matvec(const Matrix& A, const Vector& x, Vector& y) {
	for (int i = 0; i < A.numRows(); ++i)
		for (int j = 0; j < A.numCols(); ++j)
			y(i) += A(i, j) * x(j);
}

Matrix outProd(const Vector& x, const Vector& y) {
  Matrix A (x.numRows(), y.numRows());
  for (int i = 0; i < x.numRows(); ++i)
    for (int j = 0; j < y.numRows(); ++j)
      A(i, j) = x(i) * y(j);

  return A;
}

Matrix qr(const Matrix& A, Matrix& R) {
  assert(A.numRows() == A.numCols());

  // Matrix R;
  // R = A;
  Vector x (A.numRows());
  Vector v (A.numRows());
  Vector e_1 (A.numRows());
  zeroize(e_1);
  double norm_v;

  for (int k = 0; k < A.numCols(); ++k) {
    for (int i = 0; i < A.numRows(); ++i)
      x(i) = A(i,k);

    if (x(0) < 0) e_1(0) = -twoNorm(x);
    else if (x(0) >= 0) e_1(0) = twoNorm(x);
    v = e_1 + x;
    norm_v = twoNorm(v);

    for (int i = 0; i < v.numRows(); ++i)
      v(i) = v(i) / norm_v;

    for (int j = k; j < R.numRows(); ++j) {
      R(k, j) = R(k, j) - 2*v(j) * v(j) * R(k, j);
    }
  }
  return R;
}


Matrix readMatrix(std::string fileName) {
	std::string line;
	std::ifstream inputFile (fileName);

	int iRows, jCols;
	double value;
	std::string string_input;

	std::getline(inputFile, string_input);
	if (string_input.compare("MATRIX") != 0) {
		std::cout << "Error: incorrect header. Correct header: MATRIX" << std::endl;
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

	Matrix A(iRows, jCols);

	for (int i = 0; i < iRows*jCols; ++i) {
		std::getline(inputFile, string_input);
		value = std::stoi(string_input);
		A.setValue(i, value);
	}

	std::getline(inputFile, string_input);
	if (string_input.compare("END") != 0) {
		std::cout << "Error: incorrect trailer. Correct trailer: END" << std::endl;
		std::exit(-2);
	}

	inputFile.close();

	return A;
}

void writeMatrix(const Matrix& A, std::string file) {
  std::ofstream ofs;
  ofs.open (file, std::ios::out);

  ofs << "MATRIX" << std::endl;
  ofs << A.numRows() << std::endl;
  ofs << A.numCols() << std::endl;

  for (int i = 0; i < A.numRows(); ++i)
    for (int j = 0; j < A.numCols(); ++j)
      ofs << A(i, j) << std::endl;

  ofs << "END" << std::endl;
}
