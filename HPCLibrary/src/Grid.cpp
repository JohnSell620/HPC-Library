/*
 * Grid.hpp
 * Description: Grid class.
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#include <cassert>
#include <cmath>
#include <iostream>
#include "Grid.hpp"

Grid operator-(const Grid& X, const Grid& Y) {
  Grid Z(X);
  for (size_t i = 0; i < Z.numX(); ++i) {
    for (size_t j = 0; j < Z.numY(); ++j) {
      Z(i, j) = X(i, j) - Y(i, j);
    }
  }
  return Z;
}

Grid operator+(const Grid& X, const Grid& Y) {
  Grid Z(X);
  for (size_t i = 0; i < Z.numX(); ++i) {
    for (size_t j = 0; j < Z.numY(); ++j) {
      Z(i, j) = X(i, j) + Y(i, j);
    }
  }
  return Z;
}

void operator+=(Grid& Z, const Grid& X) {
  for (size_t i = 0; i < Z.numX(); ++i) {
    for (size_t j = 0; j < Z.numY(); ++j) {
      Z(i, j) += X(i, j);
    }
  }
}

void operator-=(Grid& Z, const Grid& X) {
  for (size_t i = 0; i < Z.numX(); ++i) {
    for (size_t j = 0; j < Z.numY(); ++j) {
      Z(i, j) -= X(i, j);
    }
  }
}

Grid operator*(double a, const Grid& Y) {
  Grid Z(Y);
  for (size_t i = 1; i < Z.numX() - 1; ++i) {
    for (size_t j = 1; j < Z.numY() - 1; ++j) {
      Z(i, j) = a * Y(i, j);
    }
  }
  return Z;
}

double dot(const Grid& X, const Grid& Y) {
  double sum = 0.0;
  for (size_t i = 0; i < X.numX(); ++i) {
    for (size_t j = 0; j < X.numY(); ++j) {
      sum += X(i, j) * Y(i, j);
    }
  }
  return sum;
}

double jacobiStep(const Grid& x, Grid& y) {
  assert(x.numX() == y.numX() && x.numY() == y.numY());
  double rnorm = 0.0;

  for (size_t i = 1; i < x.numX() - 1; ++i) {
    for (size_t j = 1; j < x.numY() - 1; ++j) {
      y(i, j) = (x(i - 1, j) + x(i + 1, j) + x(i, j - 1) + x(i, j + 1)) / 4.0;
      rnorm += (y(i, j) - x(i, j)) * (y(i, j) - x(i, j));
    }
  }

  return std::sqrt(rnorm);
}

void swap(Grid& x, Grid& y) { x.swap(y); }

size_t jacobi(Grid& X0, Grid& X1, size_t max_iters, double tol) {
  for (size_t iter = 0; iter < max_iters; ++iter) {
    double rnorm = jacobiStep(X0, X1);
    std::cout << "||X0-X1|| = " << rnorm << std::endl;
    if (rnorm < tol) return 0;
    swap(X0, X1);
  }
  return -1;
}
