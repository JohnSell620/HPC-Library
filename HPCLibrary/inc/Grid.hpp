/*
 * Grid.hpp
 * Description: Grid class.
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#ifndef __GRID_HPP
#define __GRID_HPP

#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>


class Grid {

public:
  explicit Grid(size_t x, size_t y) :
    xPoints(x+2), yPoints(y+2), arrayData(xPoints*yPoints) {}
  explicit Grid(size_t x, size_t y, double init) :
    xPoints(x+2), yPoints(y+2), arrayData(xPoints*yPoints, init) {}

  double &operator()(size_t i, size_t j) {
    return arrayData[i*yPoints + j];
  }
  const double &operator()(size_t i, size_t j) const {
    return arrayData[i*yPoints + j];
  }

  size_t numX() const { return xPoints; }
  size_t numY() const { return yPoints; }

  void swap(Grid &x) {
    std::swap(x.xPoints, xPoints);
    std::swap(x.yPoints, yPoints);
    arrayData.swap(x.arrayData);
  }

  void operator=(const Grid& x) {
    assert(x.xPoints == xPoints && x.yPoints == yPoints);
    std::copy(x.arrayData.begin(), x.arrayData.end(), arrayData.begin());
  }

private:
  size_t xPoints, yPoints;
  std::vector<double> arrayData;
};

Grid operator-(const Grid& X, const Grid& Y);
Grid operator+(const Grid& X, const Grid& Y);
void operator+=(Grid& Z, const Grid& X);
void operator-=(Grid& Z, const Grid& X);
Grid operator*(double a, const Grid& Y);
double dot(const Grid& X, const Grid& Y);
double jacobiStep(const Grid& x, Grid& y);
void swap(Grid& x, Grid& y);
size_t jacobi(Grid& X0, Grid& X1, size_t max_iters, double tol);

#endif // __GRID_HPP
