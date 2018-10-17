/*
 * Laplacian.hpp
 * Description: 5-point Laplacian stencil class.
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#ifndef __LAPLACIAN_HPP
#define __LAPLACIAN_HPP

#include "Grid.hpp"

class Laplacian { };

Grid operator*(const Laplacian& A, const Grid& x);
double mpiDot(const Grid& X, const Grid& Y);
size_t ir(const Laplacian& A, Grid& x, const Grid& b, size_t max_iter, double tol);
size_t cg(const Laplacian& A, Grid& x, const Grid& b, size_t max_iter, double tol);

#endif // __LAPLACIAN_HPP
