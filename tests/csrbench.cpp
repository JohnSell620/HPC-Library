/**
		csrbench.cpp
		Description: csrbench.cpp for PS5 -- to compare CSR performance with Roofline model

		@author Johnny Sellers
		@version 1.1 05/10/2017
*/
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <string.h>
#include <functional>
#include "Vector.hpp"
#include "CSR.hpp"
#include "Timer.hpp"

using namespace std;

double benchmark(size_type M, size_type N, long numruns, function<void(const CSRMatrix&, const Vector&, Vector&)>);
void runBenchmark(function<void (const CSRMatrix&, const Vector&, Vector&)>f, long maxsize);

void matvec_CSR(const CSRMatrix& A, const Vector& x, Vector& y) {
  A.matvec(x, y);
}

int main() {

  runBenchmark(matvec_CSR, 8L*32L*32L*8192L);

	return 0;
}

void runBenchmark(function<void (const CSRMatrix&, const Vector&, Vector&)>f, long maxsize) {
	if (strncmp(typeid(size_type).name(), "i", 2) == 0) cout << "integer indexes" << endl;
	else if (strncmp(typeid(size_type).name(), "m", 2) == 0) cout << "unsigned long indexes" << endl;

	cout << "N\tTime\tFLOPs" << endl;
  for (long i = 16; i <= maxsize; i *= 4) {
    long numruns = 4L*4L*1048L*1048L*1048L/(i*i) + 2;
    double t = benchmark(i, i, numruns, f);

		CSRMatrix B(i, i);
		size_type xpts = std::sqrt((double) i);
		B.piscretize(xpts, xpts);
		size_type n = B.numNonzeros();
		double flops = 2*1.e3*numruns*n;

		cout << i << "\t" << t << "\t" << flops / t << endl;
  }
}

double benchmark(size_type M, size_type N, long numruns, function<void (const CSRMatrix&, const Vector&, Vector&)>f) {
  size_type xpoints = std::sqrt((double) M);
  assert(xpoints*xpoints == M);

  CSRMatrix A(M, M);
  Vector x(M), y(M);
  A.piscretize(xpoints, xpoints);
  randomize(x);
  randomize(y);

  Timer T;
  T.start();
  for (size_type i = 0; i < numruns; ++i) {
    f(A, x, y);
  }
  T.stop();

  zeroize(y);

  return T.elapsed();
}
