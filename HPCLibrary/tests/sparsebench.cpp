/*
 * sparsebench.cpp
 * Description: Estimates and prints floating point performance
 * of the sparse matrix-vector product of the COO Matrix class.
 * @author Johnny Sellers
 * @version 0.1 05/30/2017
 */
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include "Vector.hpp"
#include "COO.hpp"
#include "Timer.hpp"

using std::cout; using std::endl; using std::string;


double benchmark_sparse(int M, int N, long numruns, std::function<void(const COOMatrix&, const Vector&, Vector&)>);
void runBenchmark_sparse(std::function<void (const COOMatrix&, const Vector&, Vector&)>f, long maxsize);

void matvec_sparse(const COOMatrix& A, const Vector& x, Vector& y) {
  A.matvec(x, y);
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
		cout << "Usage: " << argv[0] << " [matrix_type]" << endl;
    runBenchmark_sparse(matvec_sparse, 32L*4096L);
		return -1;
	}

  if (string(argv[1]) == "sparse")
    runBenchmark_sparse(matvec_sparse, 32L*4096L);
  else return -2;

  return 0;
}


void runBenchmark_sparse(std::function<void (const COOMatrix&, const Vector&, Vector&)>f, long maxsize) {
  cout << "N\tN*N\tTime\tFlops\tTperX"  << endl;
  for (long i = 16; i <= maxsize; i *= 4) {
    long numruns = 4L*1048L*1048L*1048L/(i*i) + 2;
    double t = benchmark_sparse(i, i, numruns, f);

		COOMatrix B(i, i);
		int xpts = std::sqrt((double) i);
		piscretize(B, xpts, xpts);
		int n = B.numNonzeros();
    //
    // Fill in the next line with the correct formula
    double flops = 2*1.e3*numruns*n;
    //

    cout << i << "\t" << i*i << "\t" << t << "\t" << flops / t << "\t" << t/((double)numruns)  << endl;
  }
}

// int K not used
double benchmark_sparse(int M, int N, long numruns, std::function<void (const COOMatrix&, const Vector&, Vector&)>f) {
  int xpoints = std::sqrt((double) M);
  assert(xpoints*xpoints == M);

  COOMatrix A(M, M);
  Vector x(M), y(M);
  piscretize(A, xpoints, xpoints);
  randomize(x);
  randomize(y);

  Timer T;
  T.start();
  for (int i = 0; i < numruns; ++i) {
    f(A, x, y);
  }
  T.stop();

  zeroize(y);

  return T.elapsed();
}
