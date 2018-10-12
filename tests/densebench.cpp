/*
 * densebench.cpp
 * Description: Benchmarking of dense matrix-matrix
 * multiplication.
 * @author Johnny Sellers
 * @version 1.1 04/25/2017
 */
#include <iostream>
#include <functional>
#include "Matrix.hpp"
#include "Timer.hpp"

double runBench(std::function<void (const Matrix&, const Matrix&, Matrix&)> f);
double bench(int M, int N, int K, long numruns, std::function<void(const Matrix&, const Matrix&, Matrix&)> f);

int main() {

	double clockspeed			= 3.4e9;
	double rate_h12		    = clockspeed*8/12;
	double rate_h20				= clockspeed*8/20;
	double rate_cbh12		  = clockspeed*8/12;
	double rate_cbh20			= clockspeed*8/20;

	double achieved_h 	= runBench(hoistedTiledMultiply2x2);
	double achieved_cbh = runBench(hoistedCopyBlockedTiledMultiply2x2);


	std::cout << "Routine Clock CPUID Loop-ops Scalar 2-wide 4-wide 4-wide-fma Achieved" << std::endl;

	std::cout << "hoisted " << clockspeed << " AVX2 12 " << rate_h12 << " " << 2*rate_h12 << " " << 4*rate_h12 << " " << 8*rate_h12 << " " << achieved_h << std::endl;
	std::cout << "hoisted " << clockspeed << " AVX2 20 " << rate_h20 << " " << 2*rate_h20 << " " << 4*rate_h20 << " " << 8*rate_h20 << " " << achieved_h << std::endl;
	std::cout << "copyblockhoisted " << clockspeed << " AVX2 12 " << rate_cbh12 << " " << 2*rate_cbh12 << " " << 4*rate_cbh12 << " " << 8*rate_cbh12 << " " << achieved_cbh << std::endl;
	std::cout << "copyblockhoisted " << clockspeed << " AVX2 20 " << rate_cbh20 << " " << 2*rate_cbh20 << " " << 4*rate_cbh20 << " " << 8*rate_cbh20 << " " << achieved_cbh << std::endl;

	return 0;

}


double runBench(std::function<void (const Matrix&, const Matrix&, Matrix&)> f) {
	double a = 0.0;
  for (long i = 8; i <= 4096/4; i *= 2) {
		long numruns = 8L*1048L/(i*i*i) + 2;
    double t = bench(i, i, i, numruns, f);
    double flops_per_multiply = i*i*i;
		if (a < 2.0*1.e3*numruns*flops_per_multiply/t)
			a = 2.0*1.e3*numruns*flops_per_multiply/t;
  }

	return a;
}

double bench(int M, int N, int K, long numruns, std::function<void (const Matrix&, const Matrix&, Matrix&)> f) {
  Matrix A(M, K), B(K, N), C(M, N);
  randomizeMatrix(A);
  randomizeMatrix(B);
  randomizeMatrix(C);

  Timer T;
  T.start();
  for (int i = 0; i < numruns; ++i) {
    f(A, B, C);
  }
  T.stop();

  zeroizeMatrix(C);

  return T.elapsed();
}
