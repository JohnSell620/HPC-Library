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

double runBench(
	std::function<Matrix<int> (const Matrix<int>&, const Matrix<int>&)> fptr);
double bench(
	int M,
	int N,
	int K,
	long numruns,
	std::function<Matrix<int> (const Matrix<int>&, const Matrix<int>&)> fptr);

int main() {

	double clockspeed			= 3.4e9;
	double rate_h12		    = clockspeed*8/12;
	double rate_h20				= clockspeed*8/20;
	double rate_cbh12		  = clockspeed*8/12;
	double rate_cbh20			= clockspeed*8/20;

	double achieved_h 	= runBench(&Matrix<int>::squareMultiply);
	double achieved_cbh = runBench(&Matrix<int>::strassenMultiply);

	std::cout <<
		"Routine Clock CPUID Loop-ops Scalar "
		"2-wide 4-wide 4-wide-fma Achieved" << std::endl;

	std::cout << "hoisted " << clockspeed << " AVX2 12 " 	<< rate_h12
						<< " " 				<< 2*rate_h12 << " " 					<< 4*rate_h12
						<< " " 				<< 8*rate_h12 << " " 					<< achieved_h
						<< std::endl;
	std::cout << "hoisted " << clockspeed << " AVX2 20 " 	<< rate_h20
						<< " " 				<< 2*rate_h20 << " " 					<< 4*rate_h20
						<< " " 				<< 8*rate_h20 << " " 					<< achieved_h
						<< std::endl;
	std::cout << "copyblockhoisted " << clockspeed 	 << " AVX2 12 " << rate_cbh12
						<< " " 								 << 2*rate_cbh12 << " " 				<< 4*rate_cbh12
						<< " " 								 << 8*rate_cbh12 << " " 				<< achieved_cbh
						<< std::endl;
	std::cout << "copyblockhoisted " << clockspeed   << " AVX2 20 "  << rate_cbh20
						<< " " 								 << 2*rate_cbh20 << " " 				 << 4*rate_cbh20
						<< " " 								 << 8*rate_cbh20 << " " 				 << achieved_cbh
						<< std::endl;

	return 0;

}

double runBench(
	std::function<Matrix<int> (const Matrix<int>&, const Matrix<int>&)> fptr)
{
	double exec_time = 0.0;
  for (long i = 8; i <= 4096 / 16; i *= 2) {
		long numruns = 8L * 1048L / (i * i * i) + 2;
    double t = bench(i, i, i, numruns, fptr);
    double flops_per_multiply = i * i * i;
		if (exec_time < 2.0 * 1.e3 * numruns * flops_per_multiply / t)
			exec_time = 2.0 * 1.e3 * numruns * flops_per_multiply / t;
  }
	return exec_time;
}

double bench(
	int M,
	int N,
	int K,
	long numruns,
	std::function<Matrix<int> (const Matrix<int>&, const Matrix<int>&)> fptr)
{
  Matrix<int> A(M, K), B(K, N), C(M, N);
  A.randomizeMatrix();
  B.randomizeMatrix();

  Timer T;
  T.start();

  for (int i = 0; i < numruns; ++i)
    C = fptr(std::as_const(A), std::as_const(B));

  T.stop();

  C.zeroizeMatrix();

  return T.elapsed();
}
