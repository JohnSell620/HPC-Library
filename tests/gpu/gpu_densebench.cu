/**
		gpu_densebench.cpp
		Description: Benchmarking of dense matrix-matrix multiplication on GPU.

		@author Johnny Sellers
		@version 1.1 06/11/2017
*/
#include <iostream>
#include <functional>
#include "Matrix.hpp"
#include "Timer.hpp"

__global__
void matmat_multiply(const Matrix *A, const Matrix *B, Matrix *C)
{
  multiply(A, B, C);
}


double runBench(int M, int N, int K, long numruns);

int main() {

	double clockspeed			= 3.4e9;
	double rate_h12		    = clockspeed*8/12;
	double rate_h20				= clockspeed*8/20;
	double rate_cbh12		  = clockspeed*8/12;
	double rate_cbh20			= clockspeed*8/20;

  double achieved_gpu = runBench(matmat_gpu);


	std::cout << "Routine Clock CPUID Loop-ops Scalar 2-wide 4-wide 4-wide-fma Achieved" << std::endl;

	std::cout << "hoisted " << clockspeed << " AVX2 12 " << rate_h12 << " " << 2*rate_h12 << " " << 4*rate_h12 << " " << 8*rate_h12 << " " << achieved_h << std::endl;
	std::cout << "hoisted " << clockspeed << " AVX2 20 " << rate_h20 << " " << 2*rate_h20 << " " << 4*rate_h20 << " " << 8*rate_h20 << " " << achieved_h << std::endl;
	std::cout << "copyblockhoisted " << clockspeed << " AVX2 12 " << rate_cbh12 << " " << 2*rate_cbh12 << " " << 4*rate_cbh12 << " " << 8*rate_cbh12 << " " << achieved_cbh << std::endl;
	std::cout << "copyblockhoisted " << clockspeed << " AVX2 20 " << rate_cbh20 << " " << 2*rate_cbh20 << " " << 4*rate_cbh20 << " " << 8*rate_cbh20 << " " << achieved_cbh << std::endl;

	return 0;

}


double runBench(int M, int N, int K, long numruns) {
  Matrix A(M, K), B(K, N), C(M, N);           // host copies
  Matrix *d_A(M,K), *d_B(K, N), *d_C(M, N);   // device copies

  // Allocate space for device copies
  cudaMalloc((void **)&d_A, size_A);
  cudaMalloc((void **)&d_A, size_B);
  cudaMalloc((void **)&d_A, size_C);

  randomizeMatrix(A);
  randomizeMatrix(B);
  randomizeMatrix(C);

  // Copy inputs to device
  cudaMemcpy(d_A, &A, size_A, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, &B, size_B, cudaMemcpyHostToDevice);
  cudaMemcpy(d_C, &C, size_C, cudaMemcpyHostToDevice);

	double a = 0.0;
  for (long i = 8; i <= 4096/4; i *= 2) {
		long numruns = 8L*1048L/(i*i*i) + 2;
    // double t = bench(i, i, i, numruns, f);
    Timer T;
    T.start();
    for (int i = 0; i < numruns; ++i) {
      // f(A, B, C);
      matmat_multiply<<<1,1>>>(d_A, d_B, d_C);
    }
    T.stop();

    double t = T.elapsed();
    double flops_per_multiply = i*i*i;
		if (a < 2.0*1.e3*numruns*flops_per_multiply/t)
			a = 2.0*1.e3*numruns*flops_per_multiply/t;
  }

  // Copy results back to host
  cudaMemcpy(&C, d_C, size_C, cudaMemcpyDeviceToHost);

  // Cleanup
  cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);

	return a;
}
