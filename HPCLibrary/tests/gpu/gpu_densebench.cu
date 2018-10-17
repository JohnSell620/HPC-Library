/*
 * gpu_densebench.cpp
 * Description: Benchmarking of dense matrix-matrix
 * multiplication on GPU.
 * @author Johnny Sellers
 * @version 1.1 06/11/2017
 */
#include <iostream>
#include <functional>
#include <algorithm>
#include <random>
#include "Matrix.hpp"
#include "Timer.hpp"

#define BLOCK_SIZE 16

class GPUMatrix: public Matrix {
public:
  GPUMatrix(int M, int N):Matrix(M,N) {}
  GPUMatrix():Matrix() {}

  __host__ __device__ std::vector<double>& elements() {
    return arrayData;
  }
};

// Without shared memory
__global__
void MatMulKernel(double *A, double *B, double *C, int Awidth, int Bwidth) {
  double Cvalue = 0;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  for (int k = 0; k < Awidth; ++k)
    Cvalue += A[row*Awidth+k] * B[k*Bwidth+col];
  C[row*Bwidth+col] = Cvalue;
}

double runBench(int M, int N, int K) {
  // Host copies of GPUMatrix A, B, C
  GPUMatrix h_A(M, K), h_B(K, N), h_C(M, N);
  randomizeMatrix(h_A);
  randomizeMatrix(h_B);
  randomizeMatrix(h_C);

  // Copy host matrix elements to arrays for kernel
  double h_Aarr[M*K], h_Barr[K*N], h_Carr[M*N];
  std::copy(h_A.elements().begin(), h_A.elements().end(), h_Aarr);
  std::copy(h_B.elements().begin(), h_B.elements().end(), h_Barr);

  // Allocate space for device copies of A, B elements
  double *d_Aarr;
  size_t size = h_A.numRows()*h_A.numCols()*sizeof(double);
  cudaMalloc((void **)&d_Aarr, size);
  cudaMemcpy(d_Aarr, h_Aarr, size, cudaMemcpyHostToDevice);
  double *d_Barr;
  size = K * N * sizeof(double);
  cudaMalloc((void **)&d_Barr, size);
  cudaMemcpy(d_Barr, h_Barr, size, cudaMemcpyHostToDevice);
  double *d_Carr;
  size = M * N * sizeof(double);
  cudaMalloc((void **)&d_Carr, size);

  // Block and grid dimensions for kernel
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
  dim3 dimGrid(h_B.numCols() / dimBlock.x, h_A.numRows() / dimBlock.y);

  Timer T;
  T.start();
  for (long i = 8; i <= 4096/4; i *= 2) {
    long numruns = 8L*1048L*1048L*16L/(i*i*i) + 2;
    for (int k = 0; k < numruns; ++k)
      MatMulKernel<<<dimGrid, dimBlock>>>(d_Aarr, d_Barr, d_Carr, K, N);
    T.stop();
  }
  T.stop();
  double t = T.elapsed();

  // Copy results back to host GPUMatrix
  cudaMemcpy(h_Carr, d_Carr, size, cudaMemcpyDeviceToHost);
  std::copy(h_Carr, h_Carr + size, std::back_inserter(h_C.elements()));

  // Test print 1-Norm of C
  Matrix& C = h_C;
  std::cout << "1-Norm of C: " << oneNorm(C) << std::endl;

  // Cleanup
  cudaFree(d_Aarr);
  cudaFree(d_Barr);
  cudaFree(d_Carr);

	return t;
}

int main() {
  int dimA = 256, dimB = 16, dimC = 128;
  double t = runBench(dimA, dimB, dimC);
  std::cout << "Achived clockspeed: " << t << std::endl;
  return 0;
}
