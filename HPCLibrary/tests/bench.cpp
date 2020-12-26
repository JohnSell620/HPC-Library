/*
 * bench.cpp
 * Description: Benchmarking program to compare CSC and AOS
 * implementations with COO implementation
 * @author Johnny Sellers
 * @version 0.1 04/28/2017
 */
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include "CSC.hpp"
#include "AOS.hpp"
#include "COO.hpp"
#include "Vector.hpp"
#include "Timer.hpp"

using namespace std;

double benchmark_COO(int M, int N, long numruns, function<void(const COOMatrix&, const Vector&, Vector&)>);
void runBenchmark_COO(function<void (const COOMatrix&, const Vector&, Vector&)>f, long maxsize);

double benchmark_AOS(int M, int N, long numruns, function<void(const AOSMatrix&, const Vector&, Vector&)>);
void runBenchmark_AOS(function<void (const AOSMatrix&, const Vector&, Vector&)>f, long maxsize);

double benchmark_CSC(int M, int N, long numruns, function<void(const CSCMatrix&, const Vector&, Vector&)>);
void runBenchmark_CSC(function<void (const CSCMatrix&, const Vector&, Vector&)>f, long maxsize);

void matvec_COO(const COOMatrix& A, const Vector& x, Vector& y) {
  A.matvec(x, y);
}

void matvec_AOS(const AOSMatrix& A, const Vector& x, Vector& y) {
  A.matvec(x, y);
}

void matvec_CSC(const CSCMatrix& A, const Vector& x, Vector& y) {
  A.matvec(x, y);
}

int main(int argc, char *argv[]) {

	if (argc != 2) {
		cout << "Usage: " << argv[0] << " [matrix_type] (\"COOMatrix\", \"AOSMatrix\", \"CSCMatrix\")" << endl;
		return -1;
	}

  if (string(argv[1]) == "COOMatrix")
	  runBenchmark_COO(matvec_COO, 8L*32L*32L*8192L);
	else if (string(argv[1]) == "AOSMatrix")
		runBenchmark_AOS(matvec_AOS, 8L*32L*32L*8192L);
	else if (string(argv[1]) == "CSCMatrix")
		runBenchmark_CSC(matvec_CSC, 8L*32L*32L*8192L);
	else return -2;

	return 0;
}

void runBenchmark_COO(function<void (const COOMatrix&, const Vector&, Vector&)>f, long maxsize) {
	cout << "N\tTperX" << endl;
  for (long i = 16; i <= maxsize; i *= 4) {
    long numruns = 4L*1048L*1048L*1048L/(i*i) + 2;
    double t = benchmark_COO(i, i, numruns, f);

		COOMatrix B(i, i);
		int xpts = std::sqrt((double) i);
		piscretize(B, xpts, xpts);

		cout << i << " " << t/((double)numruns) << endl;
  }
}

double benchmark_COO(int M, int N, long numruns, function<void (const COOMatrix&, const Vector&, Vector&)>f) {
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

void runBenchmark_AOS(function<void (const AOSMatrix&, const Vector&, Vector&)>f, long maxsize) {
	cout << "N\tTperX" << endl;
  for (long i = 16; i <= maxsize; i *= 4) {
    long numruns = 4L*1048L*1048L*1048L/(i*i) + 2;
    double t = benchmark_AOS(i, i, numruns, f);

		AOSMatrix B(i, i);
		int xpts = std::sqrt((double) i);
		piscretize(B, xpts, xpts);

		cout << i << " " << t/((double)numruns) << endl;
  }
}

double benchmark_AOS(int M, int N, long numruns, function<void (const AOSMatrix&, const Vector&, Vector&)>f) {
  int xpoints = std::sqrt((double) M);
  assert(xpoints*xpoints == M);

  AOSMatrix A(M, M);
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

void runBenchmark_CSC(function<void (const CSCMatrix&, const Vector&, Vector&)>f, long maxsize) {
	cout << "N\tTperX" << endl;
  for (long i = 16; i <= maxsize; i *= 4) {
    long numruns = 4L*1048L*1048L*1048L/(i*i) + 2;
    double t = benchmark_CSC(i, i, numruns, f);

		CSCMatrix B(i, i);
		int xpts = std::sqrt((double) i);
		B.piscretize(xpts, xpts);

		cout << i << " " << t/((double)numruns) << endl;
  }
}

double benchmark_CSC(int M, int N, long numruns, function<void (const CSCMatrix&, const Vector&, Vector&)>f) {
  int xpoints = std::sqrt((double) M);
  assert(xpoints*xpoints == M);

  CSCMatrix A(M, M);
  Vector x(M), y(M);
  A.piscretize(xpoints, xpoints);
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
