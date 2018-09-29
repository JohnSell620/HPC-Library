/**
		Vector.cpp
		Description: Implements Vector.hpp function prototypes.

		@author Johnny Sellers
		@version 0.1 05/10/2017
*/
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <thread>
#include <vector>
#include <mutex>
#include <future>
#include "Vector.hpp"

static std::mutex partial_mutex;

Vector operator+(const Vector& x, const Vector& y) {
  assert(x.numRows() == y.numRows());

  Vector z(x.numRows());
  for (int i = 0; i < x.numRows(); ++i)
    z(i) = x(i) + y(i);

  return z;
}

Vector operator-(const Vector& x, const Vector& y) {
  assert(x.numRows() == y.numRows());

  Vector z(x.numRows());
  for (int i = 0; i < x.numRows(); ++i)
    z(i) = x(i) - y(i);

  return z;
}

Vector operator*(const double& a, const Vector& x) {
  Vector y(x.numRows());
  for (int i = 0; i < x.numRows(); ++i)
    y(i) = a * x(i);

  return y;
}

void zeroize(Vector& v) {
  for (int i = 0; i < v.numRows(); ++i)
    v(i) = 0.0;
}

void randomize(Vector& v) {
  static std::default_random_engine generator;
  static std::uniform_real_distribution<double> distribution(2.0, 32.0);
  static auto dice = std::bind(distribution, generator);

  for (int i = 0; i < v.numRows(); ++i)
    v(i) = dice();
}

double dotProd(const Vector& x, const Vector& y) {
  assert (x.numRows() == y.numRows());
  double sum = 0.0;
  for (int i = 0; i < x.numRows(); ++i)
    sum += x(i) * y(i);

  return sum;
}

double oneNorm(const Vector& v) {
  double sum = 0.0;
  for (int i = 0; i < v.numRows(); ++i)
    sum += std::abs(v(i));

  return sum;
}

double infinityNorm(const Vector& v) {
  double d = 0.0;
  for (int i = 0; i < v.numRows(); ++i)
    d = std::max(d, std::abs(v(i)));

  return d;
}

/* will be affected by roundoff errors for v.size() > ~3000 */
double twoNorm(const Vector& v) {
  double sum = 0.0;
  for (int i = 0; i < v.numRows(); ++i)
    sum += v(i)*v(i);

  return std::sqrt(sum);
}

#ifdef _OPENMP
#include <omp.h>

double ompTwoNorm(const Vector& x) {
	std::string env = getenv("OMP_NUM_THREADS");
	int thread_count;
	if (env.compare("") == 0) thread_count = 1;
	else thread_count = stoi(env);

  double sum = 0.0;
#	pragma omp parallel for reduction(+:sum) num_threads(thread_count)
  for (int i = 0; i < x.numRows(); ++i)
    sum += x(i)*x(i);

  return std::sqrt(sum);
}

#endif

#ifdef _THREADING

void ptn_worker(const Vector& x, size_t begin, size_t end, double& partial) {
	double part_i = 0.0;
	for (size_t i = begin; i < end; ++i)
		part_i += x(i) * x(i);
	{ std::lock_guard<std::mutex> partial_guard (partial_mutex);
		partial += part_i;
	}
}

double partitionedTwoNorm(const Vector& x, size_t partitions) {
	double partial = 0.0;
	size_t part = x.numRows()/partitions;

	std::vector<std::thread> threads;
	for (size_t i = 0; i < partitions; ++i)
		threads.push_back(std::thread(ptn_worker, std::cref(x), i*part, (i+1)*part, std::ref(partial)));

	for (size_t i = 0; i < partitions; ++i)
		threads[i].join();

	return std::sqrt(partial);
}

double recursiveTwoNorm(const Vector& x, size_t levels) {
	return std::sqrt(rtn_worker(x, 0, x.numRows(), levels));
}

double rtn_worker(const Vector& x, size_t begin, size_t end, size_t level) {
	if (level == 0) {
		double sum_squares = 0.0;
		for (size_t i = begin; i < end; ++i)
			sum_squares += x(i) * x(i);

		return sum_squares;
	}
	else {
		std::future<double> sum1 = std::async(std::launch::deferred, rtn_worker, std::cref(x), begin, begin+(end-begin)/2, level-1);
		std::future<double> sum2 = std::async(std::launch::deferred, rtn_worker, std::cref(x), begin+(end-begin)/2, end, level-1);

		return sum1.get() + sum2.get();
	}
}

#endif


Vector readVector(std::istream& inputStream) {
	int iRows;
	float float_input;
	std::string string_input;

	std::getline(inputStream, string_input);
	if (string_input.compare("VECTOR") != 0) {
		std::cout << "Error: incorrect header" << std::endl;
		std::exit(-2);
	}

	std::getline(inputStream, string_input);
	iRows = std::stoi(string_input);
	if (iRows < 0) {
		std:: cout << "Error: vector size < 0?" << std::endl;
		std::exit(-2);
	}

	Vector vector (iRows);

	for (int i=0; i < iRows; ++i)
	{
		std::getline(inputStream, string_input);
		float_input = std::stof(string_input);
		vector.setValue(i, float_input);
	}

	std::getline(inputStream, string_input);

	if (string_input.compare("END") != 0) {
		std::cout << "Error: incorrect trailer" << std::endl;
		std::exit(-2);
	}

	return vector;
}

Vector readVector(std::string fileName) {
	std::string line;
	std::ifstream inputFile (fileName);
	if (inputFile.is_open())
	{
		Vector vector = readVector(inputFile);
		inputFile.close();
		return vector;
	}
	else {
		std::cout << "Error: could not read file" << std::endl;
		std::exit(-3);
	}
}

void writeVector(const Vector& x, std::ostream& os) {
	os << "VECTOR" << std::endl;
	os << x.numRows() << std::endl;

	for (int i = 0; i < x.numRows(); ++i)
		os << x(i) << std::endl;

	os << "END" << std::endl;
}

void writeVector(const Vector& x, std::string file) {
		std::ofstream ofs;
		ofs.open (file, std::ios::out);

		ofs << "VECTOR" << std::endl;
		ofs << x.numRows() << std::endl;
		for (int i = 0; i < x.numRows(); ++i)
			ofs << x(i) << std::endl;

		ofs << "END" << std::endl;

		ofs.close();
}
