/*
 * Vector.hpp
 * Description: Vector class.
 * @author Johnny Sellers
 * @version 0.1 05/10/2017
 */
#if !defined(VECTOR_HPP)
#define VECTOR_HPP

#include <vector>
#include <string>

typedef unsigned long size_t;

// Prototype for wrapper function to read environment variables
std::string getenv_to_string(const char *in);
std::string getenv(const std::string& in);

class Vector {
  typedef std::vector<double>::size_type size_type;

public:
  explicit Vector(size_type M) : iRows(M), arrayData(iRows      ) {}
  explicit Vector(size_type M, double init) : iRows(M), arrayData(iRows, init) {}

  ~Vector() {}

  double &operator()(size_type i) { return arrayData[i]; }
  const double &operator()(size_type i) const { return arrayData[i]; }

  size_type numRows() const { return arrayData.size(); }

  // Used in readVector()
  void setValue(int k, double value) {
		arrayData[k] = value;
	}

private:
  size_type iRows;
  std::vector<double> arrayData;
};

Vector operator+(const Vector& x, const Vector& y);
Vector operator-(const Vector& x, const Vector& y);
Vector operator*(const double& a, const Vector& x);
void zeroize(Vector& x);
void randomize(Vector& x);
double oneNorm(const Vector& x);
double twoNorm(const Vector& x);
double infinityNorm(const Vector& x);
double dotProd(const Vector& x, const Vector& y);
double ompTwoNorm(const Vector& x);
double partitionedTwoNorm(const Vector& x, size_t partitions);
void ptn_worker(const Vector& x, size_t begin, size_t end, double& partial);
double recursiveTwoNorm(const Vector& x, size_t levels);
double rtn_worker(const Vector& x, size_t begin, size_t end, size_t level);
Vector readVector(std::istream& inputStream);
Vector readVector(std::string fileName);
void writeVector(const Vector& x, std::ostream& os);
void writeVector(const Vector& x, std::string file);

#endif // VECTOR_HPP
