/*
* Matrix.cpp
*/
#include "Matrix.hpp"

template<typename T>
T &Matrix<T>::operator()(int i, int j) {
  return this->arrayData[i*this->jCols + j];
}

template<typename T>
const T &Matrix<T>::operator()(int i, int j) const {
  return this->arrayData[i*this->jCols + j];
}

template<typename T>
void Matrix<T>::setValue(int i, int j, T value) {
  this->arrayData[i*(this->jCols) + j] = value;
}

template<typename T>
void Matrix<T>::setValue(int k, T value) {
  this->arrayData[k] = value;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) {
  assert(B.numRows() == jCols);
  int M = B.numRows(), N = B.numCols();
  if (M == N) {
    bool isPowerOf2 = true;
    for (int i = 0; i < 31; ++i) {
      if (N ^ 1 << 1) {
        isPowerOf2 = false;
        break;
      }
    }
    if (isPowerOf2) {
      return Matrix<T>::strassenMultiply(*this, B);
    }
  }
  M = this->numRows();
  Matrix<T> C(M, N);
  for (int i = 0; i < M; ++i)
    for (int j = 0; j < N; ++j)
      for (int k = 0; k < M; ++k)
        C(i, j) += arrayData[i*jCols + k] * B(k,j);
  return C;
}

template<typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& x) {
  assert(x.size() == this->numCols());
  std::vector<T> b (this->numRows());
  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < this->numCols(); ++j)
      b[i] += this->operator()(i, j) * x[j];
  return b;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& B) {
  assert(B.numRows() == iRows && B.numCols() == jCols);
  Matrix<T> C(iRows, jCols);
  for (int i = 0; i < B.numRows(); ++i)
    for (int j = 0; j < B.numCols(); ++j)
      C(i, j) = arrayData[i*jCols + j] + B(i,j);
  return C;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& B) {
  assert(B.numRows() == iRows && B.numCols() == jCols);
  Matrix<T> C(iRows, jCols);
  for (int i = 0; i < B.numRows(); ++i)
    for (int j = 0; j < B.numCols(); ++j)
      C(i, j) = arrayData[i*jCols + j] - B(i,j);
  return C;
}

template<typename T>
void Matrix<T>::operator==(const Matrix<T>& B) {
  assert(B.numRows() == iRows && B.numCols() == jCols);
  for (unsigned i = 0; i < B.numRows(); ++i)
    for (unsigned j = 0; j < B.numCols(); ++j)
      this->arrayData[i*jCols + j] = B(i, j);
}

template<typename T>
bool Matrix<T>::isEqual(const Matrix<T>& B) {
  if (this->numRows() != B.numRows() || this->numCols() != B.numCols())
    return false;
  for (unsigned i = 0; i < iRows; ++i)
    for (unsigned j = 0; j < jCols; ++j)
      if (this->operator()(i,j) != B(i,j))
        return false;
  return true;
}

template<typename T>
Vector Matrix<T>::matvec(const Vector& x) {
  assert(this->numCols() == x.numRows());
  Vector b(this->numRows());
  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < x.numRows(); ++j)
      b(i) += this->operator()(i,j) * x(j);
  return b;
}

template<typename T>
Matrix<T> Matrix<T>::partition(
  int start_row,
  int end_row,
  int start_col,
  int end_col) const {
  assert(
    0 <= start_row &&
    end_row < this->numRows() &&
    0 <= start_col &&
    end_col < this->numCols());
  Matrix<T> A_rc(end_row - start_row + 1, end_col - start_col + 1);
  for (int i = start_row; i <= end_row; ++i)
    for (int j = start_col; j <= end_col; ++j)
      A_rc(i - start_row, j - start_col) = this->operator()(i, j);
  return A_rc;
}

/*
 * Mmay insert 0-elements into combined matrix depending on
 * input matrices' dimensions and alignment.
 */
template<typename T>
Matrix<T> Matrix<T>::combine(vector<vector< Matrix<T> >>& matrices) {
  int m = 0, n = 0;
  for (auto v : matrices) {
    int max_row_dim = 0, total_columns = 0;
    for (auto matrix : v) {
      total_columns += matrix.numCols();
      max_row_dim = matrix.numRows() > max_row_dim ?
                    matrix.numRows() : max_row_dim;
    }
    n = total_columns > n ? total_columns : n;
    m += max_row_dim;
  }
  Matrix<T> C(m, n);
  int ci_start = 0;
  for (int i = 0; i < matrices.size(); i++) {
    int cj_start = 0, max_row_dim = 0;
    for (int j = 0; j < matrices[i].size(); j++) {
      auto submatrix = matrices[i][j];
      max_row_dim = submatrix.numRows() > max_row_dim ?
                    submatrix.numRows() : max_row_dim;
      for (int ii = 0, ci = ci_start; ii < submatrix.numRows(); ii++, ci++) {
        for (int jj = 0, cj = cj_start; jj < submatrix.numCols(); jj++, cj++) {
          C(ci, cj) = submatrix(ii, jj);
        }
      }
      cj_start += submatrix.numCols();
    }
    ci_start += max_row_dim;
  }
  return C;
}

template<typename T>
Matrix<T> Matrix<T>::squareMultiply(const Matrix<T>& A, const Matrix<T>& B) {
  assert(A.numRows() == A.numCols() && A.numRows() == B.numRows());
  int n = A.numRows();
  if (n == 1) {
    Matrix<T> C(1,1);
    C(0,0) = A(0,0) * B(0,0);
    return C;
  }
  Matrix<T> A_11 = A.partition(  0, n/2 - 1,   0, n/2 - 1);
  Matrix<T> A_12 = A.partition(  0, n/2 - 1, n/2,   n - 1);
  Matrix<T> A_21 = A.partition(n/2,   n - 1,   0, n/2 - 1);
  Matrix<T> A_22 = A.partition(n/2,   n - 1, n/2,   n - 1);
  Matrix<T> B_11 = B.partition(  0, n/2 - 1,   0, n/2 - 1);
  Matrix<T> B_12 = B.partition(  0, n/2 - 1, n/2,   n - 1);
  Matrix<T> B_21 = B.partition(n/2,   n - 1,   0, n/2 - 1);
  Matrix<T> B_22 = B.partition(n/2,   n - 1, n/2,   n - 1);
  Matrix<T> C_11 = Matrix::squareMultiply(A_11, B_11) +
                   Matrix::squareMultiply(A_12, B_21);
  Matrix<T> C_12 = Matrix::squareMultiply(A_11, B_12) +
                   Matrix::squareMultiply(A_12, B_22);
  Matrix<T> C_21 = Matrix::squareMultiply(A_21, B_11) +
                   Matrix::squareMultiply(A_22, B_21);
  Matrix<T> C_22 = Matrix::squareMultiply(A_21, B_12) +
                   Matrix::squareMultiply(A_22, B_22);
  vector<vector< Matrix<T> >> C_v = {
    { C_11, C_12 },
    { C_21, C_22 }
  };
  Matrix<T> C = Matrix::combine(C_v);
  return C;
}

template<typename T>
Matrix<T> Matrix<T>::strassenMultiply(const Matrix<T>& A, const Matrix<T>& B) {
  assert(A.numRows() == A.numCols() && A.numRows() == B.numRows());
  int n = A.numRows();
  if (n == 1) {
    Matrix<T> C(1,1);
    C(0,0) = A(0,0) * B(0,0);
    return C;
  }
  Matrix<T> A_11 = A.partition(  0, n/2 - 1,   0, n/2 - 1);
  Matrix<T> A_12 = A.partition(  0, n/2 - 1, n/2,   n - 1);
  Matrix<T> A_21 = A.partition(n/2,   n - 1,   0, n/2 - 1);
  Matrix<T> A_22 = A.partition(n/2,   n - 1, n/2,   n - 1);
  Matrix<T> B_11 = B.partition(  0, n/2 - 1,   0, n/2 - 1);
  Matrix<T> B_12 = B.partition(  0, n/2 - 1, n/2,   n - 1);
  Matrix<T> B_21 = B.partition(n/2,   n - 1,   0, n/2 - 1);
  Matrix<T> B_22 = B.partition(n/2,   n - 1, n/2,   n - 1);
  Matrix<T> S_1  = B_12 - B_22;
  Matrix<T> S_2  = A_11 + A_12;
  Matrix<T> S_3  = A_21 + A_22;
  Matrix<T> S_4  = B_21 - B_11;
  Matrix<T> S_5  = A_11 + A_22;
  Matrix<T> S_6  = B_11 + B_22;
  Matrix<T> S_7  = A_12 - A_22;
  Matrix<T> S_8  = B_21 + B_22;
  Matrix<T> S_9  = A_11 - A_21;
  Matrix<T> S_10 = B_11 + B_12;
  Matrix<T> P_1  = Matrix::strassenMultiply(A_11, S_1);
  Matrix<T> P_2  = Matrix::strassenMultiply(S_2, B_22);
  Matrix<T> P_3  = Matrix::strassenMultiply(S_3, B_11);
  Matrix<T> P_4  = Matrix::strassenMultiply(A_22, S_4);
  Matrix<T> P_5  = Matrix::strassenMultiply(S_5,  S_6);
  Matrix<T> P_6  = Matrix::strassenMultiply(S_7,  S_8);
  Matrix<T> P_7  = Matrix::strassenMultiply(S_9, S_10);
  Matrix<T> C_11 = P_5 + P_4 - P_2 + P_6;
  Matrix<T> C_12 = P_1 + P_2;
  Matrix<T> C_21 = P_3 + P_4;
  Matrix<T> C_22 = P_5 + P_1 - P_3 - P_7;
  vector<vector< Matrix<T> >> C_v = { { C_11, C_12 }, { C_21, C_22 } };
  Matrix<T> C = Matrix<T>::combine(C_v);
  return C;
}

template<typename T>
double Matrix<T>::norm(char type) {
  double value = -1.0;
  switch (type) {
    case '1':
      value = this->oneNorm();
      break;
    case '2':
      value = this->twoNorm();
      break;
    case 'i':
      value = this->infinityNorm();
      break;
    default:
        value = this->twoNorm();
  }
  return value;
}

template<typename T>
Matrix<T> Matrix<T>::outerProduct(const std::vector<T>& x, const std::vector<T>& y) {
  Matrix<T> A (x.size(), y.size());
  for (int i = 0; i < x.size(); ++i)
    for (int j = 0; j < y.size(); ++j)
      A(i, j) = x[i] * y[j];
  return A;
}

// Modified Gram-Schmidt QR Factorization
template<typename T>
vector<Matrix<T>> Matrix<T>::qr() {
  assert(this->numRows() == this->numCols());

  int n = this->numCols();
  Matrix<T> Q = Matrix<T>(n, n);
  Matrix<T> R = Matrix<T>(n, n);

  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      Q(i, j) = this->operator()(i, j);

  vector<T> v (n);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned k = 0; k < v.size(); ++k)
      v[k] = Q(k, i);

    R(i, i) = norm(v);
    for (unsigned k = 0; k < v.size(); ++k)
      Q(k, i) = Q(k, i) / R(i, i);

    for (unsigned j = i + 1; j < this->numCols(); ++j) {
      T r_ij = 0.0;
      for (unsigned k = 0; k < n; ++k)
        r_ij += Q(k, i) * Q(k, j);
      R(i, j) = r_ij;

      for (unsigned k = 0; k < n; ++k)
        Q(k, j) = Q(k, j) - r_ij * Q(k, i);
    }
  }

  return { Q, R };
}

// Householder QR Factorization
template<typename T>
Matrix<T> Matrix<T>::qrHouseholder() {
  assert(this->numRows() == this->numCols());

  vector<T> x   (this->numRows());
  vector<T> v   (this->numRows());
  vector<T> e_1 (this->numRows());
  Matrix<T> R = Matrix<T>(this->numRows(), this->numCols());

  for (unsigned i = 0; i < this->numRows(); ++i)
    for (unsigned j = 0; j < this->numCols(); ++j)
      R(i, j) = this->operator()(i, j);

  double norm_v;
  for (int k = 0; k < this->numCols(); ++k) {
    for (int i = 0; i < this->numRows(); ++i)
      x[i] = this->operator()(i,k);

    if (x[0] < 0)
      e_1[0] = -norm(x);
    else if (x[0] >= 0)
      e_1[0] = norm(x);

    for (unsigned i = 0; i < v.size(); ++i)
      v[i] = e_1[i] + x[i];
    norm_v = norm(v);

    for (int i = 0; i < v.size(); ++i)
      v[i] /= norm_v;

    for (int j = k; j < R.numRows(); ++j)
      R.setValue(k, j, R(k, j) - 2* v[j] * v[j] * R(k, j));
  }
  return R;
}

template<typename T>
Matrix<T> Matrix<T>::chainMultiply(std::initializer_list<Matrix<T>> matrices) {
  vector<Matrix<T>> m;
  m.insert(m.end(), matrices.begin(), matrices.end());
  vector<int> p = { 0 };
  for (auto matrix : matrices)
    p.push_back(matrix.numCols());
  vector<vector<int>> s = matrixChainOrder(p);
  return chainMultiply(m, s, 1, m.size());
}

template<typename T>
vector<vector<int>> Matrix<T>::matrixChainOrder(vector<int>& p) {
  int n = p.size() - 1;
  vector<vector<int>> m (n, vector<int>(n));
  vector<vector<int>> s (n - 1, vector<int>(n - 1));
  for (int l = 2; l <= n; ++l) {
    for (int i = 1; i <= n - l + 1; ++i) {
      int j = i + l - 1;
      m[i-1][j-1] = INT_MAX;
      for (int k = i; k < j; ++k) {
        int q = m[i-1][k-1] + m[k][j-1] + p[i-1] * p[k] * p[j];
        if (q < m[i-1][j-1]) {
          m[i-1][j-1] = q;
          s[i-1][j-1] = k;
        }
      }
    }
  }
  return s;
}

template<typename T>
Matrix<T> Matrix<T>::chainMultiply(
  vector<Matrix<T>>& matrices,
  vector<vector<int>>& s,
  int i,
  int j)
{
  if (i == j) return matrices[i - 1];
  return chainMultiply(matrices, s, i, s[i-1][j-1]) *
         chainMultiply(matrices, s, s[i-1][j-1] + 1, j);
}

template<typename T>
void Matrix<T>::printOptimalParens(vector<vector<int>>& s, int i, int j) {
  if (i == j)
    std::cout << "A_" << i;
  else {
    if (i == 1 && j == s.size())
      std::cout << "optimal order: ";
    std::cout << "(";
    printOptimalParens(s, i, s[i-1][j-1]);
    printOptimalParens(s, s[i-1][j-1] + 1, j);
    std::cout << ")";
    if (i == 1 && j == s.size())
      std::cout << std::endl;
  }
}

template<typename T>
void Matrix<T>::randomizeMatrix() {
  const std::string t = typeid(T).name();
  auto lower_bound = t == "int" ? 2 : 2.0;
  auto upper_bound = t == "int" ? 32 : 32.0;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(lower_bound, upper_bound);
  static auto dice = std::bind(distribution, generator);

  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < this->numCols(); ++j)
      this->operator()(i, j) = dice();
}

template<typename T>
void Matrix<T>::zeroizeMatrix() {
  const std::string t = typeid(T).name();
  auto val = t == "int" ? 0 : 0.0;
  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < this->numCols(); ++j)
      this->operator()(i, j) = val;
}

template<typename T>
void Matrix<T>::overwriteMatrix(std::string fileName) {
  Matrix<T> B = Matrix<T>::readMatrix(fileName);
  this->operator==(B);
}

template<typename T>
Matrix<T> Matrix<T>::readMatrix(std::string fileName) {
  std::string line;
  std::ifstream inputFile (fileName);

  int i_rows, j_cols;
  T value;
  std::string string_input;

  std::getline(inputFile, string_input);
  if (string_input.compare("MATRIX") != 0) {
    std::cout << "Error: incorrect header. Correct header: MATRIX\n";
    std::exit(-2);
  }

  std::getline(inputFile, string_input);
  j_cols = std::stoi(string_input);
  if (j_cols <= 0) {
    std::cout << "Error: columns <= 0?\n";
    std::exit(-2);
  }

  std::getline(inputFile, string_input);
  i_rows = std::stoi(string_input);
  if (i_rows <= 0) {
    std::cout << "Error: rows <= 0?\n";
    std::exit(-2);
  }

  Matrix<T> A = Matrix<T>(i_rows, j_cols);

  for (int i = 0; i < i_rows*j_cols; ++i) {
    std::getline(inputFile, string_input);
    value = std::stod(string_input);
    A.setValue(i, value);
  }

  std::getline(inputFile, string_input);
  if (string_input.compare("END") != 0) {
    std::cout << "Error: incorrect trailer. Correct trailer: END\n";
    std::exit(-2);
  }

  inputFile.close();

  return A;
}

template<typename T>
void Matrix<T>::writeMatrix(std::string file) {
  std::ofstream ofs;
  ofs.open (file, std::ios::out);
  ofs << "MATRIX\n";
  ofs << this->numRows() << std::endl;
  ofs << this->numCols() << std::endl;
  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < this->numCols(); ++j)
      ofs << this->operator()(i, j) << std::endl;
  ofs << "END\n";
}

template<typename T>
void Matrix<T>::writeMatrix(std::ostream& os) {
  os << "MATRIX\n";
  os << this->numRows() << "\n";
  os << this->numCols() << "\n";
  int numrows = this->numRows();
  int numcols = this->numCols();
  for (int i = 0; i < numrows; ++i)
    for (int j = 0; j < numcols; ++j)
      os << this->operator()(i, j) << "\n";
  os << "END\n";
}

template<typename T>
void Matrix<T>::prettyPrint(std::string file) {
  std::ofstream ofs;
  ofs.open (file, std::ios::out);
  ofs << "MATRIX\n";
  ofs << "rows: " << this->numRows() << ", ";
  ofs << "columns: " << this->numCols() << "\n";
  for (int i = 0; i < this->numRows(); ++i) {
    for (int j = 0; j < this->numCols(); ++j) {
      T e = this->operator()(i,j);
      int spaces = 0;
      do { e /= 10; spaces++; } while (e >= 1);

      std::string s = "          "; // 10 spaces (no. of digits in INT_MAX)
      ofs << s.substr(spaces) << this->operator()(i,j) << " ";
    }
    ofs << "\n";
  }
  ofs << "END\n";
}

template<typename T>
void Matrix<T>::prettyPrint(std::ostream& os) {
  os << "MATRIX\n";
  os << "rows: " << this->numRows() << ", ";
  os << "columns: " << this->numCols() << "\n";
  for (int i = 0; i < this->numRows(); ++i) {
    for (int j = 0; j < this->numCols(); ++j) {
      T e = this->operator()(i,j);
      int spaces = 0;
      do { e /= 10; spaces++; } while (e >= 1);

      std::string s = "          "; // 10 spaces (no. of digits in INT_MAX)
      os << s.substr(spaces) << this->operator()(i,j) << " ";
    }
    os << "\n";
  }
  os << "END\n";
}

template<typename T>
double Matrix<T>::twoNorm() {
  double sum = 0.0;
  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < this->numCols(); ++j)
      sum += this->operator()(i, j)*this->operator()(i, j);
  return std::sqrt(sum);
}

template<typename T>
double Matrix<T>::oneNorm() {
  std::vector<double> v(this->numCols(), 0.0);
  for (int j = 0; j < this->numCols(); ++j)
    for (int i = 0; i < this->numRows(); ++i)
      v[j] += std::abs(this->operator()(i, j));
  return *std::max_element(v.begin(), v.end());
}

template<typename T>
double Matrix<T>::infinityNorm() {
  std::vector<double> v(this->numRows(), 0.0);
  for (int i = 0; i < this->numRows(); ++i)
    for (int j = 0; j < this->numCols(); ++j)
      v[i] += std::abs(this->operator()(i, j));
  return *std::max_element(v.begin(), v.end());
}

template<typename T>
double Matrix<T>::norm(const std::vector<T>& v) {
  double sum = 0.0;
  for (T t : v)
    sum += t * t;
  return std::sqrt(sum);
}

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
