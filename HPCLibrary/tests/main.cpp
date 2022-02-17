/*
 * main.cpp
 * Description: Example program for linking with libHPCLibrary.a.
 * @author Johnny Sellers
 * @version 0.1 05/30/2017
 */
#include "Matrix.hpp"
#include <typeinfo>

int main(int argc, char *argv[])
{
  // Create a 10-by-10 matrix, randomize, save and read it
  Matrix<int>* A = new Matrix<int>(10,10);
  A->randomizeMatrix();
  A->writeMatrix("exe/A_dense_matrix.txt");
  Matrix<int> R = Matrix<int>::readMatrix("exe/A_dense_matrix.txt");

  // Compute qr factorization and save it to file
  vector<Matrix<int>> QR = A->qr();
  QR[1].writeMatrix("exe/R_upper_triangular.txt");

  // Create a 4-by-4 matrix, randomize, save, create partition and save it
  Matrix<int> B(4, 4);
  B.randomizeMatrix();
  B.writeMatrix("exe/B.txt");
  Matrix<int> B_11 = B.partition(0, 1, 0, 1);
  B_11.writeMatrix("exe/B_11.txt");
  if (!dynamic_cast<Matrix<int>*>(&B_11))
    B_11.writeMatrix(std::cout);

  // Use of +, -, and * operators
  Matrix<int> C(4, 4), D(4, 4);
  C.randomizeMatrix();
  D.randomizeMatrix();
  C.writeMatrix("exe/C.txt");
  D.writeMatrix("exe/D.txt");
  Matrix<int> E = C + D;
  E.writeMatrix("exe/E.txt");
  Matrix<int> F = C - D;
  F.writeMatrix("exe/F.txt");
  Matrix<int> G = C * D;
  G.writeMatrix("exe/G.txt");

  // Test equality with isEqual method
  Matrix<int> Cp = C;
  if (!Cp.isEqual(C))
    std::cout << "Cp != C" << std::endl;
  for (unsigned i = 0; i < C.numRows(); i++)
    for (unsigned j = 0; j < C.numCols(); j++)
      Cp(i,j) = D(i,j);
  if (!Cp.isEqual(D))
    std::cout << "Cp != D" << std::endl;

  // Combine matrices
  vector<vector<Matrix<int>>> v = { { C, D }, { C, D } };
  Matrix<int> J = Matrix<int>::combine(v);
  J.writeMatrix("exe/J.txt");

  // Static multiplication methods
  Matrix<int> H = Matrix<int>::squareMultiply(C, D);
  H.writeMatrix("exe/H.txt");
  Matrix<int> I = Matrix<int>::strassenMultiply(C, D);
  I.writeMatrix("exe/I.txt");

  // PrettyPrint
  Matrix<int> S(4, 4);
  S.randomizeMatrix();
  S.writeMatrix("exe/S.txt");
  if (dynamic_cast<Matrix<int>*>(&S)) {
    S.prettyPrint("exe/S_pp.txt");
    S.prettyPrint(std::cout);
  }

  // Chain multiplication A = A_1 * A_2 *  ... * A_n
  Matrix<int> A_1 = B, A_2 = C, A_3 = D, A_4 = E, A_5 = F;
  vector<Matrix<int>> matrices = { A_1, A_2, A_3, A_4, A_5 };
  vector<int> p = { 0, matrices[0].numRows() };
  for (auto m : matrices)
    p.push_back(m.numCols());
  auto s = Matrix<int>::matrixChainOrder(p);
  Matrix<int>::printOptimalParens(s, 1, 5);
  Matrix<int> W = Matrix<int>::chainMultiply({ A_1, A_2, A_3, A_4, A_5 });
  W.prettyPrint(std::cout);

  return 0;
}
