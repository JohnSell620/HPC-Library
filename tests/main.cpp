/*
 * main.cpp
 * Description: Example program for linking with libHPCLibrary.a.
 * @author Johnny Sellers
 * @version 0.1 05/30/2017
 */
#include <iostream>
#include "Matrix.hpp"


int main(int argc, char *argv[])
{
  Matrix A(20,20);
  randomizeMatrix(A);
  writeMatrix(A, "obj/dense_matrix.txt");
  Matrix R = readMatrix("obj/dense_matrix.txt");
  
  qr(A,R);
  writeMatrix(R, "obj/upper_triangular.txt");

  return 0;
}
