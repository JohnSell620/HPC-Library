/*
* unittests.cpp
* Description: Unit tests for Matrix, Vector, and
* CSCMatrix, using static library libHPCLibrary.a.
* @author Johnny Sellers
* @version 0.1 07/03/2018
*/
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "Matrix.hpp"
#include "CSC.hpp"

#ifndef _THREADING
#define _THREADING
#endif


namespace {

struct MatrixTest : public ::testing::Test {
    Matrix<int> * A = new Matrix<int>(20,20);

    void SetUp() {
        A->randomizeMatrix();
        A->writeMatrix("exe/test_write_A.txt");
    }

    void TearDown() {
        delete A;
    }
};

struct CSCMatrixTest : public testing::Test {
    CSCMatrix      * B = new CSCMatrix(16,16);
    Vector         * x = new Vector(16);
    Vector         * y = new Vector(16);

    void SetUp() {
        B->piscretize(4,4);
        randomize(*x);
        randomize(*y);
    }

    void TearDown() {
        delete B;
        delete x;
        delete y;
    }
};

struct VectorTest : public testing::Test {
    Vector * x = new Vector(16);
    void SetUp() { randomize(*x); }
    void TearDown() { /* delete x; */ }
};

TEST_F(MatrixTest, writeMatrixTest) {
    // Write A std::cout to string
    std::stringstream buffer;
    std::streambuf *old = std::cout.rdbuf(buffer.rdbuf());
    A->writeMatrix(buffer);
    std::string stext = buffer.str();
    std::cout.rdbuf(old);

    // Write A to file "test_write_A.txt" then to string
    A->writeMatrix("exe/test_write_A.txt");
    std::ifstream f {"exe/test_write_A.txt"};
    std::string ftext {(std::istreambuf_iterator<char>(f)),
    std::istreambuf_iterator<char>()};

    // Compare strings
    ASSERT_STREQ(stext.c_str(), ftext.c_str());
}

TEST_F(MatrixTest, qrTest) {
    vector<Matrix<int>> QR = A->qr();
    auto R = QR[1];
    double e;
    for (int j = 0; j < A->numCols(); ++j) {
        for (int i = j+1; i < A->numRows(); ++i) {
            e = R(i,j);
            ASSERT_EQ(e, 0.0);
        }
    }
}

TEST_F(VectorTest, Vector2NormTest) {
    Vector const& const_x = *x;
    size_t partitions = 2;
    double ansp = partitionedTwoNorm(const_x, 2);
    double ansr = recursiveTwoNorm(const_x, 2);
}

TEST_F(CSCMatrixTest, MatvecWriteTest) {
    B->matvec(*x, *y);

    std::stringstream buffer;
    std::streambuf *old = std::cout.rdbuf(buffer.rdbuf());
    writeVector(*y, buffer);
    std::string stext = buffer.str();
    std::cout.rdbuf(old);

    writeVector(*y, "exe/test_write_y.txt");
    std::ifstream f {"exe/test_write_y.txt"};
    std::string ftext {
        (std::istreambuf_iterator<char>(f)),
        std::istreambuf_iterator<char>()
    };

    ASSERT_STREQ(stext.c_str(), ftext.c_str());
}

} // namespace

int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
