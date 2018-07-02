## Overview
This library consists mainly of various matrix classes with computational methods like QR factorization, matrix-vector and matrix-matrix multiplication, etc.

## Building


## Usage

## Benchmarking Results
Coordinate sparse matrix storage (array of structs) versus struct of arrays doing matrix-vector multiplication.
<img src="./graphs/AOSvsCOOcomparison.png" alt="AOSvsCOO" width="600px" />

Compressed sparse column versus coordinate sparse (array of structs) storage doing matrix-vect
or multiplication.
<img src="./graphs/CSCvsCOOcomparison.png" alt="CSCvsCOO" width="600px" />
