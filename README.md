## Overview
This library consists mainly of various matrix classes with computational methods like QR factorization, matrix-vector and matrix-matrix multiplication, etc. Testing MPI and GPU computations is done using CUDA, NVIDIA Container Runtime, and Docker Compose.

Credit to Nikyle Nguyen for the cluster implementation on Alpine Linux using Docker Compose. See project at https://github.com/NLKNguyen/alpine-mpich.

## Building Without Docker Containers
### Install OpenMPI
Some benchmarking programs depend on OpenMPI, but it's not required for most programs. Skip steps 1 through 3 if using these is not desired.
1. Download OpenMPI (recommend extracting contents in /usr/local).
2. Run the following commands.
```bash
$ tar -xzvf openmpi-x.x.x.tar.gz
$ cd openmpi-x.x.x
$ sudo ./configure --prefix=$HOME/openmpi --enable-mpi-cxx
$ sudo make all
$ sudo make install
```
3. In ~/.bashrc file, add the following lines.
```
export PATH=/path/to/openmpi/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/path/to/openmpi/lib\${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

### Building and Using the library
Run the following commands to build the HPCLibrary library.
```bash
$ clone https://github.com/JohnSell620/HPC-Library.git
$ cd HPC-Library/HPCLibrary
$ mkdir obj exe lib
$ make classes
$ ar rcs lib/libHPCLibrary.a obj/*.o
```
The following command links the static library to `main.cpp`:
```bash
$ c++ -std=c++11 -I ./inc -L ./lib -static ./tests/main.cpp -lHPCLibrary -o ./exe/libHPCLibraryClient
```
Now run `$ ./exe/libHPCLibraryClient` to see the output of `main.cpp`. To build the benchmarking tests, just run `$ make`.

Run `$ make precomp_headers` to pre-compile the .hpp files, and include these to optimize programs.

### Pulling Docker Images and Using Docker Compose Cluster
The images jsellers/gpu-mpi:latest and jsellers/gpu-mpi:onbuild can be found at https://hub.docker.com/ (imminent).

## Usage
Try these commands from the HPC-Library/HPCLibrary directory after running `make all`.
```bash
$ ./exe/bench
$ ./exe/csrbench
$ ./exe/sparsebench
```

## Benchmarking Results
Coordinate sparse matrix storage (array of structs) versus struct of arrays doing matrix-vector multiplication.
<img src="./HPCLibrary/graphs/AOSvsCOOcomparison.png" alt="AOSvsCOO" width="600px" />

Compressed sparse column versus coordinate sparse (array of structs) storage doing matrix-vect
or multiplication.
<img src="./HPCLibrary/graphs/CSCvsCOOcomparison.png" alt="CSCvsCOO" width="600px" />


## TODO
1. Fix gpu_densebench.cu timing issue.
2. Correct Matrix::qr() factorization.
3. Correct issues with MPI on Docker Compose.
