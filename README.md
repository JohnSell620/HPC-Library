## Overview
This library consists mainly of various matrix classes with computational methods like QR factorization, matrix-vector and matrix-matrix multiplication, etc. Testing of MPI and GPU computations is done using Nvidia CUDA, MPICH, and Docker Compose.

Credit to Nikyle Nguyen for a cluster implementation model on Alpine Linux using Docker Compose. See [his project here](https://github.com/NLKNguyen/alpine-mpich).

## Using Docker Images
Both Docker and Docker Compose must be installed on the host machine. Then do the following:
```bash
$ git clone https://github.com/JohnSell620/HPC-Library.git
$ cd HPC-Library
$ ./cluster.sh up [size=10]
```
This will pull the Docker images `jhnns/ubuntu-cuda-mpich:latest` and `jhnns/ubuntu-cuda-mpich:onbuild` from [Docker Hub](https://hub.docker.com/r/jhnns/ubuntu-cuda-mpich/).

Use the following command to ssh into the master node:
```bash
$ ssh -o "StrictHostKeyChecking no" -i ssh/id_rsa -p 22 mpi@$(docker inspect -f '{{range .NetworkSettings.Networks}}{{.IPAddress}}{{end}}' hpc-library_master_1)
```

To exit master node and shutdown cluster:
```bash
$ exit
$ ./cluster down
```

## Building Without Docker Containers
### Install OpenMPI
Some benchmarking programs depend on Open MPI, but it's not required for most programs. Skip steps 1 through 3 if using these is not desired.
1. Download Open MPI (extracting contents in `/usr/local` recommended).
2. Run the following command (which may require `sudo`).
```bash
$ wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.3.tar.gz && \
     tar -xzvf openmpi-* && \
     cd openmpi-* && \
     ./configure --prefix=$HOME/openmpi --enable-mpi-cxx && \
     make all && \
     make install
```
3. In ~/.bashrc file, add the following lines.
```
export PATH=/path/to/openmpi/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/path/to/openmpi/lib\${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

### Install CUDA
Installation instruction at [NVIDIA's website](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html).
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

## Usage
Try these commands from the HPC-Library/HPCLibrary directory after running `make all`.
```bash
$ export PATH=/path/to/HPC-Library/HPCLibrary/exe:$PATH
$ bench
$ csrbench
$ sparsebench
```

## Benchmarking Results
Coordinate sparse matrix storage (array of structs) versus struct of arrays doing matrix-vector multiplication.
<img src="./HPCLibrary/graphs/AOSvsCOOcomparison.png" alt="AOSvsCOO" width="600px" />

Compressed sparse column versus coordinate sparse (array of structs) storage doing matrix-vect
or multiplication.
<img src="./HPCLibrary/graphs/CSCvsCOOcomparison.png" alt="CSCvsCOO" width="600px" />


## TODO
1. Fix gpu_densebench.cu timing issue.
2. Add examples (e.g., 2-D Heat Eq.).
