#
# @author Johnny Sellers
# 06/07/2017
# Makefile for HPC-Library.
#
LANG		= -std=c++17 # -Wc++11-extensions

CXX		= c++
GPU		= nvcc
MPI 		= mpic++ # -cxx=clang++ # --enable-mpi-cxx

OPTS		= -Ofast -march=native -DNDEBUG -D_THREADING -pthread
CPPFLAGS	+= -MD -MP
PICKY		= -g -Wall # -Wextra -pedantic

CXXFLAGS	= $(LANG) $(OPTS)
UGTFLAGS	= -lgtest -lgtest_main -lpthread
GPUFLAGS	= -std=c++11
MPIFLAGS	= $(LANG) $(OPTS) # $(PICKY)

INC_PATH	= inc
LIB_PATH	= lib
OBJ_PATH	= obj
MPI_PATH 	= tests/mpi
GPU_PATH	= tests/gpu
EXE_PATH	= exe
GCH_PATH	= gch
UGT_PATH 	= tests/unittests

INCLUDES	= -I ./$(INC_PATH)

VPATH 		:= src tests obj inc
HEADERS		= $(wildcard inc/*.hpp)
TESTS		= $(shell find tests/ ! -name "mpi*.cpp" -name "*.cpp" ! -name "*unittests.cpp")
MPI_TESTS 	= $(wildcard tests/mpi/*.cpp)
GPU_TESTS	= $(wildcard tests/gpu/*.cu)
UNIT_TESTS	= $(wildcard tests/unittests/*.cpp)
SOURCES		= $(wildcard src/*.cpp)

# $(info $$TESTS is [${TESTS}])

CLASSES		= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(SOURCES)))
BENCH 		= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(TESTS)))
MPI_BENCH 	= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(MPI_TESTS)))
GPU_BENCH	= $(patsubst %.cu, $(OBJ_PATH)/%.o, $(notdir $(GPU_TESTS)))
GTEST		= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(UNIT_TESTS)))

OBJECTS		= $(CLASSES) $(BENCH) $(MPI_BENCH) $(GPU_BENCH) $(GTEST)
TARGETS		= $(notdir $(BENCH))
MPI_TARGETS	= $(notdir $(MPI_BENCH))
GPU_TARGETS	= $(notdir $(GPU_BENCH))
GTEST_TARGETS 	= $(notdir $(GTEST))
PCHS		= $(notdir $(HEADERS:=.gch))

EXECS		= $(TARGETS:.o=)
MPI_EXECS	= $(MPI_TARGETS:.o=)
GPU_EXECS	= $(GPU_TARGETS:.o=)
GTEST_EXECS	= $(GTEST_TARGETS:.o=)

all		 : $(OBJECTS) $(EXECS) $(MPI_EXECS) $(GPU_EXECS) $(GTEST_EXECS)

classes 	 : $(CLASSES)

precomp_headers : $(PCHS)

$(EXECS)	 :
		  $(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o  $(EXE_PATH)/$@

$(MPI_EXECS)	 :
		  $(MPI) $(MPIFLAGS) $(INCLUDES) $^ -o $(EXE_PATH)/$@

$(GPU_EXECS)	 :
		  $(GPU) $(GPUFLAGS) $(INCLUDES) $^ -o $(EXE_PATH)/$@

$(GTEST_EXECS)	 :
		  $(CXX) $(LANG) $(OPTS) $(INCLUDES) $^ $(UGTFLAGS) -o $(EXE_PATH)/$@

$(OBJ_PATH)/%.o : %.cpp
		  $(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(OBJ_PATH)/Laplacian.o : Laplacian.cpp
		  $(MPI) -c $(MPIFLAGS) $(INCLUDES) $< -o $@

$(OBJ_PATH)/%.o : $(MPI_PATH)/%.cpp
		  $(MPI) -c $(MPIFLAGS) $(INCLUDES) $< -o $@

$(OBJ_PATH)/%.o : $(GPU_PATH)/%.cu
		  $(GPU) -c $(GPUFLAGS) $(INCLUDES) $< -o $@

$(OBJ_PATH)/%.o : $(UGT_PATH)/%.cpp
		  $(CXX) -c $(LANG) $(OPTS) $(UGTFLAGS) $(INCLUDES) $< -o $@

%.hpp.gch	: %.hpp
		  $(CXX) $(LANG) -x c++-header $< -o $(GCH_PATH)/$@

.PHONY:	clean

clean:
	/bin/rm -f obj/* exe/* gch/*

bench			: bench.o AOS.o COO.o CSC.o Vector.o
csrbench		: csrbench.o CSR.o Vector.o
densebench		: densebench.o Matrix.o Vector.o
main			: main.o Matrix.o Vector.o
sparsebench		: sparsebench.o COO.o Vector.o
mpi2norm_driver		: mpi2norm_driver.o Vector.o
mpi2norm_timer		: mpi2norm_timer.o Vector.o
mpiHeatEq2D-CG		: mpiHeatEq2D-CG.o Grid.o Laplacian.o
mpiHeatEq2D-IR		: mpiHeatEq2D-IR.o Grid.o Laplacian.o
mpi_hello_world		: mpi_hello_world.o
gpu_densebench		: gpu_densebench.o Matrix.o Vector.o
unittests		: unittests.o CSC.o Matrix.o Vector.o

bench.o			: bench.cpp AOS.hpp CSC.hpp COO.hpp Vector.hpp
csrbench.o		: csrbench.cpp CSR.hpp Vector.hpp
densebench.o		: densebench.cpp Matrix.hpp
main.o			: main.cpp Matrix.hpp
sparsebench.o		: sparsebench.cpp COO.hpp Vector.hpp
mpi2norm_timer.o	: mpi/mpi2norm_timer.cpp Vector.hpp
mpi2norm_driver.o	: mpi/mpi2norm_driver.cpp Vector.hpp
mpiHeatEq2D-CG.o	: mpi/mpiHeatEq2D-CG.cpp Grid.hpp Laplacian.hpp
mpiHeatEq2D-IR.o	: mpi/mpiHeatEq2D-IR.cpp Laplacian.hpp
mpi_hello_world.o	: mpi/mpi_hello_world.cpp
gpu_densebench.o	: gpu/gpu_densebench.cu
unittests.o		: unittests/unittests.cpp Matrix.hpp Vector.hpp
AOS.o			: AOS.hpp Vector.hpp
COO.o			: COO.hpp Vector.hpp
CSC.o			: CSC.hpp Vector.hpp
CSR.o			: CSR.hpp Vector.hpp
Matrix.o		: Matrix.hpp Vector.hpp
Vector.o		: Vector.hpp
Laplacian.o		: Laplacian.hpp
Grid.o			: Grid.hpp
