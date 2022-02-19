#
# @author Johnny Sellers
# 06/07/2017
# Makefile for HPC-Library.
#
$(info *** Pass in MPI=1 and/or GPU=1 on machines with MPI and/or GPU enabled. ***)

# Compilers
CVER = -std=c++17# g++ (default compiler)
NVCC = nvcc
MPIC = mpic++

# Debug options
DEBUG = 0
DEBUGPICKY = 0
ifeq ($(DEBUGPICKY),1)
	DEBUG=1
endif

# Set compiler flags
CXXFLAGS = $(CVER)
ifeq ($(DEBUG),1)
CXXFLAGS += -g3
ifeq ($(DEBUGPICKY),1)
CXXFLAGS += -Wall -Wextra -pedantic
ifeq ($(CVER),-std=c++11)
CXXFLAGS += -Wc++11-extensions
endif
endif
else
	CXXFLAGS += -Ofast -march=native -DNDEBUG
endif

CXXFLAGS	+= -D_THREADING -pthread
UGTFLAGS	= -lgtest -lgtest_main -lpthread
GPUFLAGS	= -std=c++11
MPIFLAGS	= $(CVER) $(CXXFLAGS)# -cxx=clang++ # --enable-mpi-cxx

# Object targets
CLASS_OBJECTS = $(filter-out obj/Laplacian.o, $(patsubst %.cpp, obj/%.o, $(notdir $(wildcard src/*.cpp))))
CPU_OBJECTS	= $(patsubst %.cpp, obj/%.o, $(notdir $(shell find tests/ ! -name "mpi*.cpp" -name "*.cpp" ! -name "*unittests.cpp")))
MPI_OBJECTS	= $(patsubst %.cpp, obj/%.o, $(notdir $(wildcard tests/mpi/*.cpp)))
GPU_OBJECTS	= $(patsubst %.cu, obj/%.o, $(notdir $(wildcard tests/gpu/*.cu)))
UTEST_OBJECTS 	= $(patsubst %.cpp, obj/%.o, $(notdir $(wildcard tests/unittests/*.cpp)))
HEADERS 	= $(wildcard inc/*.hpp)

# Executable targets
CPU_EXECS	= $(notdir $(CPU_OBJECTS:.o=))
MPI_EXECS	:= $(notdir $(MPI_OBJECTS:.o=))
GPU_EXECS	= $(notdir $(GPU_OBJECTS:.o=))
UTEST_EXECS	= $(notdir $(UTEST_OBJECTS:.o=))
PCHS		= $(notdir $(HEADERS:=.gch))

# Is MPI enabled on machine
MPI = 0
ifeq ($(MPI),1)
MPI_OBJECTS += obj/Laplacian.o
else
MPI_OBJECTS =
MPI_EXECS =
endif

# IS GPU enabled on machine
GPU = 0
ifeq ($(GPU),0)
GPU_OBJECTS =
GPU_EXECS =
endif

OBJECTS = $(CLASS_OBJECTS) $(CPU_OBJECTS) $(UTEST_OBJECTS) $(MPI_OBJECTS) $(GPU_OBJECTS)
EXECUTABLES = $(CPU_EXECS) $(UTEST_EXECS) $(MPI_EXECS) $(GPU_EXECS)

.PHONY: all, cpu, mpi, gpu, clean
all 	: $(EXECUTABLES)
cpu 	: $(CPU_EXECS)
mpi 	: $(MPI_EXECS)
gpu 	: $(GPU_EXECS)
classes	: $(notdir $(CLASS_OBJECTS)) $(notdir $(findstring obj/Laplacian.o, $(MPI_OBJECTS)))
objects	: $(notdir $(OBJECTS))
precompiled_headers : $(PCHS)
clean	:
	/bin/rm -f obj/* exe/* gch/*

# Specify directories to search
VPATH := src tests obj inc

# Source directories
INC_PATH = inc
OBJ_PATH = obj
MPI_PATH = tests/mpi
GPU_PATH = tests/gpu
EXE_PATH = exe
UGT_PATH = tests/unittests
GCH_PATH = gch

# Executables
$(CPU_EXECS)	 :
		  $(CXX) $(CXXFLAGS) -I$(INC_PATH) $^ -o $(EXE_PATH)/$@

$(MPI_EXECS)	 :
		  $(MPIC) $(MPIFLAGS) -I$(INC_PATH) $^ -o $(EXE_PATH)/$@

$(GPU_EXECS)	 :
		  $(NVCC) $(GPUFLAGS) -I$(INC_PATH) $^ -o $(EXE_PATH)/$@

# Googletest related flags $(UGTFLAGS), -lgtest, -lgtest_main, -lpthread flags must be added last.
$(UTEST_EXECS)	 :
		  $(CXX) $(CXXFLAGS) -I$(INC_PATH) $^ $(UGTFLAGS) -o $(EXE_PATH)/$@

# Objects
$(CLASS_OBJECTS) :
	$(CXX) -c $(CXXFLAGS) -I$(INC_PATH) $< -o $@

$(CPU_OBJECTS) :
		  $(CXX) -c $(CXXFLAGS) -I$(INC_PATH) $< -o $@

$(MPI_OBJECTS) :
		  $(MPIC) -c $(MPIFLAGS) -I$(INC_PATH) $< -o $@

$(GPU_OBJECTS) :
		  $(NVCC) -c $(GPUFLAGS) -I$(INC_PATH) $< -o $@

$(UTEST_OBJECTS) :
		  $(CXX) -c $(CXXFLAGS) $(UGTFLAGS) -I$(INC_PATH) $< -o $@

# Precompiled headers
%.hpp.gch	: %.hpp
		  $(CXX) $(CVER) -x c++-header $< -o $(GCH_PATH)/$@

# Executable targets
bench			: bench.o AOS.o COO.o CSC.o Vector.o
csrbench		: csrbench.o CSR.o Vector.o
densebench		: densebench.o Matrix.o Vector.o
main			: main.o  Matrix.o
sparsebench		: sparsebench.o COO.o Vector.o
mpi2norm_driver		: mpi2norm_driver.o Vector.o
mpi2norm_timer		: mpi2norm_timer.o Vector.o
mpiHeatEq2D-CG		: mpiHeatEq2D-CG.o Grid.o Laplacian.o
mpiHeatEq2D-IR		: mpiHeatEq2D-IR.o Grid.o Laplacian.o
mpi_hello_world		: mpi_hello_world.o
gpu_densebench		: gpu_densebench.o Matrix.o Vector.o
unittests		: unittests.o CSC.o Matrix.o Vector.o

# Object targets
bench.o		: bench.cpp AOS.cpp CSC.cpp COO.cpp Vector.cpp
csrbench.o		: csrbench.cpp CSR.cpp Vector.cpp
densebench.o	: densebench.cpp Matrix.cpp Vector.cpp
main.o		: main.cpp Matrix.cpp
sparsebench.o	: sparsebench.cpp COO.cpp Vector.cpp
mpi2norm_timer.o	: $(MPI_PATH)/mpi2norm_timer.cpp Vector.cpp
mpi2norm_driver.o	: $(MPI_PATH)/mpi2norm_driver.cpp Vector.cpp
mpiHeatEq2D-CG.o	: $(MPI_PATH)/mpiHeatEq2D-CG.cpp Grid.cpp Laplacian.cpp
mpiHeatEq2D-IR.o	: $(MPI_PATH)/mpiHeatEq2D-IR.cpp Laplacian.cpp
mpi_hello_world.o	: $(MPI_PATH)/mpi_hello_world.cpp
gpu_densebench.o	: $(GPU_PATH)/gpu_densebench.cu
unittests.o		: $(UGT_PATH)/unittests.cpp CSC.cpp Matrix.cpp Vector.cpp
AOS.o		: AOS.cpp
COO.o		: COO.cpp
CSC.o		: CSC.cpp
CSR.o		: CSR.cpp
Matrix.o		: Matrix.cpp
Vector.o		: Vector.cpp
Laplacian.o		: Laplacian.cpp
Grid.o		: Grid.cpp
