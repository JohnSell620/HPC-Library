#
# @author Johnny Sellers
# 06/07/2017
# Makefile for HPC-Library.
#
$(info *** Pass in MPIC=1 and/or CUDA=1 on machines with MPI and/or GPU enabled. ***)

# Compilers
GPU		= nvcc
MPI 		= mpic++

# g++ (default) compiler version
CVER		= -std=c++17

# Debug options
DEBUG		= 0
DEBUGPICKY	= 0
ifeq ($(DEBUGPICKY),1)
	DEBUG=1
endif

# Set compiler flags
CXXFLAGS	= $(CVER)
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

# Source directories
INC_PATH	= inc
LIB_PATH	= lib
OBJ_PATH	= obj
MPI_PATH 	= tests/mpi
GPU_PATH	= tests/gpu
EXE_PATH	= exe
GCH_PATH	= gch
UGT_PATH 	= tests/unittests

# Specify directories to search
VPATH 		:= src tests obj inc

# Object file targets
CLASS_TARGETS = $(filter-out Laplacian.o, $(patsubst %.cpp, %.o, $(notdir $(wildcard src/*.cpp))))
CPU_TARGETS	= $(patsubst %.cpp, %.o, $(notdir $(shell find tests/ ! -name "mpi*.cpp" -name "*.cpp" ! -name "*unittests.cpp")))
MPI_TARGETS	= $(patsubst %.cpp, %.o, $(notdir $(wildcard tests/mpi/*.cpp)))
GPU_TARGETS	= $(patsubst %.cu, %.o, $(notdir $(wildcard tests/gpu/*.cu)))
UTEST_TARGETS 	= $(patsubst %.cpp, %.o, $(notdir $(wildcard tests/unittests/*.cpp)))

# Specify executable targets
CPU_EXECS	= $(notdir $(CPU_TARGETS:.o=))
MPI_EXECS	= $(notdir $(MPI_TARGETS:.o=))
GPU_EXECS	= $(notdir $(GPU_TARGETS:.o=))
UTEST_EXECS	= $(notdir $(UTEST_TARGETS:.o=))
PCHS		= $(notdir $(HEADERS:=.gch))

OBJECTS		= $(CPU_TARGETS) $(UTEST_TARGETS)
EXECUTABLES = $(CPU_EXECS) $(UTEST_EXECS)

# If MPI is enabled on machine
MPIC		= 0
ifeq ($(MPIC),1)
	CLASS_TARGETS += Laplacian.o
	OBJECTS		+= $(CLASS_TARGETS) $(MPI_TARGETS)
	EXECUTABLES	 += $(MPI_EXECS)
else
	OBJECTS		+= $(filter-out Laplacian.o, $(CLASS_TARGETS))
endif

# If CUDA is enabled on machine
CUDA		= 0
ifeq ($(CUDA),1)
	OBJECTS		+= $(GPU_TARGETS)
	EXECUTABLES	+= $(GPU_EXECS)
endif

.PHONY: all, clean
all		 : $(EXECUTABLES)
classes : $(CLASS_TARGETS)
precompiles_headers : $(PCHS)
objects	 : $(OBJECTS)
clean 	 :
	/bin/rm -f obj/* exe/* gch/*

# Executable rules
$(CPU_EXECS)	 :
		  $(CXX) $(CXXFLAGS) -I$(INC_PATH) $^ -o $(EXE_PATH)/$@

$(MPI_EXECS)	 :
		  $(MPI) $(MPIFLAGS) -I$(INC_PATH) $^ -o $(EXE_PATH)/$@

$(GPU_EXECS)	 :
		  $(GPU) $(GPUFLAGS) -I$(INC_PATH) $^ -o $(EXE_PATH)/$@

# Googletest related flags $(UGTFLAGS), -lgtest, -lgtest_main, -lpthread flags must be added last.
$(UTEST_EXECS)	 :
		  $(CXX) $(CXXFLAGS) -I$(INC_PATH) $^ $(UGTFLAGS) -o $(EXE_PATH)/$@

# Object rules
$(OBJ_PATH)/%.o : %.cpp
		  $(CXX) -c $(CXXFLAGS) -I$(INC_PATH) $< -o $@

$(OBJ_PATH)/Laplacian.o : Laplacian.cpp
		  $(MPI) -c $(MPIFLAGS) -I$(INC_PATH) $< -o $@

$(OBJ_PATH)/%.o : $(MPI_PATH)/%.cpp
		  $(MPI) -c $(MPIFLAGS) -I$(INC_PATH) $< -o $@

$(OBJ_PATH)/%.o : $(GPU_PATH)/%.cu
		  $(GPU) -c $(GPUFLAGS) -I$(INC_PATH) $< -o $@

$(OBJ_PATH)/%.o : $(UGT_PATH)/%.cpp
		  $(CXX) -c $(CXXFLAGS) $(UGTFLAGS) -I$(INC_PATH) $< -o $@

# Precompiled headers
%.hpp.gch	: %.hpp
		  $(CXX) $(LANG) -x c++-header $< -o $(GCH_PATH)/$@

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
obj/bench.o		: bench.cpp AOS.hpp CSC.hpp COO.hpp Vector.hpp
obj/csrbench.o		: csrbench.cpp CSR.hpp Vector.hpp
obj/densebench.o	: densebench.cpp Matrix.hpp Vector.hpp
obj/main.o		: main.cpp Matrix.hpp
obj/sparsebench.o	: sparsebench.cpp COO.hpp Vector.hpp
obj/mpi2norm_timer.o	: mpi/mpi2norm_timer.cpp Vector.hpp
obj/mpi2norm_driver.o	: mpi/mpi2norm_driver.cpp Vector.hpp
obj/mpiHeatEq2D-CG.o	: mpi/mpiHeatEq2D-CG.cpp Grid.hpp Laplacian.hpp
obj/mpiHeatEq2D-IR.o	: mpi/mpiHeatEq2D-IR.cpp Laplacian.hpp
obj/mpi_hello_world.o	: mpi/mpi_hello_world.cpp
obj/gpu_densebench.o	: gpu/gpu_densebench.cu
obj/unittests.o		: unittests/unittests.cpp Matrix.hpp Vector.hpp
obj/AOS.o		: AOS.hpp Vector.hpp
obj/COO.o		: COO.hpp Vector.hpp
obj/CSC.o		: CSC.hpp Vector.hpp
obj/CSR.o		: CSR.hpp Vector.hpp
obj/Matrix.o		: Matrix.hpp Vector.hpp
obj/Vector.o		: Vector.hpp
obj/Laplacian.o		: Laplacian.hpp
obj/Grid.o		: Grid.hpp
