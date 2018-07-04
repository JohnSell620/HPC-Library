#
# @author Johnny Sellers
# 06/07/2017
# Makefile for HPC-Library.
#
LANG		= -std=c++11 # -Wc++11-extensions

CXX		= c++
MPI 		= mpic++ # -cxx=clang++ # --enable-mpi-cxx

OPTS		= -Ofast -march=native -DNDEBUG
CPPFLAGS	+= -MD -MP
PICKY		= -g -Wall # -Wextra -pedantic

CXXFLAGS	= $(LANG) $(OPTS)
MPIFLAGS	= $(LANG) # $(PICKY)

INC_PATH	= inc
LIB_PATH	= lib
OBJ_PATH	= obj
MPI_PATH 	= tests/mpi
EXE_PATH	= exe
GCH_PATH	= gch

INCLUDES	= -I ./$(INC_PATH)

VPATH 		:= src tests obj inc
HEADERS		= $(wildcard inc/*.hpp)
TESTS		= $(shell find tests/ ! -name "mpi*.cpp" -name "*.cpp")
MPI_TESTS 	= $(wildcard tests/mpi/*.cpp)
SOURCES		= $(wildcard src/*.cpp)

# $(info $$TESTS is [${TESTS}])

CLASSES		= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(SOURCES)))
BENCH 		= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(TESTS)))
MPI_BENCH 	= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(MPI_TESTS)))

OBJECTS		= $(CLASSES) $(BENCH) $(MPI_BENCH)
TARGETS		= $(notdir $(BENCH))
MPI_TARGETS	= $(notdir $(MPI_BENCH))
PCHS		= $(notdir $(HEADERS:=.gch))
EXECS		= $(TARGETS:.o=)
MPI_EXECS	= $(MPI_TARGETS:.o=)

all		: $(OBJECTS) $(EXECS) $(MPI_EXECS)

classes 	: $(CLASSES)

precomp_headers : $(PCHS)

$(EXECS)	:
			$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $(EXE_PATH)/$@

$(MPI_EXECS)	:
			$(MPI) $(MPIFLAGS) $(INCLUDES) $^ -o $(EXE_PATH)/$@

$(OBJ_PATH)/%.o : %.cpp
			$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

$(OBJ_PATH)/%.o : $(MPI_PATH)/%.cpp
			$(MPI) -c $(MPIFLAGS) $(INCLUDES) $< -o $@

%.hpp.gch	: %.hpp
			$(CXX) $(LANG) -x c++-header $< -o $(GCH_PATH)/$@

.PHONY:	clean

clean:
	/bin/rm -f obj/* exe/* gch/*

main			: main.o Matrix.o Vector.o
bench			: bench.o AOS.o COO.o CSC.o Vector.o
csrbench		: csrbench.o CSR.o Vector.o
densebench		: densebench.o Matrix.o Vector.o
sparsebench		: sparsebench.o COO.o Vector.o
mpi2norm_driver		: mpi2norm_driver.o Vector.o
mpi2norm_timer		: mpi2norm_timer.o Vector.o
bench.o			: bench.cpp AOS.hpp CSC.hpp COO.hpp Vector.hpp
csrbench.o		: csrbench.cpp CSR.hpp Vector.hpp
densebench.o		: densebench.cpp Matrix.hpp
sparsebench.o		: sparsebench.cpp COO.hpp Vector.hpp
mpi2norm_timer.o	: mpi/mpi2norm_timer.cpp Vector.hpp
mpi2norm_driver.o	: mpi/mpi2norm_driver.cpp Vector.hpp
AOS.o			: AOS.hpp Vector.hpp
COO.o			: COO.hpp Vector.hpp
CSC.o			: CSC.hpp Vector.hpp
CSR.o			: CSR.hpp Vector.hpp
Matrix.o		: Matrix.hpp Vector.hpp
Vector.o		: Vector.hpp
main.o			: main.cpp Matrix.hpp Vector.hpp
