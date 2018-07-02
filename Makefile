#
# @author Johnny Sellers
# 06/07/2017
# Makefile for HPC-Library.
#
LANG		= -std=c++11 # -Wc++11-extensions

CXX		= c++
MPI 	= mpic++ # --enable-mpi-cxx # -cxx=clang++

OPTS		= -Ofast -march=native -DNDEBUG
PICKY		= -g -Wall # -Wextra -pedantic
CPPFLAGS	+= -MD -MP

CXXFLAGS	= $(LANG) $(OPTS)
MPIFLAGS	= $(LANG) $(PICKY)

INC_PATH	= inc
LIB_PATH	= lib
OBJ_PATH	= obj
EXE_PATH	= exe

INCLUDES	= -I ./$(INC_PATH)

HEADERS		= $(wildcard inc/*.hpp)
TESTS		= $(wildcard tests/*.cpp)
MPI_TESTS = $(wildcard tests/mpi*.cpp)
SOURCES		= $(wildcard src/*.cpp)
# $(info $$var is [${MPI_TESTS}])

VPATH 	:= src tests
CLASSES		= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(SOURCES)))
BENCH 	= $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(notdir $(TESTS) $(MPI_TESTS)))
BENCH		= $(TESTS:.cpp=)
PCHS		= $(HEADERS:=.gch)
EXES		= $(OBJECTS:.o=)

all		: $(EXES)

classes 	: $(CLASSES)

bench	: $(OBJ_PATH)/AOS.o $(OBJ_PATH)/COO.o $(OBJ_PATH)/CSC.o
			$(CXX) $(CXXFLAGS) $^ -o $(EXE_PATH)/$@

csrbench	: $(OBJ_PATH)/CSR.o
			$(CXX) $(CXXFLAGS) $^ -o $(EXE_PATH)/$@

densebench	: $(OBJ_PATH)/Matrix.o
			$(CXX) $(CXXFLAGS) $^ -o $(EXE_PATH)/$@

sparsebench	: $(OBJ_PATH)/COO.o
			$(CXX) $(CXXFLAGS) $^ -o $(EXE_PATH)/$@

$(OBJ_PATH)/%.o : $(SOURCES) $(INC_PATH)/libhpc.h
		  $(CXX) -c $(CXXFLAGS) $< $(INCLUDES) -o $@

$(OBJ_PATH)/%.o 	: $(MPI_TESTS) $(INC_PATH)/libhpc.h
		  $(MPI) -c $(MPIFLAGS) $< $(INCLUDES) -o $@

.PHONY:	clean

clean:
	/bin/rm -f obj/* exe/*
