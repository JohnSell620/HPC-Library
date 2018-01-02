#
# @author Johnny Sellers
# 06/07/2017
# Makefile for...
#
LANG			= -std=c++11 #-Wc++11-extensions
PICKY			= -g -Wall -Wextra -pedantic

CXX				= c++ #mpic++ -cxx=clang++
OPTS			= -Ofast -march=native -DNDEBUG
CPPFLAGS	+= -MD -MP

CXXFLAGS	= $(LANG) $(OPTS) #$(PICKY)

HEADERS		= Grid.hpp mpiStencil.hpp Final.hpp
FINAL 		= Final.cpp mpiStencil.cpp
TESTS			= ir-mpi.cpp cg-mpi.cpp
SOURCES		= Grid.cpp

OBJECTS		= $(SOURCES:.cpp=.o) $(TESTS:.cpp=.o) $(FINAL:.cpp=.o)
TARGETS		= $(TESTS:.cpp=)
PCHS			= $(HEADERS:=.gch)

all				: $(TARGETS)

%.o				: %.cpp
		  $(CXX) -c $(CXXFLAGS) $< -o $@

ir-mpi		: ir-mpi.o Grid.o mpiStencil.o Final.o
		  $(CXX) $(CXXFLAGS) $^ -o $@

cg-mpi		: cg-mpi.o Grid.o mpiStencil.o Final.o
		  $(CXX) $(CXXFLAGS) $^ -o $@

main			: main.o Matrix.o Vector.o
		  $(CXX) $(CXXFLAGS) $^ -o $@

hpcalc		: hpcalc.o Matrix.o Vector.o
		  $(CXX) $(CXXFLAGS) $^ -o $@

clean:
	/bin/rm -f ir-mpi ir-mpi.o cg-mpi cg-mpi.o Final.o mpiStencil.o Grid.o Matrix.o Vector.o main.o hpcalc.o main *.txt

# Grid.o: Grid.hpp
# mpiStencil.o: mpiStencil.hpp
# Final.o: Final.hpp
# ir-mpi.o: Grid.hpp mpiStencil.hpp Final.hpp
# cg-mpi.o: Grid.hpp mpiStencil.hpp Final.hpp
Matrix.o: Matrix.hpp Vector.hpp
Vector.o: Vector.hpp
main.o: main.cpp
hpcalc: hpcalc.cpp
