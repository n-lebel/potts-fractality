CC = $(CXX)
CXXFLAGS = -std=c++17 -Wall -lcairo 
OBJDIR = obj

all: main

main: utilities.o hull1.o CMGrid.o key.o bond.o
#main: hull1.o CMGrid.o key.o bond.o

utilities.o: utilities.cpp hull1.h CMGrid.h
hull1.o: hull1.cpp hull1.h CMGrid.h
CMGrid.o: CMGrid.cpp bond.h CMGrid.h
bond.o: bond.cpp bond.h key.h
key.o: key.cpp key.h

clean:
	  $(RM) *.o $(objects) main
